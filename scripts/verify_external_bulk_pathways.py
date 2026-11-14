"""Independent source-level check of the two proliferative whole-tissue studies.

Only cohort-level quantities leave this process. Ranked scores and clinical
labels stay private. Arithmetic uses sorted midranks and exact fractions.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from decimal import Decimal
from fractions import Fraction
import gzip
import hashlib
from itertools import combinations
import json
from pathlib import Path
import re
from statistics import mean, median
from math import comb, prod
import struct
import zlib

import numpy as np

from swiss_epithelial_mapping import LOCK_FILE, LOCK_SHA256

PINS = {
    "GPL570.annot.gz": "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394",
    "series.txt.gz": "b96db29937984402dc4fd13efdf6c609d1b58c007d11336a533fef1c86ccf7da",
    "GSE153739_family.soft.gz": "d26e97f9c15e62a8141874daca4cb33305a452152d1dbb59b845be49052b49f1",
    "GSE153739_GT_SO_5238_Transcript_Expression_Matrix.txt.gz": "627634eb5df06434cec03b12f177b7c90e4ea54438eb468b0fc50c0df20a9f4d",
    "wgEncodeGencodeAttrsV19.txt.gz": "d464efa58f6b72f0fcfea646addfe121c9b79a93cc1c02148562a841e492c941",
    "hgnc_complete_set_2026-10-06.tsv": "8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d",
    "locked_hgnc_projection.tsv": "6609ef4bb35608c7f791391b54344e0b0e9f77a44ed3004d95925aed2779a0c2",
    LOCK_FILE.name: LOCK_SHA256,
}


def pinned(path):
    with path.open("rb") as fp:
        actual = hashlib.file_digest(fp, "sha256").hexdigest()
    if actual != PINS[path.name]:
        raise ValueError("Unexpected source checksum: " + path.name)
    return path


def geo_annotation(path):
    with gzip.open(pinned(path), "rt") as fp:
        for line in fp:
            if line.rstrip() == "!platform_table_begin":
                break
        frame = csv.DictReader(fp, delimiter="\t")
        mapping = {}
        for r in frame:
            probe = r["ID"]
            if probe == "!platform_table_end":
                break
            if probe in mapping:
                raise ValueError("Repeated physical probe")
            mapping[probe] = 0
            if (not probe.startswith("AFFX") and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", r["Gene symbol"])
                    and re.fullmatch(r"[1-9][0-9]*", r["Gene ID"])):
                mapping[probe] = int(r["Gene ID"])
    if len(mapping) != 54675:
        raise ValueError("Unexpected physical probe denominator")
    return mapping


def us_data(source, annotation):
    with gzip.open(pinned(source), "rt") as fp:
        tab = csv.reader(fp, delimiter="\t")
        status = samples = None
        for r in tab:
            if r and r[0] == "!Sample_characteristics_ch1":
                status = r[1:]
            if r and r[0] == "!Sample_geo_accession":
                samples = r[1:]
            if r == ["!series_matrix_table_begin"]:
                break
        if not samples or not status or len(status) != len(samples) or next(tab) != ["ID_REF", *samples]:
            raise ValueError("Incomplete array sample metadata")
        if len(samples) != 37 or len(set(samples)) != 37:
            raise ValueError("Unexpected sample header")
        use = [j for j, s in enumerate(status, 1) if s in ("Proliferative Phase Normal", "Proliferative Phase Endometriosis")]
        case = [status[j-1] == "Proliferative Phase Endometriosis" for j in use]
        if (len(case), sum(case)) != (11, 6):
            raise ValueError("Unexpected proliferative donor design")
        ids, rows = [], []
        for r in tab:
            if r == ["!series_matrix_table_end"]:
                break
            ids.append(r[0])
            rows.append([float(r[i]) for i in use])
    if set(ids) != set(annotation) or len(ids) != len(annotation):
        raise ValueError("Expression and platform axes disagree")
    matrix = np.asarray(rows)
    if not np.isfinite(matrix).all():
        raise ValueError("Incomplete physical rank universe")
    scan_days = defaultdict(list)
    for j, diagnosis in zip(use, case):
        part = (source.parent / (samples[j-1] + ".header.gzpart")).read_bytes()
        header = zlib.decompressobj(31).decompress(part, 8192)
        if len(part) != 2048 or struct.unpack_from("<ii", header) != (64, 4):
            raise ValueError("Unexpected partial CEL format")
        size, = struct.unpack_from("<I", header, 20)
        if len(header) < size+24:
            raise ValueError("Truncated CEL header")
        dat = [line for line in header[24:24+size].splitlines() if line.startswith(b"DatHeader")]
        if len(dat) != 1:
            raise ValueError("Expected one native DAT header")
        days = re.findall(rb"\b[0-9]{2}/[0-9]{2}/[0-9]{2}(?= [0-9]{2}:[0-9]{2}:[0-9]{2})", dat[0])
        if len(days) != 1:
            raise ValueError("Unexpected native scan date")
        scan_days[days[0]].append(diagnosis)
    count = prod(comb(len(group), sum(group)) for group in scan_days.values())
    return np.array([annotation[p] for p in ids]), matrix, case, count


def transcript_mapping(source, hgnc):
    by_tx, by_ensembl, by_hgnc = defaultdict(set), defaultdict(set), defaultdict(set)
    with gzip.open(pinned(source), "rt") as fp:
        for r in csv.reader(fp, delimiter="\t"):
            gene, transcript = r[0].partition(".")[0], r[4].partition(".")[0]
            if re.fullmatch(r"ENST\d{11}", transcript) and re.fullmatch(r"ENSG\d{11}", gene):
                by_tx[transcript].add(gene)
    with pinned(hgnc).open() as fp:
        for r in csv.DictReader(fp, delimiter="\t"):
            if r["status"] == "Approved" and re.fullmatch(r"[1-9]\d*", r["entrez_id"]):
                entrez = int(r["entrez_id"])
                by_hgnc[entrez].add(r["hgnc_id"])
                if re.fullmatch(r"ENSG\d{11}", r["ensembl_gene_id"]):
                    by_ensembl[r["ensembl_gene_id"]].add(entrez)
    return by_tx, by_ensembl, by_hgnc


def mumbai_labels(path):
    labels, fields = {}, {}
    with gzip.open(pinned(path), "rt") as fp:
        for line in fp:
            key, sep, value = line.rstrip().partition(" = ")
            if key.startswith("^"):
                fields = {}
            elif key == "!Sample_description":
                fields["column"] = value
            elif key == "!Sample_characteristics_ch1":
                name, _, val = value.partition(": ")
                fields[name] = val
            if sep and set(fields) >= {"column", "phase of menstrual cycle", "group"}:
                if fields["phase of menstrual cycle"] != "Mid- proliferative" or fields["group"] not in {"Endometriosis", "Control"}:
                    raise ValueError("Unexpected biopsy population")
                labels[fields["column"]] = fields["group"] == "Endometriosis"
    if len(labels) != 7 or sum(labels.values()) != 4:
        raise ValueError("Unexpected RNA-seq donor design")
    return labels


def mumbai_data(directory, hgnc):
    labels = mumbai_labels(directory / "GSE153739_family.soft.gz")
    by_tx, by_ensembl, by_hgnc = transcript_mapping(directory / "wgEncodeGencodeAttrsV19.txt.gz", hgnc)
    source = directory / "GSE153739_GT_SO_5238_Transcript_Expression_Matrix.txt.gz"
    gene_rows, observed_ids, gene_sources = defaultdict(list), set(), defaultdict(set)
    values = []
    with gzip.open(pinned(source), "rt") as fp:
        records = csv.reader(fp, delimiter="\t")
        field, *names = next(records)
        if field != "Transcript_ID" or set(names) != labels.keys() or len(names) != len(labels):
            raise ValueError("Transcript header and clinical records disagree")
        for r in records:
            if not r:
                continue
            j = len(values)
            tx = r[0]
            if tx in observed_ids or len(r) != len(names)+1:
                raise ValueError("Duplicate transcript or unexpected record size")
            observed_ids.add(tx)
            values.append([Decimal(v) for v in r[1:]])
            identities = by_tx.get(tx, ())
            if len(identities) == 1:
                gene, = identities
                matches = by_ensembl.get(gene, ())
                if len(matches) == 1:
                    entrez, = matches
                    if len(by_hgnc[entrez]) == 1:
                        gene_rows[entrez].append(j)
                        gene_sources[entrez].add(gene)
    gene_rows = {g: rs for g, rs in gene_rows.items() if len(gene_sources[g]) == 1}
    if len(values) != 145630 or len(gene_rows) != 24977:
        raise ValueError("Unexpected transcript or mapped-gene denominator")
    return gene_rows, np.asarray(values, dtype=object), [labels[n] for n in names]


def ranks(x):
    order = np.argsort(x, kind="stable")
    ordered = x[order]
    starts = np.r_[0, np.flatnonzero(ordered[1:] != ordered[:-1])+1]
    ends = np.r_[starts[1:], len(order)]
    numerators = np.empty(len(order), dtype=np.int64)
    numerators[order] = np.repeat(starts+ends, ends-starts)
    return numerators  # percentile = 50 × numerator / len(x)


def cohort_gene_scores(matrix, measurement, common_ids, *, transcripts):
    """Map all reference ranks to target genes. Decimal sums preserve printed FPKM ties."""
    primary, alternate = [], []
    gene_ids = list(measurement)
    for values in matrix.T:
        row_rank = ranks(values)
        num = len(values)
        average = {g: Fraction(50*int(row_rank[measurement[g]].sum()), num*len(measurement[g]))
                   for g in common_ids}
        if transcripts:
            sums = np.array([sum(values[j] for j in measurement[g]) for g in gene_ids], dtype=object)
            gene_rank = ranks(sums)
            score = {g: Fraction(50*int(r), len(sums)) for g, r in zip(gene_ids, gene_rank) if g in common_ids}
            primary.append(score)
            alternate.append(average)
        else:
            primary.append(average)
            alternate.append({g: Fraction(50, num)*Fraction(median([int(x) for x in row_rank[measurement[g]]]))
                              for g in common_ids})
    return primary, alternate


def pathway_scores(rows, sets):
    scores = {}
    for key, ids in sets.items():
        scores[key] = [mean(s[g] for g in ids) for s in rows]
    scores["GLOBAL"] = [mean(row) for row in zip(*scores.values())]
    return scores


def contrast(score, assigned):
    case = [x for x, use in zip(score, assigned) if use]
    control = [x for x, use in zip(score, assigned) if not use]
    return mean(case) - mean(control)


def compatibility_set(scores, labels, tolerance=Fraction(1, 10**10)):
    """Invert all affine absolute tail inequalities at their rational intersections.

    Every assignment can change tail status only at one of the four signed
    intersections below. Extra intersections from inactive signs are harmless.
    """
    obs = contrast(scores, labels)
    mask = [Fraction(int(x)) for x in labels]
    threshold = Fraction(1, 20)
    functions = []
    boundaries = set()
    for chosen in combinations(range(len(labels)), sum(labels)):
        use = [i in chosen for i in range(len(labels))]
        a, b = contrast(scores, use), contrast(mask, use)
        functions.append((a, b))
        for p in (-1, 1):
            for q in (-1, 1):
                if p*b != q:
                    boundaries.add((p*a + tolerance - q*obs) / (p*b - q))
    points = sorted(boundaries)

    def accepted(delta):
        obs_adjusted = abs(obs - delta)
        tail = sum(abs(a - b*delta) + tolerance >= obs_adjusted for a, b in functions)
        return Fraction(tail, len(functions)) > threshold

    # The end cells are unbounded and every observed labeling is included.
    probes = [points[0]-1, *(mean((a, b)) for a, b in zip(points, points[1:])), points[-1]+1]
    finite = [p for p in points if accepted(p)]
    if not finite:
        raise ValueError("Compatibility set unexpectedly contains no boundary")
    lo = None if accepted(probes[0]) else min(finite)
    hi = None if accepted(probes[-1]) else max(finite)
    if any(accepted(p) != ((lo is None or lo <= p) and (hi is None or p <= hi)) for p in points+probes):
        raise ValueError("Compatibility set is disconnected")
    check_points = {"zero": Fraction(0), "lower": lo, "upper": hi}
    exact_tails = {name: sum(abs(a-b*p)+tolerance >= abs(obs-p) for a, b in functions)
                   for name, p in check_points.items() if p is not None}
    return {"bounds_pp": [float(x) if x is not None else None for x in (lo, hi)],
            "minimum_extreme_for_inclusion": int(len(functions)*threshold) + 1,
            "assignments": len(functions), "tail_checks": exact_tails}


def summarize(primary, alt, swiss, labels, exact):
    n = len(labels)
    perms = set(combinations(range(n), sum(labels))) if exact else ()
    result = {}
    for key, x in primary.items():
        obs = contrast(x, labels)
        deletes = [contrast(x[:i]+x[i+1:], labels[:i]+labels[i+1:]) for i in range(n)]
        row = dict(difference=float(obs), loo_low=float(min(deletes)), loo_high=float(max(deletes)),
                   alternate=float(contrast(alt[key], labels)), swiss_members=float(contrast(swiss[key], labels)))
        if exact:
            ex = sum(abs(contrast(x, [k in ix for k in range(n)]))+Fraction(1, 10**10) >= abs(obs)
                     for ix in perms)
            row.update(extreme=ex, p=ex/len(perms))
        result[key] = row
    if exact:
        paths = sorted((k for k in result if k != "GLOBAL"), key=lambda k: result[k]["p"])
        previous = 0
        for j, k in enumerate(paths):
            previous = max(previous, (len(paths)-j)*result[k]["extreme"])
            result[k]["holm"] = min(1, previous/len(perms))
        result["GLOBAL"]["compatibility"] = compatibility_set(primary["GLOBAL"], labels)
    return result


def verify(data, annotation, lock):
    nomination = json.loads(pinned(lock).read_text())["candidates"]
    original = {c["id"]: set(c["entrez_ids"]) for c in nomination if c["family"] == "pathway"}
    us_ids, us_matrix, us_case, us_allocations = us_data(data/"gse6364/series.txt.gz", geo_annotation(annotation))
    m_rows, m_matrix, m_case = mumbai_data(data/"gse153739", data/"gse266265/hgnc_complete_set_2026-10-06.tsv")
    us_rows = defaultdict(list)
    for j, g in enumerate(us_ids):
        if g:
            us_rows[int(g)].append(j)
    common = {k: g & us_rows.keys() & m_rows.keys() for k, g in original.items()}
    if len(common) != 5 or any(len(g) < 15 or len(g)*5 < len(original[k])*4 for k, g in common.items()):
        raise ValueError("One or more pathways outside the coverage rule")
    swiss_file = LOCK_FILE.parents[1] / "external-bulk-2026-11-14/locked_hgnc_projection.tsv"
    with pinned(swiss_file).open() as fp:
        swiss_ids = {int(r["entrez_id"]) for r in csv.DictReader(fp, delimiter="\t") if r["feature_label"]}
    swiss_sets = {k: v & swiss_ids for k, v in common.items()}
    out = {"us_scan_date_assignments": us_allocations,
           "membership": {key: {"primary": len(members), "swiss": len(swiss_sets[key]), "locked": len(original[key]),
                                "entrez_ids": sorted(members), "swiss_entrez_ids": sorted(swiss_sets[key])}
                          for key, members in common.items()}}
    for name, measurement, values, labels, tx in (
            ("GSE6364", us_rows, us_matrix, us_case, False),
            ("GSE153739", m_rows, m_matrix, m_case, True)):
        base, alternative = cohort_gene_scores(values, measurement, set.union(*common.values()), transcripts=tx)
        out[name] = summarize(pathway_scores(base, common), pathway_scores(alternative, common),
                              pathway_scores(base, swiss_sets), labels, exact=tx)
    return out


def compare(observed, primary):
    """Fail on disagreements in projections, estimates, sensitivities or tests."""
    if observed["us_scan_date_assignments"] != primary["us_scan_date_assignments"]:
        raise ValueError("Scan-conditional design size disagrees")
    if set(observed["membership"]) != set(primary["membership"]):
        raise ValueError("Different nominee IDs")
    for path, counts in observed["membership"].items():
        if any(primary["membership"][path][k] != v for k, v in counts.items()):
            raise ValueError("Different measured gene membership for " + path)
    errors = []
    for cohort in ("GSE6364", "GSE153739"):
        if set(primary[cohort]) != set(observed[cohort]):
            raise ValueError("Unexpected endpoint coverage")
        for key, row in observed[cohort].items():
            original = primary[cohort][key]
            if original.keys() != row.keys() - {"compatibility"}:
                raise ValueError("Unexpected primary result fields")
            for field in original:
                if field in {"extreme", "p", "holm"} and original[field] != row[field]:
                    raise ValueError("Exact randomization or adjusted p disagrees")
                errors.append(abs(row[field] - original[field]))
    independently_inverted = observed["GSE153739"]["GLOBAL"]["compatibility"]
    interval = primary["mumbai_global_interval"]
    for field in ("tail_checks", "assignments", "minimum_extreme_for_inclusion"):
        if interval[field] != independently_inverted[field]:
            raise ValueError("Interval boundary, assignment count or tail threshold disagrees")
    for bound, independent in zip(interval["bounds_pp"], independently_inverted["bounds_pp"], strict=True):
        if (bound is None) != (independent is None):
            raise ValueError("Interval boundedness differs")
        if bound is not None:
            errors.append(abs(bound - independent))
    if not all(error <= 1e-10 for error in errors):
        raise ValueError("Different primary, sensitivity or compatibility estimates")
    mumbai = observed["GSE153739"]
    gate = mumbai["GLOBAL"]["p"] <= .05
    supported = sorted(k for k, v in mumbai.items() if k != "GLOBAL" and gate and v["holm"] <= .05)
    if gate != primary["mumbai_global_gate_passed"] or supported != sorted(primary["mumbai_supported_components"]):
        raise ValueError("Prespecified inferential decision disagrees")
    return {"distinct_endpoints": 12, "exact_tails_verified": 6, "holm_values_verified": 5,
            "all_membership_ids_verified": True, "us_scan_date_assignments": observed["us_scan_date_assignments"],
            "mumbai_interval_and_global_decision_verified": True,
            "maximum_absolute_difference": max(errors)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/private"))
    parser.add_argument("--annotation", type=Path, default=Path("data/raw/GPL570.annot.gz"))
    parser.add_argument("--lock", type=Path, default=LOCK_FILE)
    parser.add_argument("--compare", type=Path, help="Report agreement with the primary aggregate JSON")
    options = parser.parse_args()
    results = verify(options.data, options.annotation, options.lock)
    if options.compare:
        results = compare(results, json.loads(options.compare.read_text()))
    print(json.dumps(results, indent=2, allow_nan=False))
