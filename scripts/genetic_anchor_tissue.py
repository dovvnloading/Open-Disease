"""Paired lesion contrasts and secondary uterine checks for five genetically mapped genes."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from itertools import product
from pathlib import Path

import numpy as np
from scipy.stats import binomtest, rankdata

from analyze_melbourne_genes import SOURCE_HASH, candidate_probes
from bulk_gse153739 import load_mumbai
from cross_cohort_discovery import assignment_design, holm
from genome_rome import load_rome, read_annotation
from rome_scan import audit as audit_rome_scan
from open_disease.geo_microarray import sha256_file
from swiss_epithelial_mapping import HGNC_FILE

ROOT = Path(__file__).resolve().parents[1]
GENES = {
    "LNC-LBCS": 100506885,
    "TRA2A": 29896,
    "VEZT": 55591,
    "SRP14": 6727,
    "HOXB9": 3219,
}
SERIES = {"GSE7305": ("GPL570", 10, 54675), "GSE11691": ("GPL96", 9, 22283)}
PINS = {
    "GSE7305_series_matrix.txt.gz": "7dab7f8918ab0ac09f98ea2cba2f5baa378b0444537e957061a914306768cc90",
    "GSE11691_series_matrix.txt.gz": "7a2fea992101d6a04360d7becc954b2cbfc62400854fe533072776f7039d45f6",
    "GPL570.annot.gz": "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394",
    "GPL96.annot.gz": "88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4",
    "discovery_gene.tsv": "9df65b29615ec56767e9abaaf3844f8ebde5e21baf8cdab193ac24857a170da5",
}


def checked(path):
    if sha256_file(path) != PINS[path.name]:
        raise ValueError(f"Unexpected input digest for {path.name}")
    return path


def one(header, name):
    (row,) = header[name]
    return row


def pair_columns(series, header):
    """Resolve only the identities that appear explicitly in the deposited fields."""
    title = one(header, "!Sample_title")
    source = one(header, "!Sample_source_name_ch1")
    desc = one(header, "!Sample_description")
    pairs = defaultdict(dict)
    for i, (name, origin, description) in enumerate(
        zip(title, source, desc, strict=True)
    ):
        if series == "GSE7305":
            match = re.fullmatch(r"Female - caucasian - ([EN])-(\d+-[MG])", origin)
            if not match:
                raise ValueError("Unknown source-key format")
            letter, key = match.groups()
            tissue = "lesion" if letter == "E" else "uterus"
            prefix = (
                "Endometrium/Ovary-Disease" if letter == "E" else "Endometrium-Normal"
            )
            if not re.fullmatch(prefix + r" \d+", name):
                raise ValueError("Source tissue and title disagree")
        elif series == "GSE11691":
            match = re.fullmatch(r"(Endometriosis|Endometrium)-(\d+)", name)
            if not match:
                raise ValueError("Unknown source title")
            tissue_name, key = match.groups()
            tissue = "lesion" if tissue_name == "Endometriosis" else "uterus"
            suffix = {
                "eutopic endometrium": "U",
                "endometriosis, peritoneal": "P",
                "endometriosis, broard ligament": "BL",
            }.get(origin)
            if (
                not suffix
                or (suffix == "U") != (tissue == "uterus")
                or description != key + suffix
            ):
                raise ValueError("Tissue, donor title and sample description disagree")
        else:
            raise ValueError("Unrecognized series")
        if tissue in pairs[key]:
            raise ValueError("Duplicate tissue from the same donor")
        pairs[key][tissue] = i
    if len(pairs) != SERIES[series][1] or any(
        set(p) != {"lesion", "uterus"} for p in pairs.values()
    ):
        raise ValueError("Incomplete or unexpected pair set")
    lesion = [p["lesion"] for p in pairs.values()]
    uterus = [p["uterus"] for p in pairs.values()]
    if series == "GSE7305" and any(
        desc[i] != desc[j] for i, j in zip(lesion, uterus, strict=True)
    ):
        raise ValueError("Matched phases disagree")
    return lesion, uterus


def paired_summary(delta):
    delta = np.asarray(delta, dtype=float)
    if delta.ndim != 1 or not 2 <= len(delta) <= 16 or not np.isfinite(delta).all():
        raise ValueError("Expected 2 to 16 finite paired differences")
    signs = np.array(list(product([-1, 1], repeat=len(delta))))
    effect = float(delta.mean())
    null = signs @ delta / len(delta)
    up = int(sum(null >= effect - 1e-10))
    down = int(sum(null <= effect + 1e-10))
    two = int(sum(abs(null) >= abs(effect) - 1e-10))
    positive, negative = int(sum(delta > 1e-10)), int(sum(delta < -1e-10))
    sign_p = (
        binomtest(positive, positive + negative).pvalue if positive + negative else 1
    )
    omitted = (delta.sum() - delta) / (len(delta) - 1)
    return dict(
        effect_pp=effect,
        pairs=len(delta),
        higher=positive,
        lower=negative,
        ties=len(delta) - positive - negative,
        assignments=len(null),
        upper=up,
        lower_tail=down,
        two_sided=two,
        p_sign=float(sign_p),
        omit_low=float(omitted.min()),
        omit_high=float(omitted.max()),
    )


def read_paired_arrays(directory, series):
    platform, pairs, expected_features = SERIES[series]
    mapping = read_annotation(checked(directory / f"{platform}.annot.gz"))
    path = checked(directory / f"{series}_series_matrix.txt.gz")
    with gzip.open(path, "rt", newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        header = defaultdict(list)
        for row in rows:
            if row == ["!series_matrix_table_begin"]:
                break
            if row:
                header[row[0]].append(row[1:])
        cols = one(header, "!Sample_geo_accession")
        if (
            next(rows) != ["ID_REF", *cols]
            or len(set(cols)) != 2 * pairs
            or set(one(header, "!Sample_platform_id")) != {platform}
        ):
            raise ValueError("Matrix columns or platform disagree")
        lesion, uterus = pair_columns(series, header)
        probes, signals = [], []
        for row in rows:
            if row == ["!series_matrix_table_end"]:
                break
            if len(row) != 1 + len(cols):
                raise ValueError("Unexpected sample count in feature row")
            probes.append(row[0])
            signals.append([float(x) for x in row[1:]])
        else:
            raise ValueError("Incomplete matrix")
    if (
        len(probes) != len(set(probes))
        or set(probes) != mapping.keys()
        or len(probes) != expected_features
    ):
        raise ValueError("Incomplete platform or duplicate features")
    values = np.array(signals)
    if not np.isfinite(values).all():
        raise ValueError("Audited complete reference changed")
    return probes, values, lesion, uterus, mapping


def paired_arrays(directory, series):
    probes, values, lesion, uterus, mapping = read_paired_arrays(directory, series)
    ranks = 100 * (rankdata(values, axis=0) - 0.5) / len(values)
    mapped = np.array([mapping[p] for p in probes])
    out = {}
    for symbol, gene in GENES.items():
        selected = mapped == gene
        if not selected.any():
            out[symbol] = None
            continue
        diff = ranks[selected][:, lesion] - ranks[selected][:, uterus]
        result = paired_summary(diff.mean(axis=0))
        result["probes"] = {
            probes[j]: float(d)
            for j, d in zip(np.flatnonzero(selected), diff.mean(axis=1), strict=True)
        }
        # MAS5 is deposited as positive linear-scale signal. Hever's RMA
        # submission does not explicitly document the final log/linear units.
        if series == "GSE11691" and np.all(values[selected] > 0):
            log_values = np.log2(values[selected])
            result["mean_log2_signal_ratio"] = float(
                (log_values[:, lesion] - log_values[:, uterus]).mean()
            )
        out[symbol] = result
    return out


def common_direction(studies):
    if any(x is None for x in studies):
        return 1.0
    p_up = max(x["upper"] / x["assignments"] for x in studies)
    p_down = max(x["lower_tail"] / x["assignments"] for x in studies)
    return min(1.0, 2 * min(p_up, p_down))


def case_control(genes, scores, cases, strata):
    obs, permutations, _ = assignment_design(cases, strata)
    result = {}
    for symbol, gene in GENES.items():
        pos = np.flatnonzero(genes == gene)
        if not len(pos):
            result[symbol] = None
            continue
        values = scores[:, pos.item()]
        effect = float(obs @ values)
        tail = int(sum(abs(permutations @ values) >= abs(effect) - 1e-10))
        result[symbol] = dict(
            effect_pp=effect,
            p=tail / len(permutations),
            extreme=tail,
            assignments=len(permutations),
        )
    return result


def analyze(data):
    import openpyxl

    directory = data / "genetics-expression"
    pairs = {study: paired_arrays(directory, study) for study in SERIES}
    joint = {s: common_direction([study[s] for study in pairs.values()]) for s in GENES}
    adjusted = dict(zip(GENES, holm(list(joint.values())), strict=True))
    path = ROOT / "evidence/cross-cohort-2026-11-14/discovery_gene.tsv"
    with checked(path).open() as stream:
        fi = {
            int(r["id"]): dict(effect_pp=float(r["effect_pp"]), p=float(r["p_raw"]))
            for r in csv.DictReader(stream, delimiter="\t")
        }
    finland = {s: fi.get(i) for s, i in GENES.items()}
    genes, scores, cases, strata = load_rome(data / "gse25628")
    rome = case_control(genes, scores, cases, strata)
    days, receipt = audit_rome_scan(data / "gse25628")
    if (
        receipt["ordered_sha256"]
        != "360626a37e13529bb37331a8f389a9ebe0e597ba1da43cbd1500e0843d439825"
    ):
        raise ValueError("Native Rome scan-day source changed")
    days = np.asarray(days)
    keep = np.array(
        [np.any(cases[days == day]) and np.any(~cases[days == day]) for day in days]
    )
    scan_matched = case_control(genes, scores[keep], cases[keep], days[keep])
    scan_ps = [v["p"] if v else 1.0 for v in scan_matched.values()]
    scan_holm = dict(zip(GENES, holm(scan_ps), strict=True))
    ids, scores, cases, _ = load_mumbai(
        data / "gse153739", data / "gse266265" / HGNC_FILE
    )
    mumbai = case_control(ids, scores, cases, np.zeros(len(cases)))
    bulk = {}
    for s in GENES:
        a, b = finland[s], rome[s]
        bulk[s] = (
            max(a["p"], b["p"])
            if a and b and a["effect_pp"] * b["effect_pp"] > 0
            else 1.0
        )
    workbook = data / "fung-2018/Fung_supplementary_tables.xlsx"
    if sha256_file(workbook) != SOURCE_HASH:
        raise ValueError("Unexpected Melbourne workbook")
    book = openpyxl.load_workbook(workbook, read_only=True, data_only=True)
    try:
        rows = candidate_probes(
            iter(book["TableS20"].values),
            targets={s: (i, "") for s, i in GENES.items()},
        )
    finally:
        book.close()
    melbourne = {}
    for s in GENES:
        probes = {r["probe_id"]: r["two_sided_p"] for r in rows if r["symbol"] == s}
        melbourne[s] = dict(
            probes=probes,
            p_bonferroni=(
                min(1.0, len(probes) * min(probes.values())) if probes else None
            ),
        )
    return dict(
        paired=pairs,
        common_direction_p=joint,
        common_direction_holm=adjusted,
        finland=finland,
        rome=rome,
        rome_scan_matched=scan_matched,
        rome_scan_matched_holm=scan_holm,
        mumbai=mumbai,
        melbourne=melbourne,
        uterine_conjunction_p=bulk,
        uterine_conjunction_holm=dict(
            zip(GENES, holm(list(bulk.values())), strict=True)
        ),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze(args.data), indent=2, allow_nan=False))
