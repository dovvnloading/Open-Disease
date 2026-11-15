"""Independent offline check of all paired genetic-anchor results. No participant output.

Ranks use binary searches of sorted columns. Integer subset sums enumerate every null
configuration and preserve ties exactly. No functions from the analysis are imported.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

FAMILY = {
    "LNC-LBCS": "100506885",
    "TRA2A": "29896",
    "VEZT": "55591",
    "SRP14": "6727",
    "HOXB9": "3219",
}
STUDIES = {"GSE7305": ("GPL570", 10, 54675), "GSE11691": ("GPL96", 9, 22283)}
PINS = {
    "GPL570.annot.gz": "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394",
    "GPL96.annot.gz": "88e0b22362bac779eb220b3b185c80faa6510a92b9358eaad159a561ab4351c4",
    "GSE7305_series_matrix.txt.gz": "7dab7f8918ab0ac09f98ea2cba2f5baa378b0444537e957061a914306768cc90",
    "GSE11691_series_matrix.txt.gz": "7a2fea992101d6a04360d7becc954b2cbfc62400854fe533072776f7039d45f6",
}


def parse_geo(path, start, stop):
    with path.open("rb") as file:
        if hashlib.file_digest(file, "sha256").hexdigest() != PINS[path.name]:
            raise ValueError("Input does not match the source digest")
    metadata = defaultdict(list)
    with gzip.open(path, "rt", newline="") as source:
        lines = csv.reader(source, delimiter="\t")
        for row in lines:
            if row == [start]:
                break
            if row and row[0].startswith("!Sample_"):
                metadata[row[0]].append(row[1:])
        columns = next(lines)
        data = []
        for row in lines:
            if row == [stop]:
                return metadata, columns, data
            if len(row) != len(columns):
                raise ValueError("Unexpected GEO row width")
            data.append(row)
    raise ValueError("Missing GEO table end")


def donor_pairs(series, metadata):
    (titles,) = metadata["!Sample_title"]
    (origin,) = metadata["!Sample_source_name_ch1"]
    (notes,) = metadata["!Sample_description"]
    pairs = defaultdict(dict)
    for i, title in enumerate(titles):
        if series == "GSE7305":
            # The source field contains tissue code + shared pair key.
            label = origin[i].split(" - ")[-1]
            side, key = label.split("-", 1)
            normal = title.startswith("Endometrium-Normal ")
            if side != ("N" if normal else "E"):
                raise ValueError("Hever source and title conflict")
        else:
            # Independently use the description's anatomical suffix and numeric prefix.
            match = re.fullmatch(r"(\d+)(U|P|BL)", notes[i])
            if match is None:
                raise ValueError("Unrecognized Hull sample description")
            key, site = match.groups()
            normal = site == "U"
            expected = "Endometrium-" if normal else "Endometriosis-"
            if (
                title != expected + key
                or (origin[i] == "eutopic endometrium") != normal
            ):
                raise ValueError("Hull source, description and title conflict")
        if normal in pairs[key]:
            raise ValueError("Duplicate tissue in a pair")
        pairs[key][normal] = i
    if any(set(p) != {True, False} for p in pairs.values()):
        raise ValueError("Incomplete pair")
    both = [(p[False], p[True]) for p in pairs.values()]
    if series == "GSE7305" and any(notes[a] != notes[b] for a, b in both):
        raise ValueError("Paired Hever phase metadata disagree")
    return both


def exact_pair_test(delta, scale):
    # delta = integer sums of doubled empirical-CDF midpoint ranks.
    n = len(delta)
    total = sum(delta)
    subset_sums = [0]
    for d in delta:
        subset_sums += [s + d for s in subset_sums]
    random_totals = [total - 2 * s for s in subset_sums]
    positive = sum(d > 0 for d in delta)
    negative = sum(d < 0 for d in delta)
    m = positive + negative
    tail = sum(math.comb(m, j) for j in range(min(positive, negative) + 1))
    omitted = [scale * (total - d) / (n - 1) for d in delta]
    return dict(
        effect_pp=scale * total / n,
        pairs=n,
        higher=positive,
        lower=negative,
        ties=n - m,
        assignments=2**n,
        upper=sum(s >= total for s in random_totals),
        lower_tail=sum(s <= total for s in random_totals),
        two_sided=sum(abs(s) >= abs(total) for s in random_totals),
        p_sign=min(1, 2 * tail / 2**m),
        omit_low=min(omitted),
        omit_high=max(omitted),
    )


def calculate_study(directory, name):
    platform, count, features = STUDIES[name]
    _, col, annotation = parse_geo(
        directory / (platform + ".annot.gz"),
        "!platform_table_begin",
        "!platform_table_end",
    )
    namecol, symcol, idcol = (col.index(k) for k in ("ID", "Gene symbol", "Gene ID"))
    target_features = defaultdict(set)
    for row in annotation:
        if row[idcol] in FAMILY.values() and not row[namecol].startswith("AFFX"):
            if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", row[symcol]):
                target_features[row[idcol]].add(row[namecol])
    metadata, columns, data = parse_geo(
        directory / (name + "_series_matrix.txt.gz"),
        "!series_matrix_table_begin",
        "!series_matrix_table_end",
    )
    pairs = donor_pairs(name, metadata)
    probes = [r[0] for r in data]
    if len(set(probes)) != features or len(data) != features or len(pairs) != count:
        raise ValueError("Unexpected reference or donor count")
    (accession,) = metadata["!Sample_geo_accession"]
    if columns != ["ID_REF", *accession] or len(set(accession)) != 2 * count:
        raise ValueError("Matrix column identity mismatch")
    vals = np.asarray([row[1:] for row in data], dtype=np.float64)
    if not np.isfinite(vals).all():
        raise ValueError("Expected complete numeric reference")
    order = np.sort(vals, axis=0)
    by_probe = {p: i for i, p in enumerate(probes)}
    output = {}
    for symbol, gene_id in FAMILY.items():
        selected = sorted(target_features[gene_id])
        if not selected:
            output[symbol] = None
            continue
        probe_delta = {}
        ratios = []
        for probe in selected:
            target = vals[by_probe[probe]]
            midpoint = []
            for j, value in enumerate(target):
                column = order[:, j]
                midpoint.append(
                    int(
                        np.searchsorted(column, value, "left")
                        + np.searchsorted(column, value, "right")
                    )
                )
            probe_delta[probe] = [midpoint[a] - midpoint[b] for a, b in pairs]
            if name == "GSE11691":
                ratios.extend(math.log2(target[a] / target[b]) for a, b in pairs)
        delta = list(map(sum, zip(*probe_delta.values(), strict=True)))
        one_probe = 50 / features
        out = exact_pair_test(delta, one_probe / len(selected))
        out["probes"] = {p: sum(d) * one_probe / count for p, d in probe_delta.items()}
        if ratios:
            out["mean_log2_signal_ratio"] = sum(ratios) / len(ratios)
        output[symbol] = out
    return output


def compare(actual, reference):
    if isinstance(actual, dict):
        if set(actual) != set(reference):
            raise ValueError("Aggregate fields disagree")
        for key in actual:
            compare(actual[key], reference[key])
    elif actual is None:
        if reference is not None:
            raise ValueError("Coverage differs")
    elif not math.isclose(actual, reference, abs_tol=1e-10, rel_tol=1e-12):
        raise ValueError("Numerical mismatch")


def verify(directory, expected):
    data = json.loads(expected.read_text())
    results = {name: calculate_study(directory, name) for name in STUDIES}
    compare(results, data["paired"])
    p = {}
    for gene in FAMILY:
        if any(study[gene] is None for study in results.values()):
            p[gene] = Fraction(1)
            continue
        worst = [
            max(
                Fraction(study[gene][tail], study[gene]["assignments"])
                for study in results.values()
            )
            for tail in ("upper", "lower_tail")
        ]
        p[gene] = min(Fraction(1), 2 * min(worst))
    ordered = sorted(p, key=p.get)
    holm = {}
    for i, key in enumerate(ordered):
        holm[key] = min(1, max((len(p) - j) * p[ordered[j]] for j in range(i + 1)))
    compare(p, data["common_direction_p"])
    compare(holm, data["common_direction_holm"])
    return {
        "paired_study_gene_estimates_verified": sum(
            v is not None for s in results.values() for v in s.values()
        ),
        "absent_study_gene_comparisons_verified": sum(
            v is None for s in results.values() for v in s.values()
        ),
        "directional_conjunctions_verified": len(p),
        "holm_adjustments_verified": len(holm),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("aggregates", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.data, args.aggregates), indent=2))
