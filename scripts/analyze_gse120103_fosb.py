"""Run the prospectively specified GSE120103 FOSB percentile comparison.

Reads two pinned local GEO files and prints only aggregate results.
"""
from __future__ import annotations

import argparse
from array import array
from collections import defaultdict
import csv
import gzip
import json
import math
from pathlib import Path
import re

import numpy as np
from scipy.stats import t

from open_disease.geo_microarray import sha256_file


PROBE = "A_23_P429998"
MATRIX_HASH = "0ab4900ba460232b1c76ca25028b3a034dd50491d3a7e63d4bf2fa2c4afbc16e"
ANNOTATION_HASH = "aacec81b2d95f183f0827b43e769388a92917811b81920c13ecc8fe148c5a3ab"
MIN_COMMON_FEATURES = 20_000
ARRAYS = 36
# GEO labels the source for affected infertile women "Fertile" but gives code 2B.
# Check the entire submitted string, together with the title and characteristic.
GROUPS = {
    "fertile_control": ("Disease free Endometrium of fertile women",
                        "Fertile Secretory without Endometriosis (Group 1A)"),
    "fertile_case": ("Endometrium  from Stage IV Ovarian Endometriosis of fertile women",
                     "Fertile Secretory with Stage IV Ovarian Endometriosis (Group 2A)"),
    "infertile_control": ("Disease free Endometrium of Infertile women",
                          "Infertile Secretory without Endometriosis (Group 1B)"),
    "infertile_case": ("Endometrium  from Stage IV Ovarian Endometriosis of Infertile women",
                       "Fertile Secretory with Stage IV Ovarian Endometriosis (Group 2B)"),
}


def check_annotation(path):
    matches = []
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        for row in reader:
            if row == ["!platform_table_begin"]:
                break
        columns = next(reader, [])
        required = ("ID", "Gene symbol", "Gene ID")
        if any(columns.count(key) != 1 for key in required):
            raise ValueError("Missing or ambiguous probe annotation header")
        idx, symbol_col, gene_col = (columns.index(key) for key in required)
        for row in reader:
            if row == ["!platform_table_end"]:
                break
            if len(row) != len(columns):
                raise ValueError("Malformed platform annotation row")
            feature, symbol, gene = row[idx], row[symbol_col], row[gene_col]
            if (feature == PROBE or "FOSB" in re.split(r"\s*///\s*", symbol)
                    or "2354" in re.split(r"\s*///\s*", gene)):
                matches.append((feature, symbol, gene))
        else:
            raise ValueError("Missing platform table terminator")
        if any(any(row) for row in reader):
            raise ValueError("Unexpected trailing platform content")
    if matches != [(PROBE, "FOSB", "2354")]:
        raise ValueError("FOSB probe mapping differs from the frozen annotation")


def checked_groups(header):
    if (header["!Series_geo_accession"] != [["GSE120103"]]
            or header["!Series_platform_id"] != [["GPL6480"]]):
        raise ValueError("Unexpected GEO series or platform")

    def vector(key):
        rows = header[key]
        if len(rows) != 1 or len(rows[0]) != ARRAYS:
            raise ValueError("Missing, repeated or misaligned GEO sample field")
        return rows[0]

    ids = vector("!Sample_geo_accession")
    if len(set(ids)) != ARRAYS or not all(re.fullmatch(r"GSM[1-9][0-9]*", item) for item in ids):
        raise ValueError("Expected 36 distinct GEO array accessions")
    if vector("!Sample_platform_id") != ["GPL6480"] * ARRAYS:
        raise ValueError("Sample platform mismatch")
    if vector("!Sample_channel_count") != ["1"] * ARRAYS:
        raise ValueError("Expected only single-channel arrays")

    titles = vector("!Sample_title")
    sources = vector("!Sample_source_name_ch1")
    rows = header["!Sample_characteristics_ch1"]
    if not rows or any(len(row) != ARRAYS for row in rows):
        raise ValueError("Missing or misaligned source characteristics")
    group_by_label = {label: group for group, (label, _) in GROUPS.items()}
    replicates = defaultdict(list)
    groups = []
    for title, source, entries in zip(titles, sources, zip(*rows), strict=True):
        fields = {}
        for entry in entries:
            key, sep, value = entry.partition(": ")
            if not sep or key in fields:
                raise ValueError("Malformed or duplicated sample characteristic")
            fields[key] = value
        label = fields.get("sample group")
        group = group_by_label.get(label)
        if group is None or fields.get("tissue") != "Endometrium":
            raise ValueError("Unexpected clinical group or tissue")
        if source != "human endometrium, " + GROUPS[group][1]:
            raise ValueError("Secretory source name or group code disagrees with group")
        match = re.fullmatch(re.escape(label) + r", replicate ([1-9])", title)
        if match is None:
            raise ValueError("Sample title or replicate disagrees with group")
        replicates[group].append(int(match[1]))
        groups.append(group)
    if any(sorted(replicates[group]) != list(range(1, 10)) for group in GROUPS):
        raise ValueError("Expected replicates 1 through 9 once each in every group")
    return ids, groups


def read_percentiles(path):
    seen = set()
    complete = array("d")
    target = None
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        rows = csv.reader(stream, delimiter="\t")
        header = defaultdict(list)
        for row in rows:
            if row == ["!series_matrix_table_begin"]:
                break
            if row and row[0].startswith(("!Series_", "!Sample_")):
                header[row[0]].append(row[1:])
        else:
            raise ValueError("Missing start of matrix")
        ids, groups = checked_groups(header)
        if next(rows, None) != ["ID_REF", *ids]:
            raise ValueError("Matrix columns disagree with sample order")
        for row in rows:
            if row == ["!series_matrix_table_end"]:
                break
            if len(row) != ARRAYS + 1 or not row[0].strip() or row[0] in seen:
                raise ValueError("Invalid feature row or duplicate probe")
            seen.add(row[0])
            try:
                values = [math.nan if token.strip().casefold() in {"", "na", "null"}
                          else float(token) for token in row[1:]]
            except ValueError:
                raise ValueError("Malformed submitted intensity token") from None
            if all(map(math.isfinite, values)):
                complete.extend(values)
                if row[0] == PROBE:
                    target = np.array(values)
        else:
            raise ValueError("Missing end of matrix")
        if any(any(row) for row in rows):
            raise ValueError("Unexpected trailing matrix content")

    if target is None:
        raise ValueError("Frozen probe is missing or incomplete")
    feature_values = np.frombuffer(complete, dtype=np.float64).reshape(-1, ARRAYS)
    common_count = len(feature_values)
    if common_count < MIN_COMMON_FEATURES:
        raise ValueError("Insufficient complete features for the common reference")
    below = (feature_values < target).sum(axis=0)
    tied = (feature_values == target).sum(axis=0)
    percentiles = 100 * (below + tied / 2) / common_count
    cells = {group: percentiles[np.array(groups) == group] for group in GROUPS}
    return cells, len(seen), common_count


def summarize(cells):
    differences = {stratum: cells[stratum + "_case"].mean() - cells[stratum + "_control"].mean()
                   for stratum in ("fertile", "infertile")}
    estimate = float(np.mean(list(differences.values())))
    variances = {group: values.var(ddof=1) / (4 * len(values)) for group, values in cells.items()}
    total_variance = sum(variances.values())
    if not np.isfinite(total_variance) or total_variance <= 0:
        raise ValueError("The prespecified contrast variance cannot be estimated")
    df = total_variance**2 / sum(v**2 / (len(cells[g]) - 1) for g, v in variances.items())
    se = math.sqrt(total_variance)
    margin = t.ppf(0.975, df) * se
    low, high = estimate - margin, estimate + margin
    p_value = 2 * t.sf(abs(estimate) / se, df)
    loo = []
    for group, values in cells.items():
        signed_weight = 0.5 if group.endswith("_case") else -0.5
        loo.extend(estimate + signed_weight * (values.mean() - values) / (len(values) - 1))
    decision = "inconclusive"
    if low > 0:
        decision = "supports_expected_direction"
    elif high < 0:
        decision = "supports_opposite_direction"
    return {
        "stratum_differences_percentile_points": differences,
        "difference_percentile_points": estimate, "ci95_percentile_points": [low, high],
        "standard_error_percentile_points": se, "welch_df": df, "two_sided_p": p_value,
        "leave_one_out_range_percentile_points": [min(loo), max(loo)], "decision": decision,
    }


def analyze(matrix, annotation):
    if sha256_file(matrix) != MATRIX_HASH or sha256_file(annotation) != ANNOTATION_HASH:
        raise ValueError("Input SHA-256 differs from the audited snapshot")
    check_annotation(annotation)
    cells, total, common = read_percentiles(matrix)
    return {
        "series": "GSE120103", "platform": "GPL6480", "probe": PROBE,
        "counts": {group: len(values) for group, values in cells.items()},
        "features_in_matrix": total, "complete_features": common,
        "excluded_features": total - common, "features_tested": 1, **summarize(cells),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(analyze(args.matrix, args.annotation), indent=2, allow_nan=False))
