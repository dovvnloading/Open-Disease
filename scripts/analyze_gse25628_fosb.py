"""Reproduce the prespecified GSE25628 comparison from two pinned local files.

Print one aggregate result. No network, participant identifiers or array values
are emitted. See docs/research-notes/GSE25628_FOSB_PLAN_2026-11-13.md.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import re

import numpy as np
from scipy.stats import ttest_ind

from open_disease.geo_microarray import sha256_file

PROBE = "202768_at"
MATRIX_HASH = "3d77aae9769afb1f4bb1d03e2bc73fbab4ccc52ac5cd0f0945e4c0224fcca089"
ANNOTATION_HASH = "9f5a482f9f3d478b47fc80124600319e569f1e9ab447f3f256d61b1fbacce4fa"
GROUP_CODES = {
    ("patient", "Normal (eutopic)"): "case",
    ("patient", "Pathological (ectopic)"): "lesion",
    ("healthy donor", "Normal (control)"): "control",
}
COUNTS = {"case": 8, "control": 6, "lesion": 8}


def checked_groups(header):
    if header["!Series_geo_accession"] != [["GSE25628"]] or header["!Series_platform_id"] != [["GPL571"]]:
        raise ValueError("Unexpected GEO series or platform")
    ids = header["!Sample_geo_accession"]
    if len(ids) != 1 or len(ids[0]) != 22 or len(set(ids[0])) != 22 or not all(
        re.fullmatch(r"GSM[1-9][0-9]*", obs) for obs in ids[0]
    ):
        raise ValueError("Expected exactly 22 distinct array accessions")
    if header["!Sample_platform_id"] != [["GPL571"] * 22]:
        raise ValueError("Sample platform mismatch")
    processing = header["!Sample_data_processing"]
    if len(processing) != 1 or len(processing[0]) != 22 or not all(
        "RMA + quantile normalization" in entry for entry in processing[0]
    ):
        raise ValueError("Unexpected processing")
    rows = header["!Sample_characteristics_ch1"]
    if not rows or any(len(row) != 22 for row in rows):
        raise ValueError("Misaligned source characteristics")
    groups = []
    for column in zip(*rows):
        fields = {}
        for entry in column:
            key, sep, value = entry.partition(": ")
            if not sep or key in fields:
                raise ValueError("Malformed or repeated source characteristic")
            fields[key] = value
        group = GROUP_CODES.get((fields.get("donor"), fields.get("disease state")))
        if fields.get("tissue") != "endometrium" or group is None:
            raise ValueError("Unrecognized tissue or clinical group")
        groups.append(group)
    if Counter(groups) != COUNTS:
        raise ValueError("Unexpected case, control or lesion count")
    return ids[0], groups


def check_annotation(source):
    with gzip.open(source, "rt", newline="") as stream:
        for line in stream:
            if line.strip() == "!platform_table_begin":
                break
        rows = csv.DictReader(stream, delimiter="\t")
        target = [(r["ID"], r["Gene symbol"], r["Gene ID"]) for r in rows
                  if r["ID"] == PROBE or r["Gene symbol"] == "FOSB"]
    if target != [(PROBE, "FOSB", "2354")]:
        raise ValueError("FOSB probe mapping differs from the frozen annotation")


def read_target(source):
    with gzip.open(source, "rt", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = defaultdict(list)
        for row in reader:
            if row == ["!series_matrix_table_begin"]:
                break
            if row and row[0].startswith(("!Series_", "!Sample_")):
                header[row[0]].append(row[1:])
        ids, groups = checked_groups(header)
        if next(reader, None) != ["ID_REF", *ids]:
            raise ValueError("Expression columns disagree with the audited sample order")
        probes = set()
        observed = None
        for row in reader:
            if row == ["!series_matrix_table_end"]:
                break
            if len(row) != 23 or not row[0] or row[0] in probes:
                raise ValueError("Invalid feature row or duplicate probe")
            probes.add(row[0])
            if row[0] == PROBE:
                observed = {"case": [], "control": []}
                for group, value in zip(groups, row[1:]):
                    if group != "lesion":
                        try:
                            number = float(value)
                        except ValueError:
                            raise ValueError("Missing or invalid eligible FOSB intensity") from None
                        if not np.isfinite(number):
                            raise ValueError("Nonfinite eligible FOSB intensity")
                        observed[group].append(number)
        else:
            raise ValueError("Missing end of matrix")
        if any(any(row) for row in reader):
            raise ValueError("Unexpected trailing matrix content")
    if observed is None:
        raise ValueError("Frozen probe is absent")
    return observed, len(probes)


def summarize(case, control):
    case, control = np.asarray(case, dtype=float), np.asarray(control, dtype=float)
    difference = float(case.mean() - control.mean())
    test = ttest_ind(case, control, equal_var=False)
    low, high = map(float, test.confidence_interval())
    if not np.isfinite([difference, low, high, test.pvalue, test.df]).all():
        raise ValueError("The prespecified Welch contrast cannot be estimated")
    leave_one_out = np.concatenate([
        (case.sum() - case) / (len(case) - 1) - control.mean(),
        case.mean() - (control.sum() - control) / (len(control) - 1),
    ])
    decision = "inconclusive"
    if low > 0:
        decision = "supports_expected_direction"
    elif high < 0:
        decision = "supports_opposite_direction"
    return {
        "difference_log2": difference, "ci95_log2": [low, high],
        "ratio_geometric_means": float(2 ** difference),
        "welch_df": float(test.df), "two_sided_p": float(test.pvalue),
        "leave_one_out_range_log2": [float(leave_one_out.min()), float(leave_one_out.max())],
        "decision": decision,
    }


def analyze(matrix, annotation):
    if sha256_file(matrix) != MATRIX_HASH or sha256_file(annotation) != ANNOTATION_HASH:
        raise ValueError("Input SHA-256 differs from the audited snapshot")
    check_annotation(annotation)
    observed, feature_count = read_target(matrix)
    return {
        "series": "GSE25628", "platform": "GPL571", "probe": PROBE,
        "counts": COUNTS, "features_in_matrix": feature_count, "features_tested": 1,
        **summarize(observed["case"], observed["control"]),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(analyze(args.matrix, args.annotation), indent=2, allow_nan=False))
