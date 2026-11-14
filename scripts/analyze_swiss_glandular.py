"""Run the prospectively defined GSE266265 glandular pathway analysis.

Takes three cached public sources. All intermediate feature, sample and pathway
tables are ephemeral. The durable JSON contains only aggregate estimates.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np

from audit_gse266265_epithelial import verify, write_sample_map
from cross_cohort_discovery import assignment_design, holm
from open_disease.geo_microarray import sha256_file
from swiss_epithelial_mapping import load_projection, write_projection

PLAN = "7f380a743491e98036b200bd6c2d097d7bd489cd"
PROJECTION_SHA256 = "6609ef4bb35608c7f791391b54344e0b0e9f77a44ed3004d95925aed2779a0c2"
ENDPOINTS = ("GLOBAL", "R-HSA-9008059", "R-HSA-5654741", "R-HSA-163560",
             "R-HSA-5654710", "R-HSA-432722")
RSCRIPTS = Path(__file__).resolve().parent / "r"


def read_scores(path):
    """Read only a complete, internally consistent donor-by-endpoint design."""
    records, design = {}, {}
    with path.open(encoding="utf-8", newline="") as source:
        table = csv.DictReader(source, delimiter="\t")
        if table.fieldnames != ["donor_index", "capture", "case", "endpoint", "primary", "cell_first"]:
            raise ValueError("Unexpected score fields")
        for row in table:
            donor, endpoint = row["donor_index"], row["endpoint"]
            case = {"TRUE": True, "FALSE": False, "1": True, "0": False}.get(row["case"])
            batch = row["capture"]
            if (not donor or not batch or case is None or endpoint not in ENDPOINTS
                    or (donor, endpoint) in records):
                raise ValueError("Invalid donor, endpoint or design label")
            if design.setdefault(donor, (batch, case)) != (batch, case):
                raise ValueError("Inconsistent donor design")
            records[donor, endpoint] = float(row["primary"]), float(row["cell_first"])
    donors = list(design)
    if not donors or len(records) != len(donors) * len(ENDPOINTS):
        raise ValueError("Missing donor or endpoint")
    values = np.array([[records[d, e] for e in ENDPOINTS] for d in donors])
    if not np.isfinite(values).all() or np.any((values < 0) | (values > 100)):
        raise ValueError("Scores outside finite percentile range")
    if not np.allclose(values[:, 0], values[:, 1:].mean(axis=1), rtol=0, atol=1e-10):
        raise ValueError("Global score differs from its five components")
    return (np.array([design[d][1] for d in donors]), np.array([design[d][0] for d in donors]),
            values[:, :, 0], values[:, :, 1])


def summarize(case, capture, primary, cell_first):
    observed, assignments, groups = assignment_design(case, capture)
    if len(groups) < 2:
        raise ValueError("At least two mixed batches required")
    differences = observed @ primary
    tail = (abs(assignments @ primary) >= abs(differences) - 1e-10).sum(axis=0)
    n = len(assignments)
    if not np.all((tail >= 1) & (tail <= n)):
        raise ValueError("Invalid exact tail count")
    p_values = tail / n
    adjusted = holm(p_values[1:])
    cell_differences = observed @ cell_first
    batch_differences = np.array([primary[ids][case[ids]].mean(axis=0)
                                  - primary[ids][~case[ids]].mean(axis=0) for ids, _ in groups])
    omitted = (batch_differences.sum(axis=0) - batch_differences) / (len(groups) - 1)
    endpoints = {}
    for i, key in enumerate(ENDPOINTS):
        endpoints[key] = {
            "difference": float(differences[i]), "extreme": int(tail[i]), "p": float(p_values[i]),
            "holm": None if i == 0 else float(adjusted[i - 1]),
            "supported": bool(p_values[0] <= 0.05 and (i == 0 or adjusted[i - 1] <= 0.05)),
            "batch_omit_low": float(omitted[:, i].min()), "batch_omit_high": float(omitted[:, i].max()),
            "cell_first": float(cell_differences[i]),
        }
    return {"protocol": PLAN, "donors": len(case), "cases": int(case.sum()),
            "controls": int((~case).sum()), "batches": len(groups), "assignments": n,
            "endpoints": endpoints}


def analyze(rds, supplement, annotation_directory, rscript):
    verify(rds)
    env = dict(os.environ, R_DEFAULT_PACKAGES="NULL", OPENBLAS_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(dir=rds.parent) as scratch:
        private = Path(scratch)
        labels, projection = private / "features.txt", private / "projection.tsv"
        sample_map, scores = private / "design.tsv", private / "scores.tsv"
        write_sample_map(supplement, sample_map)
        subprocess.run([rscript, RSCRIPTS / "extract_swiss_feature_labels.R", rds, labels],
                       env=env, check=True)
        rows, _ = load_projection(annotation_directory, labels_path=labels)
        write_projection(projection, rows)
        if sha256_file(projection) != PROJECTION_SHA256:
            raise ValueError("Swiss gene projection differs from the pre-outcome projection")
        measured = subprocess.run([rscript, RSCRIPTS / "score_swiss_glandular.R", rds,
                                   sample_map, projection, scores], env=env, check=True,
                                  stdout=subprocess.PIPE, text=True)
        cell_counts = re.findall(r"^glandular_cells=(\d+)$", measured.stdout, re.MULTILINE)
        if len(cell_counts) != 1 or int(cell_counts[0]) < 13 * 20:
            raise ValueError("Expected one aggregate glandular-cell count")
        result = summarize(*read_scores(scores))
        if (result["donors"], result["cases"], result["controls"], result["batches"],
                result["assignments"]) != (13, 6, 7, 4, 144):
            raise ValueError("Donor design differs from the frozen study")
        result["cells"] = int(cell_counts[0])
    return result


def compare_independent(result, other):
    """Require agreement in exact counts and all prespecified continuous estimates."""
    for key in ("donors", "cases", "controls", "batches", "assignments", "cells"):
        if result[key] != other[key]:
            raise ValueError(f"Independent design mismatch for {key}")
    if set(other["endpoints"]) != set(ENDPOINTS):
        raise ValueError("Independent endpoint family mismatch")
    for endpoint, row in result["endpoints"].items():
        check = other["endpoints"][endpoint]
        for key in ("extreme", "holm", "p", "difference", "batch_omit_low", "batch_omit_high", "cell_first"):
            if row[key] is None or check[key] is None or key == "extreme":
                matches = row[key] == check[key]
            else:
                matches = np.isclose(row[key], check[key], atol=1e-10, rtol=0)
            if not matches:
                raise ValueError(f"Independent results differ for {endpoint}, {key}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rds", type=Path)
    parser.add_argument("supplement", type=Path)
    parser.add_argument("annotation_directory", type=Path)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--compare", type=Path, help="Optional JSON from the independent source verifier")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.rds, args.supplement, args.annotation_directory, args.rscript)
    if args.compare:
        compare_independent(result, json.loads(args.compare.read_text(encoding="utf-8")))
    content = json.dumps(result, indent=2, allow_nan=False) + "\n"
    with tempfile.TemporaryDirectory(dir=args.output.parent) as scratch:
        temporary = Path(scratch) / "aggregate.json"
        temporary.write_text(content, encoding="utf-8")
        temporary.replace(args.output)


if __name__ == "__main__":
    main()
