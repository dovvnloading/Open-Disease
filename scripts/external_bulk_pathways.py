"""Test the pre-outcome five-pathway endpoint in two additional biopsy series."""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from bulk_gse153739 import load_mumbai
from bulk_gse6364 import load_us_both
from cross_cohort_discovery import align, assignment_design, holm
from open_disease.geo_microarray import sha256_file
from swiss_epithelial_mapping import HGNC_FILE, LOCK_FILE, LOCK_SHA256
from swiss_global_interval import allocation_contrasts, direct_tail, shift_limits

PROTOCOL = "fe541f1e113de0b41cfc283b9e16ae2456a28f5c"
SWISS_MAPPING = "6609ef4bb35608c7f791391b54344e0b0e9f77a44ed3004d95925aed2779a0c2"


def locked_members(swiss_projection):
    if sha256_file(LOCK_FILE) != LOCK_SHA256 or sha256_file(swiss_projection) != SWISS_MAPPING:
        raise ValueError("Unexpected original nominations or Swiss projection")
    rows = json.loads(LOCK_FILE.read_text())["candidates"]
    locked = {r["id"]: r["entrez_ids"] for r in rows if r["family"] == "pathway"}
    swiss = {k: set() for k in locked}
    with swiss_projection.open(encoding="utf-8") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["family"] == "pathway" and row["row_index"]:
                swiss[row["id"]].add(int(row["entrez_id"]))
    return locked, swiss


def project(locked, swiss, measured):
    common, restricted = {}, {}
    for k, members in locked.items():
        common[k] = sorted(measured.intersection(members))
        if len(common[k]) < 15 or len(common[k]) * 5 < 4 * len(members):
            raise ValueError("Insufficient common measured membership for " + k)
        restricted[k] = sorted(swiss[k].intersection(common[k]))
        if not restricted[k]:
            raise ValueError("No shared Swiss genes in " + k)
    return common, restricted


def endpoints(genes, scores, membership):
    paths = np.column_stack([align(genes, scores, members).mean(axis=1) for members in membership.values()])
    return np.column_stack([paths, paths.mean(axis=1)])


def difference(scores, cases):
    return scores[cases].mean(axis=0) - scores[~cases].mean(axis=0)


def evaluate(dataset, common, restricted, *, inference):
    genes, primary, cases, alternative = dataset
    if (cases.dtype != np.bool_ or min(sum(cases), sum(~cases)) < 2
            or primary.shape != alternative.shape or primary.shape != (len(cases), len(genes))
            or not np.isfinite(primary).all() or not np.isfinite(alternative).all()):
        raise ValueError("Invalid aligned study design or scores")
    labels = [*common, "GLOBAL"]
    scores = endpoints(genes, primary, common)
    effect = difference(scores, cases)
    alternate = difference(endpoints(genes, alternative, common), cases)
    swiss = difference(endpoints(genes, primary, restricted), cases)
    omit = []
    for donor in range(len(cases)):
        keep = np.arange(len(cases)) != donor
        omit.append(difference(scores[keep], cases[keep]))
    low, high = np.min(omit, axis=0), np.max(omit, axis=0)
    rows = {k: {"difference": float(effect[i]), "loo_low": float(low[i]), "loo_high": float(high[i]),
                "alternate": float(alternate[i]), "swiss_members": float(swiss[i])}
            for i, k in enumerate(labels)}
    interval = None
    if inference:
        strata = np.zeros(len(cases), dtype=int)
        _, weights, _ = assignment_design(cases, strata)
        extreme = np.sum(np.abs(weights @ scores) >= np.abs(effect) - 1e-10, axis=0)
        if (extreme < 1).any():
            raise ValueError("Observed assignment absent from permutation tail")
        p = extreme / len(weights)
        adjusted = holm(p[:-1])
        for i, k in enumerate(labels):
            rows[k].update(extreme=int(extreme[i]), p=float(p[i]))
            if k != "GLOBAL":
                rows[k]["holm"] = float(adjusted[i])
        # The interval inverts the very same rank-score reference. Exact rational
        # arithmetic here is for the donor scores computed by this implementation.
        outcomes = [Fraction(float(x)) for x in scores[:, -1]]
        disease = cases.astype(int).tolist()
        blocks = strata.tolist()
        observed, allocations = allocation_contrasts(outcomes, disease, blocks)
        left, right, required = shift_limits(observed, allocations)
        checks = {k: direct_tail(outcomes, disease, blocks, x) for k, x in
                  (("zero", Fraction(0)), ("lower", left), ("upper", right)) if x is not None}
        if checks["zero"] != rows["GLOBAL"]["extreme"] or any(
                n < required for k, n in checks.items() if k != "zero"):
            raise ValueError("Interval failed its exact boundary or null check")
        interval = {"bounds_pp": [None if x is None else float(x) for x in (left, right)],
                    "minimum_extreme_for_inclusion": required, "tail_checks": checks,
                    "assignments": len(weights)}
    return rows, interval


def analyze(data, annotation):
    locked, swiss = locked_members(LOCK_FILE.parents[1] / "external-bulk-2026-11-14/locked_hgnc_projection.tsv")
    mumbai = load_mumbai(data / "gse153739", data / "gse266265" / HGNC_FILE)
    original_genes = sorted(set().union(*(set(v) for v in locked.values())))
    us, scan = load_us_both(data / "gse6364", annotation, original_genes)
    if (len(us[2]), int(us[2].sum()), len(mumbai[2]), int(mumbai[2].sum())) != (11, 6, 7, 4):
        raise ValueError("Different number of eligible case and control specimens")
    if scan["conditional_assignments"] != 1 or scan["mixed_scan_days"] != 0:
        raise ValueError("Source scan-date design differs from the frozen U.S. assessment")
    common, restricted = project(locked, swiss, set(us[0]) & set(mumbai[0]))
    us_rows, _ = evaluate(us, common, restricted, inference=False)
    mumbai_rows, interval = evaluate(mumbai, common, restricted, inference=True)
    return {
        "protocol": PROTOCOL,
        "membership": {k: {"primary": len(common[k]), "swiss": len(restricted[k]),
                           "locked": len(locked[k]), "entrez_ids": common[k],
                           "swiss_entrez_ids": restricted[k]} for k in common},
        "GSE6364": us_rows, "GSE153739": mumbai_rows, "mumbai_global_interval": interval,
        "us_scan_date_assignments": scan["conditional_assignments"],
        "mumbai_global_gate_passed": bool(mumbai_rows["GLOBAL"]["p"] <= 0.05),
        "mumbai_supported_components": [k for k in common if mumbai_rows["GLOBAL"]["p"] <= 0.05
                                       and mumbai_rows[k]["holm"] <= 0.05],
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=Path("data/private"))
    ap.add_argument("--gpl570", type=Path, default=Path("data/raw/GPL570.annot.gz"))
    args = ap.parse_args()
    print(json.dumps(analyze(args.data, args.gpl570), indent=2, allow_nan=False))
