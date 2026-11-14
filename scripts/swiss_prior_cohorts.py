"""Compare full vs Swiss-measured pathway membership on the original two arrays.

Reuses both pinned cohort loaders, their platform-specific donor rank universes
and their original donor contrast. Writes only the five aggregate comparisons.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from cross_cohort_discovery import align
from genome_endometdb import load_finland
from genome_rome import load_rome
from open_disease.geo_microarray import sha256_file
from swiss_epithelial_mapping import LOCK_FILE, load_projection

VALIDATION = LOCK_FILE.with_name("validation.tsv")
VALIDATION_SHA256 = "fb79bdef1268bbb50d610574bef0d0ecbe182e61f7562e806980dd9061e130e5"


def donor_contrast(scores, cases, strata):
    """Average the case-minus-control donor mean differences equally across strata."""
    contrasts = []
    for level in np.unique(strata):
        block = strata == level
        if not (np.any(block & cases) and np.any(block & ~cases)):
            raise ValueError("Prior cohort stratum lacks group overlap")
        contrasts.append(scores[block & cases].mean(axis=0) - scores[block & ~cases].mean(axis=0))
    return np.mean(contrasts, axis=0)


def pathway_effects(dataset, memberships):
    genes, scores, cases, strata = dataset
    features = np.column_stack([
        align(genes, scores, np.asarray(ids, dtype=int)).mean(axis=1) for ids in memberships
    ])
    return donor_contrast(features, cases, strata)


def run(data):
    projected, _ = load_projection(data / "gse266265")
    if sha256_file(VALIDATION) != VALIDATION_SHA256:
        raise ValueError("Published earlier full-membership estimates differ")
    with VALIDATION.open(encoding="utf-8", newline="") as source:
        baseline = {r["id"]: r for r in csv.DictReader(source, delimiter="\t") if r["family"] == "pathway"}
    full, measured = {}, {}
    for family, key, gene, _, index in projected:
        if family != "pathway":
            continue
        full.setdefault(key, []).append(int(gene))
        measured.setdefault(key, [])
        if index:
            measured[key].append(int(gene))
    if len(full) != 5 or full.keys() != baseline.keys():
        raise ValueError("Projection differs from prior five pathways")
    keys = list(full)
    report = {
        key: {"full_members": len(full[key]), "matched_members": len(measured[key])}
        for key in keys
    }
    # Compute full and matched estimates in the same run, to verify the cached
    # arrays, phenotype contrasts and feature ranks against the published results.
    memberships = [full[key] for key in keys] + [measured[key] for key in keys]
    for cohort, loader, path, expected_ncase, expected_ncontrol, field in (
        ("finland", load_finland, "gse141549", 12, 13, "finnish_effect_pp"),
        ("rome", load_rome, "gse25628", 8, 6, "rome_effect_pp"),
    ):
        dataset = loader(data / path)
        cases = dataset[2]
        if cases.sum() != expected_ncase or (~cases).sum() != expected_ncontrol:
            raise ValueError(f"{cohort} donor counts differ from the original cohorts")
        effects = pathway_effects(dataset, memberships)
        for i, key in enumerate(keys):
            full_effect, matched_effect = float(effects[i]), float(effects[len(keys) + i])
            if not np.isclose(full_effect, float(baseline[key][field]), atol=1e-9, rtol=0):
                raise ValueError(f"{cohort} full effect differs from the published result")
            report[key][cohort + "_full_pp"] = full_effect
            report[key][cohort + "_matched_pp"] = matched_effect
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path, help="Parent directory of all three cohorts")
    parser.add_argument("--out", type=Path, required=True, help="Output JSON of aggregate estimates")
    args = parser.parse_args()
    report = run(args.data)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
