"""Describe the locked ten candidates by Finnish cycle phase, retaining the completed Rome contrast."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from cross_cohort_discovery import align, write_tsv
from genome_endometdb import load_finland
from open_disease.geo_microarray import sha256_file

LOCK_HASH = "ba140f2bb026b0416e1c705ab28f616f430cc4295f94d174feaf5d15947c1e31"
VALIDATION_HASH = "fb79bdef1268bbb50d610574bef0d0ecbe182e61f7562e806980dd9061e130e5"


def phase_effects(values, cases, strata):
    """Equal-weight separately within phases, retaining both groups inside every run."""
    effects = {"proliferative": [], "secretory": []}
    for cell in np.unique(strata):
        _, phase = cell.split("-")
        selected = strata == cell
        if phase not in effects or not cases[selected].any() or cases[selected].all():
            raise ValueError("Unknown phase or stratum without case-control overlap")
        delta = values[selected & cases].mean(axis=0) - values[selected & ~cases].mean(axis=0)
        effects[phase].append(delta)
    return {phase: (len(x), np.mean(x, axis=0)) for phase, x in effects.items() if x}


def analyze(finland, lock_path, validation_path):
    if sha256_file(lock_path) != LOCK_HASH or sha256_file(validation_path) != VALIDATION_HASH:
        raise ValueError("The published nomination or Rome validation table changed")
    candidates = json.loads(lock_path.read_text())["candidates"]
    with validation_path.open(newline="") as fh:
        validation = list(csv.DictReader(fh, delimiter="\t"))
    if [r["id"] for r in validation] != [r["id"] for r in candidates]:
        raise ValueError("Validation rows do not match the locked candidates")
    genes, scores, cases, strata = load_finland(finland)
    columns = []
    for row in candidates:
        ids = [int(row["id"])] if row["family"] == "gene" else row["entrez_ids"]
        columns.append(align(genes, scores, ids).mean(axis=1))
    effect = phase_effects(np.column_stack(columns), cases, strata)
    n_prolif, prolif = effect["proliferative"]
    n_sec, secretory = effect["secretory"]
    if (n_prolif, n_sec) != (4, 3):
        raise ValueError("The frozen Finnish cohort's run-by-phase structure changed")
    combined = (n_prolif * prolif + n_sec * secretory) / (n_prolif + n_sec)
    rows = []
    for i, (candidate, val) in enumerate(zip(candidates, validation, strict=True)):
        direction = "higher" if combined[i] > 0 else "lower"
        if (abs(combined[i] - candidate["effect_pp"]) > 1e-9 or direction != candidate["direction"]
                or direction != val["direction"] or abs(float(val["finnish_effect_pp"]) - combined[i]) > 1e-9):
            raise ValueError("Phase results do not reconstruct the original Finnish effect")
        rome = float(val["rome_effect_pp"])
        rows.append({
            "family": candidate["family"], "id": candidate["id"],
            "proliferative_strata": n_prolif, "finnish_proliferative_pp": float(prolif[i]),
            "secretory_strata": n_sec, "finnish_secretory_pp": float(secretory[i]),
            "finnish_pooled_pp": float(combined[i]), "rome_proliferative_pp": rome,
            "proliferative_concordant_with_rome": bool(np.sign(prolif[i]) == np.sign(rome)),
            "secretory_concordant_with_rome": bool(np.sign(secretory[i]) == np.sign(rome)),
        })
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--finland", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--validation", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.finland, args.lock, args.validation)
    write_tsv(args.output, result[0].keys(), result)
    pathways = [r for r in result if r["family"] == "pathway"]
    reversed_pathways = [r for r in pathways if r["finnish_pooled_pp"] * r["rome_proliferative_pp"] < 0]
    print(json.dumps({
        "candidates": len(result), "reversed_pathways": len(reversed_pathways),
        "positive_finnish_proliferative_pathways": sum(r["finnish_proliferative_pp"] > 0 for r in reversed_pathways),
    }))
