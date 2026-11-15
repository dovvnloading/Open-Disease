"""Outcome-informed Finnish tests of the Melbourne reference patterns, with a phase control."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from openpyxl import load_workbook

from analyze_melbourne_genes import SOURCE_HASH
from cross_cohort_discovery import align, assignment_design, holm
from genome_endometdb import load_finland
from genome_rome import load_rome
from melbourne_pattern_transfer import (
    FINNISH_TABLE_HASH, SHEETS, exact_correlations, signature, standardize, symbol_bridge,
)
from open_disease.geo_microarray import sha256_file


def within_strata(scores, factors):
    """Standardize features by their exact-null SD with the project's stratified design."""
    centered = np.empty_like(scores)
    variance = np.zeros(scores.shape[1])
    for ids, factor in factors:
        values = scores[ids] - scores[ids[0]]
        centered[ids] = values - values.mean(axis=0)
        variance += values.var(axis=0, ddof=1) * factor
    if np.any(variance <= 0):
        raise ValueError("Nonvariable reference gene in a stratified comparison")
    return centered / np.sqrt(variance)


def evaluate(scores, labels, strata, ref):
    observed, weights, factors = assignment_design(labels, strata)
    z = within_strata(scores, factors)
    statistic = observed @ z
    r = float(np.corrcoef(statistic, ref)[0, 1])
    dist = exact_correlations(z, ref, weights)
    extreme = int((np.abs(dist) >= abs(r) - 1e-12).sum())
    if not np.any(np.isclose(dist, r, rtol=0, atol=1e-12)):
        raise ValueError("Observed correlation omitted from null distribution")
    return dict(cases_or_secretory=int(labels.sum()), controls_or_proliferative=int((~labels).sum()),
                strata=len(factors), genes=len(ref), assignments=len(weights), r=r, extreme=extreme,
                p_two_sided=extreme / len(weights), null_central_95=np.quantile(dist, [.025,.975]).tolist()), statistic


def include_mixed(labels, strata):
    keep = np.zeros(len(labels), dtype=bool)
    for st in np.unique(strata):
        group = strata == st
        if np.any(labels[group]) and np.any(~labels[group]):
            keep |= group
    return keep


def partial_r(x, y, axis):
    intercept_and_axis = np.column_stack([np.ones(len(axis)), axis])
    matrix = np.column_stack([x, y])
    residuals = matrix - intercept_and_axis @ np.linalg.lstsq(intercept_and_axis, matrix, rcond=None)[0]
    return float(np.corrcoef(residuals.T)[0, 1])


def run(rome_dir, workbook, finnish, native):
    if sha256_file(workbook) != SOURCE_HASH or sha256_file(finnish) != FINNISH_TABLE_HASH:
        raise ValueError("Unexpected input checksum")
    with finnish.open(newline="") as f:
        fin_published = {int(row["id"]): float(row["statistic"]) for row in csv.DictReader(f, delimiter="\t")}
    fg, fs, fc, strata = native
    if len(fc) != 25 or fc.sum() != 12 or len(fin_published) != 12044:
        raise ValueError("Unexpected Finnish design or shared universe")
    runs, phases = np.array([str(x).split("-") for x in strata]).T
    if set(phases) != {"secretory", "proliferative"}:
        raise ValueError("Unknown Finnish histology")
    rg, rs, rc, _ = load_rome(rome_dir)
    bridge = symbol_bridge(rome_dir)
    book = load_workbook(workbook, read_only=True, data_only=True)
    try:
        refs = {s: signature(book[s], bridge, set(fin_published))[:2] for s in SHEETS}
    finally:
        book.close()
    results = {}
    full = {}
    for name, (ids, t) in refs.items():
        score = align(fg, fs, ids)
        full[name], stat = evaluate(score, fc, strata, t)
        stored = np.array([fin_published[i] for i in ids])
        if not np.allclose(stat, stored, rtol=0, atol=1e-11):
            raise ValueError("Reconstructed Finnish statistics disagree with the published screen")
        full[name]["max_published_abs_diff"] = float(np.max(np.abs(stat - stored)))
    for row, p in zip(full.values(), holm([row["p_two_sided"] for row in full.values()]), strict=True):
        row["p_holm_two_references"] = float(p)
    results["disease_contrast"] = full
    phase_ids, phase_ref = refs["TableS3"]
    panel = align(fg, fs, phase_ids)
    splits = {}
    for label in ("proliferative","secretory"):
        keep = phases == label
        splits[label], _ = evaluate(panel[keep], fc[keep], strata[keep], phase_ref)
    for row, p in zip(splits.values(), holm([row["p_two_sided"] for row in splits.values()]), strict=True):
        row["p_holm_two_phase_splits"] = float(p)
    results["disease_contrast_by_phase"] = splits

    secretory = phases == "secretory"
    control_strata = np.array([f"{run}-{'case' if affected else 'control'}"
                               for run, affected in zip(runs, fc, strict=True)])
    keep = include_mixed(secretory, control_strata)
    results["positive_control_secretory_vs_proliferative"], _ = evaluate(
        panel[keep], secretory[keep], control_strata[keep], phase_ref)
    results["positive_control_secretory_vs_proliferative"]["from_donors"] = len(fc)
    results["positive_control_secretory_vs_proliferative"]["included_cases"] = int(fc[keep].sum())
    results["positive_control_secretory_vs_proliferative"]["included_controls"] = int((~fc[keep]).sum())
    disease_ids, disease_ref = refs["TableS20"]
    lookup = dict(zip(disease_ids, disease_ref, strict=True))
    ds = np.array([lookup[i] for i in phase_ids])
    fs3 = np.array([fin_published[i] for i in phase_ids])
    r3 = standardize(align(rg, rs, phase_ids))
    rome_contrast = r3[rc].mean(axis=0) - r3[~rc].mean(axis=0)
    results["same_gene_descriptive_sensitivity"] = dict(
        genes=len(phase_ids),
        fin_vs_melbourne_disease_r=float(np.corrcoef(fs3, ds)[0, 1]),
        rome_vs_melbourne_disease_r=float(np.corrcoef(rome_contrast, ds)[0, 1]),
        fin_vs_rome_disease_r=float(np.corrcoef(fs3, rome_contrast)[0, 1]),
        fin_vs_rome_partial_r_after_linear_phase_axis=partial_r(fs3, rome_contrast, phase_ref),
        melbourne_disease_vs_phase_r=float(np.corrcoef(ds, phase_ref)[0, 1]))
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rome", type=Path, required=True)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--finnish", type=Path, required=True)
    parser.add_argument("--finland", type=Path, required=True)
    args = parser.parse_args()
    out = run(args.rome, args.workbook, args.finnish, load_finland(args.finland))
    print(json.dumps(out, indent=2))
