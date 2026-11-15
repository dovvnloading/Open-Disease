"""Decompose a frozen expression score with a separate, batch-adjusted cycle reference."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from bulk_gse153739 import load_mumbai
from cross_cohort_discovery import align, assignment_design
from genome_rome import load_rome
from genome_signature_transfer import weights_and_candidates
from melbourne_signature_phase import MATRIX, META_SHA, clinical, clinical_digest, read_profiles
from rome_scan_design import scan_design
from swiss_epithelial_mapping import HGNC_FILE
from swiss_genome_signature import load_glandular_profiles


def design(stages, batches, reference):
    phase = [k for k in sorted(set(stages)) if k != reference]
    columns = [batches == b for b in sorted(set(batches))] + [stages == k for k in phase]
    return np.column_stack(columns).astype(float), phase


def cycle_components(x, stages, batches, weights):
    """Estimate six batch-adjusted stage contrasts and project the original weights."""
    model, levels = design(stages, batches, 1)
    coef, _, rank, _ = np.linalg.lstsq(model, x, rcond=None)
    if rank != model.shape[1] or len(levels) != 6:
        raise ValueError("Stage and batch are aliased or stages missing")
    contrasts = coef[-len(levels):]
    _, singular, right = np.linalg.svd(contrasts, full_matrices=False)
    if np.linalg.matrix_rank(contrasts) != 6:
        raise ValueError("Stage reference has fewer than six independent contrasts")
    projected = right.T @ (right @ weights)
    residual = weights - projected
    norm = np.abs(weights).sum()
    components = np.column_stack([weights, projected, residual]) / norm
    # Change baseline to stage 7 and use a QR projection to test parameterization invariance.
    alternate, _ = design(stages, batches, 7)
    other = np.linalg.lstsq(alternate, x, rcond=None)[0][-len(levels):]
    q = np.linalg.qr(other.T, mode="reduced")[0]
    if not np.allclose(q @ (q.T @ weights), projected, rtol=0, atol=1e-9):
        raise ValueError("Projection depends on stage parameterization")
    if np.max(np.abs(contrasts @ components[:, 2])) > 1e-9:
        raise ValueError("Residual retains an in-sample modeled stage contrast")
    mid = contrasts[levels.index(6)] - contrasts[levels.index(3)]
    return components, mid, singular


def summarize(contrast, component):
    full, cycle, rest = map(float, contrast @ component)
    if not np.isclose(full, cycle + rest, rtol=0, atol=1e-10):
        raise ValueError("Original score failed the reconstruction check")
    return dict(full=full, cycle_aligned=cycle, residual=rest)


def analyze(data, rds, workbook, rscript):
    hgnc = data / "gse266265" / HGNC_FILE
    records = clinical(data / "gse234354" / "meta")
    if clinical_digest(records) != META_SHA:
        raise ValueError("Melbourne sample-stage mapping changed")
    mel_genes, mel_ranks, stages, batches, _ = read_profiles(data / "gse234354" / MATRIX, hgnc, records)
    rg, rx, rc, _ = load_rome(data / "gse25628")
    scan_cases, dates, _ = scan_design(data / "gse25628")
    if not np.array_equal(rc, scan_cases):
        raise ValueError("Inconsistent Rome cases")
    keep = np.isin(dates, np.intersect1d(dates[rc], dates[~rc]))
    mg, mx, mc, _ = load_mumbai(data / "gse153739", hgnc)
    swiss = load_glandular_profiles(rds, workbook, hgnc, rscript)
    cohorts = {"Rome": (rg, rx[keep], rc[keep], dates[keep]),
               "Mumbai": (mg, mx, mc, np.zeros(len(mc), dtype=int)), "Swiss_glandular": swiss}
    reference, weights, _ = weights_and_candidates()
    common = np.array(sorted(set(reference).intersection(mel_genes, *(r[0] for r in cohorts.values()))))
    w = weights[np.searchsorted(reference, common)]
    matrix = align(mel_genes, mel_ranks, common)
    component, mid, singular = cycle_components(matrix, stages, batches, w)
    output = {}
    for name, (gene, x, case, block) in cohorts.items():
        x = align(gene, x, common)
        obs, _, _ = assignment_design(case, block)
        gc = obs @ x
        if not np.allclose(gc @ component, obs @ (x @ component), rtol=0, atol=1e-10):
            raise ValueError("Gene-first and donor-first component contrasts disagree")
        output[name] = summarize(gc, component)
    return dict(common_genes=len(common), original_genes=len(reference),
                coverage_absolute_weight=float(abs(w).sum()/abs(weights).sum()),
                stage_contrast_singular_values=singular.tolist(),
                cycle_fraction_squared_weight_norm=float(np.sum(component[:,1]**2)/np.sum(component[:,0]**2)),
                Melbourne_stage6_minus_stage3_full_model=summarize(mid, component), contrasts=output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("rds", type=Path)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    print(json.dumps(analyze(args.data, args.rds, args.workbook, args.rscript), indent=2, allow_nan=False))
