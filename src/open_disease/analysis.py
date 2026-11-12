"""Exploratory donor-level baseline: stratified Welch tests of log2(1+CPM).

Not a replacement for carefully modeled edgeR/DESeq2 analyses, confounder
adjustment, or independent clinical research.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from .schema import DONOR, validate_pseudobulk

STRATA = ["study_id", "phase", "cell_type"]


def bh_adjust(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    if p.ndim != 1 or not np.isfinite(p).all() or ((p < 0) | (p > 1)).any():
        raise ValueError("P values must be finite and between zero and one")
    if not len(p):
        return p.copy()
    order = np.argsort(p, kind="stable")
    ranks = np.arange(1, len(p) + 1)
    corrected = np.minimum.accumulate((p[order] * len(p) / ranks)[::-1])[::-1]
    q = np.empty(len(p), dtype=float)
    q[order] = np.clip(corrected, 0, 1)
    return q


def analyze_pseudobulk(data: pd.DataFrame, *, min_donors_per_group: int = 3) -> pd.DataFrame:
    if min_donors_per_group < 3:
        raise ValueError("Baseline requires at least three donors per diagnosis group")
    df = validate_pseudobulk(data)
    results = []
    for (study, phase, cell_type), group in df.groupby(STRATA, sort=True):
        donors = group[DONOR + ["diagnosis"]].drop_duplicates()
        n_case = int((donors.diagnosis == "case").sum())
        n_control = int((donors.diagnosis == "control").sum())
        gene_sizes = group.groupby("gene_id").size()
        if (gene_sizes != len(donors)).any():
            raise ValueError(f"Incomplete gene matrix in {study}/{phase}/{cell_type}; add explicit zero rows, not missing=0")
        for gene, rows in group.groupby("gene_id", sort=True):
            values = np.log2(1 + rows.counts.to_numpy(dtype=float) * 1_000_000 /
                             rows.library_size.to_numpy(dtype=float))
            cases = values[(rows.diagnosis == "case").to_numpy()]
            controls = values[(rows.diagnosis == "control").to_numpy()]
            effect = float(cases.mean() - controls.mean()) if len(cases) and len(controls) else float("nan")
            status, p = "ok", 1.0
            if min(n_case, n_control) < min_donors_per_group:
                status = "insufficient_donors"
            elif np.var(cases) == 0 and np.var(controls) == 0:
                status = "zero_variance"
            else:
                with np.errstate(all="ignore"):
                    p = float(ttest_ind(cases, controls, equal_var=False).pvalue)
                if not np.isfinite(p):
                    status, p = "undefined_test", 1.0
            results.append(dict(study_id=study, phase=phase, cell_type=cell_type,
                                gene_id=gene, n_case=n_case, n_control=n_control,
                                effect_log2_cpm=effect, p_value=p, status=status,
                                method="welch_log2_cpm_v0"))
    if not results:
        raise ValueError("No genes were available for analysis")
    result = pd.DataFrame(results)
    result["q_value"] = np.nan
    for _, indices in result.groupby(STRATA).groups.items():
        result.loc[indices, "q_value"] = bh_adjust(result.loc[indices, "p_value"].to_numpy())
    return result.sort_values(STRATA + ["gene_id"]).reset_index(drop=True)
