"""Preflight audit for donor-level pseudobulk tables.

This audits structural suitability for an exploratory software baseline.
Passing does NOT establish unconfounded biology, cohort independence,
external clinical validity, or suitability for clinical decisions.
"""
from __future__ import annotations

import pandas as pd

from .schema import DONOR, validate_pseudobulk


def cohort_audit(data: pd.DataFrame, *, min_donors_per_group: int = 3,
                 min_cells_per_sample: int = 20, min_genes: int = 20) -> pd.DataFrame:
    """Summarize donor counts and block structurally unsuitable strata.

    The complete gene-matrix check detects missing donor/gene rows: absent
    rows may not silently be interpreted as biological zeros. The script
    cannot verify that an upstream export used *all* transcriptome features.
    """
    if min_donors_per_group < 3 or min_cells_per_sample < 1 or min_genes < 1:
        raise ValueError("min_donors_per_group>=3, min_cells_per_sample>=1, min_genes>=1 required")
    df = validate_pseudobulk(data)
    output = []
    for (study, phase, cell_type), group in df.groupby(
            ["study_id", "phase", "cell_type"], sort=True):
        donors = group[DONOR + ["diagnosis", "n_cells", "library_size"]].drop_duplicates()
        if donors.duplicated(DONOR).any():
            raise ValueError("Conflicting donor sample values within analysis stratum")
        n_case = int((donors["diagnosis"] == "case").sum())
        n_control = int((donors["diagnosis"] == "control").sum())
        n_total = len(donors)
        gene_n = group.groupby("gene_id")["donor_id"].nunique()
        n_genes = len(gene_n)
        matrix_complete = bool(gene_n.eq(n_total).all())
        min_cells = int(donors["n_cells"].min())
        min_library = int(donors["library_size"].min())
        flags = []
        if n_case < min_donors_per_group or n_control < min_donors_per_group:
            flags.append("insufficient_donors")
        if not matrix_complete:
            flags.append("incomplete_gene_matrix")
        if min_cells < min_cells_per_sample:
            flags.append("low_cell_count")
        if n_genes < min_genes:
            flags.append("low_gene_coverage")
        # Cannot infer independence, covariate balance, licensing, or consent
        # merely from donor IDs or an expression table.
        output.append({
            "study_id": study, "phase": phase, "cell_type": cell_type,
            "n_case_donors": n_case, "n_control_donors": n_control,
            "n_total_donors": n_total, "n_genes": n_genes,
            "min_cells_per_sample": min_cells, "min_library_size": min_library,
            "gene_matrix_complete": matrix_complete,
            "baseline_eligible": not flags, "blocking_flags": "|".join(flags),
            "external_independence_review": "required",
            "covariate_review": "required",
            "data_rights_review": "required",
        })
    return pd.DataFrame(output).sort_values(
        ["study_id", "phase", "cell_type"]).reset_index(drop=True)
