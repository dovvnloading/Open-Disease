"""Transparent hypothesis rankings with strictly separated holdout validation.

Heuristic scores are NOT causal probabilities, patient predictions or treatment
recommendations. Holdout results are never allowed to select or rerank targets.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

EVIDENCE_TYPES = {"fine_mapping", "colocalization", "integrated_multiomics", "locus_only"}
EVIDENCE_COLS = ["gene_id", "evidence_type", "support_score", "source_id", "safety_review"]


def required(df: pd.DataFrame, cols: list[str]) -> None:
    missing = sorted(set(cols) - set(df.columns))
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))


def rank_targets(de: pd.DataFrame, evidence: pd.DataFrame, *, discovery_study: str,
                 phase: str, cell_type: str, fdr: float = 0.10, top_n: int = 25) -> pd.DataFrame:
    required(de, ["study_id", "phase", "cell_type", "gene_id", "status", "q_value",
                  "effect_log2_cpm", "n_case", "n_control"])
    required(evidence, EVIDENCE_COLS)
    if not (0 < fdr < 1) or top_n < 1:
        raise ValueError("fdr must be in (0,1) and top_n must be positive")
    if evidence.gene_id.duplicated().any():
        raise ValueError("Resolve duplicate evidence assessments per gene before ranking")
    if evidence[EVIDENCE_COLS].isna().any().any():
        raise ValueError("Genetic evidence cannot contain missing required fields")
    if not set(evidence.evidence_type) <= EVIDENCE_TYPES:
        raise ValueError("Unknown evidence_type")
    if not set(evidence.safety_review) <= {"not_assessed", "needs_review", "known_concern"}:
        raise ValueError("Invalid safety_review status")
    scores = pd.to_numeric(evidence.support_score, errors="coerce")
    if scores.isna().any() or not np.isfinite(scores.to_numpy(dtype=float)).all() or ((scores < 0) | (scores > 1)).any():
        raise ValueError("support_score must be a declared heuristic in [0,1], not a causal probability")
    if evidence.source_id.astype(str).str.strip().eq("").any():
        raise ValueError("Missing evidence source_id")
    subset = de[(de.study_id == discovery_study) & (de.phase == phase) &
                (de.cell_type == cell_type)].copy()
    if subset.empty:
        raise ValueError("Selected discovery study, phase and cell type are absent")
    if subset.gene_id.duplicated().any():
        raise ValueError("Duplicate discovery gene results")
    subset = subset[(subset.status == "ok") & (subset.q_value <= fdr) &
                    subset.effect_log2_cpm.notna()]
    joined = subset.merge(evidence, on="gene_id", how="inner", validate="one_to_one")
    joined = joined[joined.evidence_type != "locus_only"].copy()
    joined["discovery_study"] = discovery_study
    joined["score"] = (0.45 * joined.support_score.astype(float) +
                       0.35 * (joined.effect_log2_cpm.abs() / 2).clip(upper=1) +
                       0.20 * (1 - joined.q_value / fdr).clip(lower=0, upper=1) -
                       0.25 * (joined.safety_review == "known_concern").astype(float))
    joined = joined.sort_values(["score", "gene_id"], ascending=[False, True]).head(top_n)
    joined.insert(0, "rank", range(1, len(joined) + 1))
    columns = ["rank", "gene_id", "discovery_study", "phase", "cell_type", "score",
               "effect_log2_cpm", "q_value", "n_case", "n_control", "evidence_type",
               "support_score", "source_id", "safety_review"]
    return joined.reindex(columns=columns).reset_index(drop=True)


def validate_targets(ranked: pd.DataFrame, de: pd.DataFrame, *, validation_study: str,
                     fdr: float = 0.10) -> pd.DataFrame:
    required(ranked, ["rank", "gene_id", "discovery_study", "phase", "cell_type", "effect_log2_cpm"])
    required(de, ["study_id", "phase", "cell_type", "gene_id", "effect_log2_cpm",
                  "q_value", "status", "n_case", "n_control"])
    if not (0 < fdr < 1):
        raise ValueError("fdr must be in (0,1)")
    if (ranked.discovery_study == validation_study).any():
        raise ValueError("DATA LEAKAGE: validation_study equals discovery_study")
    if ranked.duplicated(["gene_id", "phase", "cell_type"]).any():
        raise ValueError("Duplicate ranked target")
    holdout = de[de.study_id == validation_study].copy()
    if holdout.duplicated(["gene_id", "phase", "cell_type"]).any():
        raise ValueError("Duplicate holdout result")
    holdout = holdout.rename(columns={
        "effect_log2_cpm": "validation_effect_log2_cpm", "q_value": "validation_q_value",
        "status": "holdout_test_status", "n_case": "validation_n_case",
        "n_control": "validation_n_control"
    })
    result = ranked.merge(holdout[["gene_id", "phase", "cell_type", "validation_effect_log2_cpm",
                                   "validation_q_value", "holdout_test_status",
                                   "validation_n_case", "validation_n_control"]],
                          on=["gene_id", "phase", "cell_type"], how="left", validate="one_to_one")
    result["replicated"] = ((result.holdout_test_status == "ok") &
                            (result.validation_q_value <= fdr) &
                            (np.sign(result.effect_log2_cpm) == np.sign(result.validation_effect_log2_cpm)) &
                            (result.effect_log2_cpm != 0))
    result["validation_status"] = np.where(
        result.holdout_test_status.isna(), "not_measured",
        np.where(result.replicated, "replicated_exploratory", "not_replicated"))
    result.insert(0, "validation_study", validation_study)
    return result.sort_values("rank").reset_index(drop=True)
