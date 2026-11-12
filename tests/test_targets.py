import pandas as pd
import pytest

from open_disease.analysis import analyze_pseudobulk
from open_disease.synthetic import generate_synthetic
from open_disease.targets import rank_targets, validate_targets


def inputs(tmp_path):
    pb, path = generate_synthetic(tmp_path)
    results = analyze_pseudobulk(pd.read_csv(pb))
    evidence = pd.read_csv(path)
    ranked = rank_targets(results, evidence, discovery_study="synthetic-discovery",
                          phase="proliferative", cell_type="epithelial")
    return results, ranked, evidence


def test_rank_only_mapped_genes_from_discovery(tmp_path):
    _, ranked, _ = inputs(tmp_path)
    assert "SYN000" in set(ranked.gene_id)
    assert "SYN003" not in set(ranked.gene_id)
    assert (ranked.evidence_type != "locus_only").all()
    assert ranked["rank"].tolist() == list(range(1, len(ranked) + 1))


def test_holdout_replication_of_programmed_signal(tmp_path):
    results, ranked, _ = inputs(tmp_path)
    validation = validate_targets(ranked, results, validation_study="synthetic-holdout")
    outcome = validation[validation.gene_id == "SYN000"].iloc[0]
    assert bool(outcome.replicated)
    assert outcome.validation_status == "replicated_exploratory"


def test_cannot_reuse_discovery_as_holdout(tmp_path):
    results, ranked, _ = inputs(tmp_path)
    with pytest.raises(ValueError, match="DATA LEAKAGE"):
        validate_targets(ranked, results, validation_study="synthetic-discovery")


def test_support_score_outside_unit_interval_rejected(tmp_path):
    results, _, evidence = inputs(tmp_path)
    evidence.loc[0, "support_score"] = 5.0
    with pytest.raises(ValueError, match="support_score"):
        rank_targets(results, evidence, discovery_study="synthetic-discovery",
                     phase="proliferative", cell_type="epithelial")


def test_synthetic_generation_is_reproducible(tmp_path):
    generate_synthetic(tmp_path / "a", seed=48)
    generate_synthetic(tmp_path / "b", seed=48)
    assert (tmp_path / "a/pseudobulk.csv").read_bytes() == (tmp_path / "b/pseudobulk.csv").read_bytes()
