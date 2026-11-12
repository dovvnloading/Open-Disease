import pandas as pd
import pytest

from open_disease.qc import cohort_audit
from open_disease.synthetic import generate_synthetic


def synthetic(tmp_path):
    path, _ = generate_synthetic(tmp_path)
    return pd.read_csv(path)


def test_synthetic_cohorts_are_structurally_eligible_not_clinically_valid(tmp_path):
    report = cohort_audit(synthetic(tmp_path))
    assert len(report) == 4
    assert report.baseline_eligible.all()
    assert set(report.n_case_donors) == {6}
    assert set(report.n_control_donors) == {6}
    assert set(report.external_independence_review) == {"required"}
    assert set(report.covariate_review) == {"required"}


def test_missing_gene_is_not_interpreted_as_zero(tmp_path):
    df = synthetic(tmp_path).drop(index=0)
    report = cohort_audit(df)
    assert not report.baseline_eligible.all()
    assert report.blocking_flags.str.contains("incomplete_gene_matrix").any()


def test_weak_sample_size_is_blocked(tmp_path):
    report = cohort_audit(synthetic(tmp_path), min_donors_per_group=7)
    assert not report.baseline_eligible.any()
    assert report.blocking_flags.str.contains("insufficient_donors").all()


def test_low_gene_coverage_is_blocked(tmp_path):
    report = cohort_audit(synthetic(tmp_path), min_genes=60)
    assert not report.baseline_eligible.any()
    assert report.blocking_flags.str.contains("low_gene_coverage").all()


def test_low_cell_count_is_blocked(tmp_path):
    report = cohort_audit(synthetic(tmp_path), min_cells_per_sample=100)
    assert not report.baseline_eligible.any()
    assert report.blocking_flags.str.contains("low_cell_count").all()


def test_bad_threshold_is_rejected(tmp_path):
    with pytest.raises(ValueError, match="min_donors_per_group"):
        cohort_audit(synthetic(tmp_path), min_donors_per_group=1)
