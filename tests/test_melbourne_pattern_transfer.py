"""Regression checks for correlated genes, orientation and post-hoc phase contrasts."""
from itertools import combinations
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import numpy as np
import pytest

from cross_cohort_discovery import assignment_design
from finnish_cycle_alignment import evaluate, include_mixed, partial_r, within_strata
from melbourne_pattern_transfer import exact_correlations, standardize


def test_gram_method_against_independent_direct_enumeration_and_orientation():
    rng = np.random.default_rng(123)
    scores = rng.normal(size=(7, 20))
    scores[:, 5:10] += scores[:, :5] * 4
    ref = rng.normal(size=20)
    labels = np.arange(7) < 3
    _, weights, _ = assignment_design(labels, np.zeros(7))
    z = standardize(scores)
    calculated = exact_correlations(z, ref, weights)
    direct = []
    for combo in combinations(range(7), 3):
        mask = np.isin(np.arange(7), combo)
        contrast = z[mask].mean(axis=0) - z[~mask].mean(axis=0)
        direct.append(np.corrcoef(ref, contrast)[0, 1])
    np.testing.assert_allclose(calculated, direct, atol=1e-12)
    np.testing.assert_allclose(calculated, -exact_correlations(z, -ref, weights), atol=1e-12)


def test_stratified_null_scale_matches_full_distribution():
    rng = np.random.default_rng(28)
    scores = rng.normal(size=(8, 9))
    labels = np.array([1, 1, 0, 1, 0, 1, 0, 0], dtype=bool)
    strata = np.array([1, 1, 1, 2, 2, 3, 3, 3])
    ref = np.arange(9)
    observed, weights, factors = assignment_design(labels, strata)
    z = within_strata(scores, factors)
    direct = weights @ scores
    np.testing.assert_allclose(weights @ z, direct / np.std(direct, axis=0), atol=1e-12)
    row, stat = evaluate(scores, labels, strata, ref)
    expected = np.array([np.corrcoef(ref, t)[0, 1] for t in direct / np.std(direct, axis=0)])
    actual = np.corrcoef(ref, observed @ z)[0, 1]
    assert row["assignments"] == 18
    assert row["extreme"] == sum(np.abs(expected) >= abs(actual) - 1e-12)
    assert np.isclose(row["r"], actual)
    np.testing.assert_allclose(stat, observed @ z)


def test_positive_control_eligibility_and_residualization():
    labels = np.array([1, 1, 1, 0, 0], dtype=bool)
    np.testing.assert_array_equal(include_mixed(labels, np.array([0, 0, 1, 1, 2])),
                                  [False, False, True, True, False])
    axis = np.arange(10)
    x = np.array([1, 5, 6, 8, 0, 4, 8, 4, 2, 9])
    assert np.isclose(partial_r(x + 2 * axis, -x + 33 * axis, axis), -1)


def test_zero_information_rejected():
    with pytest.raises(ValueError, match="zero variation"):
        standardize(np.zeros((5, 10)))
    with pytest.raises(ValueError, match="Constant external"):
        exact_correlations(np.eye(5), np.ones(5), np.eye(5))
