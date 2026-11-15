"""Regression coverage for synchronized profile projection and cohort inclusion."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from cross_cohort_discovery import assignment_design
from finland_rome_transport import score_transport, overlap
from verify_genome_transport import direct_dot


def test_projection_equals_explicit_whole_donor_relabeling():
    cases = np.array([True, False, True, False, True])
    x = np.array([[2, 5, 1], [1, 2, 7], [7, 3, 4], [2, 1, 5], [9, 4, 3]], dtype=float)
    reference = np.array([1.4, .7, -.8])
    design = assignment_design(cases, ['all'] * len(cases))
    result, z = score_transport(x, reference, cases, design)
    expected, extreme, explicit_z = direct_dot(x, cases, reference)
    assert result['mean_product'] == pytest.approx(expected)
    assert result['absolute_extreme'] == extreme
    np.testing.assert_allclose(z, explicit_z)
    # A duplicated full gene panel must give identical inference, even when
    # every gene has a perfectly correlated duplicate.
    twice, _ = score_transport(np.tile(x, (1, 2)), np.tile(reference, 2), cases, design)
    assert twice['exact_two_sided_p'] == result['exact_two_sided_p']
    assert twice['mean_product'] == pytest.approx(expected)


def test_reference_reversal_and_zero_variance_guard():
    cases = np.array([True, True, False, False])
    x = np.array([[1., 9.], [2., 8.], [8., 2.], [9., 1.]])
    ref = np.array([1., -1.])
    design = assignment_design(cases, ['a'] * 4)
    a, _ = score_transport(x, ref, cases, design)
    b, _ = score_transport(x, -ref, cases, design)
    assert a['mean_product'] == pytest.approx(-b['mean_product'])
    assert a['exact_two_sided_p'] == b['exact_two_sided_p']
    assert a['pairs_with_higher_case_score'] == 0
    assert b['pairs_with_higher_case_score'] == 4
    with pytest.raises(ValueError, match='constant'):
        score_transport(np.column_stack((x[:, 0], np.ones(len(x)))), ref, cases, design)
    paired = assignment_design(cases, ['a', 'b', 'a', 'b'])
    with pytest.raises(ValueError, match='single'):
        score_transport(x, ref, cases, paired)


def test_top_overlap_resolves_equal_statistic_by_entrez():
    ids = np.arange(1, 101)[::-1]
    ref = np.ones(len(ids))
    target = np.where(ids == 1, -10., 1.)
    row = overlap(ids, ref, target)[0]
    assert (row['n'], row['shared_in_top'], row['shared_opposite_direction']) == (1, 1, 1)


def test_scan_date_weighting_excludes_unmixed_date_and_preserves_within_date_counts():
    from rome_scan_transport import compare_mixed_dates
    cases = np.array([True, False, True, True, False, True, True])
    dates = np.array(['a', 'a', 'b', 'b', 'b', 'c', 'c'])
    profiles = np.array([[1,4], [3,2], [1,8], [2,7], [5,4], [6,2], [7,3]], dtype=float)
    result = compare_mixed_dates(profiles, np.array([1.,-1.]), cases, dates)
    assert result['assignments'] == 2 * 3
    assert result['excluded_arrays'] == 2
    assert result['pairs_case_higher'] == 0
    assert result['case_control_pairs'] == 3
    assert result['mean_product'] == pytest.approx(np.mean([r['mean_product'] for r in result['date_specific']]))
