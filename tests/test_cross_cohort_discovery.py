"""Regression tests for stratified maxT, held-out direction and nonredundant nomination."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import cross_cohort_discovery as screen


def test_synchronous_maxT_including_ties_and_constant_columns():
    case = np.array([True, False, True, False])
    strata = np.array(["a", "a", "b", "b"])
    # Two equal-size paired strata. The two columns have different directions
    # in the second stratum. Some assignment maximizes one of them every time.
    x = np.array([[3, 3, 7], [1, 1, 7], [6, 1, 8], [4, 3, 8]], float)
    design = screen.assignment_design(case, strata)
    e, z, raw, mx, variable = screen.test_family(x, design, chunk_size=1)
    np.testing.assert_array_equal(e, [2, 0, 0])
    np.testing.assert_allclose(z, [np.sqrt(2), 0, 0], atol=1e-14)
    np.testing.assert_array_equal(raw, [2, 4, 4])
    np.testing.assert_array_equal(mx, [4, 4, 4])
    np.testing.assert_array_equal(variable, [True, True, False])
    p = screen.test_family(x, design, directions=[1, -1, 1])[2]
    np.testing.assert_array_equal(p, [1, 3, 4])


def test_stratified_variance_matches_empirical_randomization_null():
    cases = np.array([True, False, True, True, False, False, False])
    strata = np.array(["a"] * 3 + ["b"] * 4)
    x = np.array([[9,1], [3,2], [4,7], [6,9], [5,3], [5,6], [1,1]], float)
    design = screen.assignment_design(cases, strata)
    _, weights, _ = design
    obs = (x[cases & (strata == "a")].mean(0) - x[~cases & (strata == "a")].mean(0)
           + x[cases & (strata == "b")].mean(0) - x[~cases & (strata == "b")].mean(0)) / 2
    all_effects = weights @ x
    scale = all_effects.std(0)
    z = obs / scale
    effect, statistic, raw, adjusted, _ = screen.test_family(x, design)
    np.testing.assert_allclose(effect, obs)
    np.testing.assert_allclose(statistic, z)
    np.testing.assert_array_equal(raw, (abs(all_effects/scale) >= abs(z) - screen.TOL).sum(0))
    maxima = abs(all_effects/scale).max(1)
    np.testing.assert_array_equal(adjusted, (maxima[:, None] >= abs(z) - screen.TOL).sum(0))


def test_holm_direction_and_pathway_overlap_filter():
    np.testing.assert_allclose(screen.holm([.03, .005, .04, .03]), [.09, .02, .09, .09])
    rows = [dict(id=f"R-HSA-{i}", statistic=z, effect_pp=np.sign(z)) for i,z in [(9, -7), (1, 7), (4, 3), (3, 0)]]
    members = {"R-HSA-1": ("first", [1,2]), "R-HSA-9": ("skip", [1,2,3,4]),
               "R-HSA-4": ("second", [4,5])}
    assert [r["id"] for r in screen.nominate(rows, members)] == ["R-HSA-1", "R-HSA-4"]
    assert screen.nominate([rows[0]], None)[0]["direction"] == "lower"
    with pytest.raises(ValueError, match="lacks cases or controls"):
        screen.assignment_design(np.array([True, False]), np.array(["a", "b"]))


def test_modified_lock_is_rejected_before_any_holdout_read(tmp_path):
    changed = tmp_path / "nominations.json"
    changed.write_text('{"candidates":[]}')
    with pytest.raises(ValueError, match="public nomination"):
        screen.validation(tmp_path / "missing", changed, tmp_path)
