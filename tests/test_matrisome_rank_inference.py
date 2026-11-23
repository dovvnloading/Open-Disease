from collections import Counter
from itertools import combinations, product
from math import comb
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from matrisome_rank_inference import (
    exact_rank_test, exact_test_from_distribution, holm_three_targets, ranked_sum_distribution,
)


def enumerate_wins(cases, comparators):
    """Enumerate whole-array assignments and count all cross-group comparisons directly."""
    ranks = set(range(cases + comparators))
    wins = []
    for label_assignment in combinations(sorted(ranks), cases):
        controls = ranks.difference(label_assignment)
        wins.append(sum(case > control for case in label_assignment for control in controls))
    return wins


@pytest.mark.parametrize("weights", [(1, 1), (3, 2)])
def test_joint_null_matches_label_enumeration_at_every_observation(weights):
    dimensions = [(2, 2), (1, 3)]
    support = [enumerate_wins(*group) for group in dimensions]
    null = [sum(w * u for w, u in zip(weights, scores)) for scores in product(*support)]
    assert ranked_sum_distribution(dimensions, weights) == Counter(null)
    for win_counts in product(*(sorted(set(values)) for values in support)):
        strata = [(*group, wins) for group, wins in zip(dimensions, win_counts)]
        observed = sum(w * u for w, u in zip(weights, win_counts))
        mean = sum(null) / len(null)
        extreme = sum(abs(value - mean) >= abs(observed - mean) for value in null)
        assert exact_rank_test(strata, weights) == {
            "observed": observed, "maximum": max(null), "extreme": extreme,
            "assignments": len(null), "two_sided_p": extreme / len(null),
        }


def test_out_of_support_weighted_sum_is_rejected():
    even_scores = ranked_sum_distribution([(3, 3)], [2])
    with pytest.raises(ValueError, match="outside"):
        exact_test_from_distribution(even_scores, 1)


def test_many_joint_assignments_do_not_overflow_machine_integers():
    count = comb(10, 5) ** 10
    test = exact_rank_test([(5, 5, 0)] * 10)
    assert test["assignments"] == count
    assert test["extreme"] == 2
    assert test["two_sided_p"] == 2 / count


@pytest.mark.parametrize(("strata", "weights"), [
    ([], None), ([(1, 1, 1)], []), ([(1, 0, 0)], None),
    ([(1, 1, 2)], None), ([(1, 1, -1)], None),
    ([(1, 1, 1)], [0]), ([(1, 1, 0.5)], None), ([(True, 1, 0)], None),
])
def test_invalid_inputs(strata, weights):
    with pytest.raises(ValueError):
        exact_rank_test(strata, weights)


def test_holm_preserves_names_and_enforces_three_targets():
    adjusted = holm_three_targets({"Mumbai": .09, "US": .08, "Rome": .03})
    assert adjusted == pytest.approx({"US": .16, "Rome": .09, "Mumbai": .16})
    with pytest.raises(ValueError):
        holm_three_targets({"US": .05, "Mumbai": .1})
    with pytest.raises(ValueError):
        holm_three_targets({"US": .1, "Rome": float("nan"), "Mumbai": .11})
