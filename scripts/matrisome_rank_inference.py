"""Exact permutation tests for the matrisome benchmark's fixed score ordering."""

from collections import Counter
from math import isfinite
from numbers import Integral

from cross_cohort_discovery import holm


def _rank_counts(cases, comparators):
    """Count binary label orderings by their number of case-over-comparator wins."""
    # Each column holds the distribution for a given number of comparators.
    # An ordering ending with the highest-score comparator gives no new wins.
    # One ending with the highest-score case adds a win over every comparator.
    row = [[1] for _ in range(comparators + 1)]
    for _ in range(cases):
        for n in range(1, comparators + 1):
            highest_comparator = row[n - 1]
            counts = [0] * n + row[n]
            for wins, frequency in enumerate(highest_comparator):
                counts[wins] += frequency
            row[n] = counts
    return row[-1]


def _integer(value):
    if not isinstance(value, Integral) or isinstance(value, bool):
        raise ValueError("Counts, wins and weights must be integers")
    return int(value)


def ranked_sum_distribution(blocks, weights=None):
    """Count all within-block label assignments by weighted U in a strict ordering."""
    blocks = list(blocks)
    weights = [1] * len(blocks) if weights is None else list(weights)
    if not blocks or len(blocks) != len(weights):
        raise ValueError("Provide at least one stratum and one weight per stratum")

    joint = {0: 1}
    for (cases, comparators), weight in zip(blocks, weights, strict=True):
        cases, comparators, weight = map(_integer, (cases, comparators, weight))
        if min(cases, comparators, weight) < 1:
            raise ValueError("Each stratum needs both groups and a positive weight")
        combined = Counter()
        for u, frequency in enumerate(_rank_counts(cases, comparators)):
            for previous, previous_frequency in joint.items():
                combined[previous + weight * u] += previous_frequency * frequency
        joint = combined
    return joint


def exact_test_from_distribution(null, observed):
    """Calculate the inclusive two-sided tail around the exact permutation mean."""
    observed = _integer(observed)
    if observed not in null:
        raise ValueError("The observed statistic is outside the permutation support")
    assignments = sum(null.values())
    # Scale by the assignment count so a fractional mean never needs rounding.
    center = sum(value * n for value, n in null.items())
    distance = abs(assignments * observed - center)
    extreme = sum(n for value, n in null.items() if abs(assignments * value - center) >= distance)
    return {"observed": observed, "maximum": max(null), "extreme": extreme,
            "assignments": assignments, "two_sided_p": extreme / assignments}


def exact_rank_test(strata, weights=None):
    """Test a weighted sum of within-stratum U statistics in a strict rank order.

    Each triple gives case count, comparator count, and observed case wins.
    Ranks must be fixed independently of labels, including any tie-break.
    Positive integer weights multiply each U before summing.
    """
    strata = list(strata)
    weights = [1] * len(strata) if weights is None else list(weights)
    blocks = []
    observed = 0
    for (cases, comparators, wins), weight in zip(strata, weights, strict=True):
        cases, comparators, wins, weight = map(_integer, (cases, comparators, wins, weight))
        if not 0 <= wins <= cases * comparators:
            raise ValueError("The number of wins must lie between zero and the number of pairs")
        observed += weight * wins
        blocks.append((cases, comparators))
    return exact_test_from_distribution(ranked_sum_distribution(blocks, weights), observed)


def holm_three_targets(p_values):
    """Adjust the three primary target comparisons, preserving their names."""
    if len(p_values) != 3 or not all(isfinite(p) and 0 <= p <= 1 for p in p_values.values()):
        raise ValueError("The primary family requires three finite target p-values in [0, 1]")
    return dict(zip(p_values, map(float, holm(list(p_values.values()))), strict=True))
