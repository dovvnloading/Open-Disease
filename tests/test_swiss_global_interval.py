"""Post-outcome compatibility bounds on invented donor scores only."""
from fractions import Fraction
from itertools import combinations
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from swiss_global_interval import allocation_contrasts, shift_limits, TOLERANCE


def test_inversion_matches_direct_reallocation_at_every_transition():
    outcomes = list(map(Fraction, [0, 1, 3, 5, 8]))
    labels = [1, 1, 0, 0, 0]
    batches = [0] * 5
    observed, contrasts = allocation_contrasts(outcomes, labels, batches)
    lower, upper, required = shift_limits(observed, contrasts, alpha=Fraction(3, 10))
    assert required == 4
    assert lower < observed < upper
    # Independently derived crossings of the two signed equality branches,
    # then direct labeling on both sides of every crossing, including boundaries.
    roots = {observed, lower, upper}
    for intercept, slope in contrasts:
        for perm_sign in (-1, 1):
            for obs_sign in (-1, 1):
                denominator = perm_sign*slope - obs_sign
                if denominator:
                    roots.add((perm_sign*intercept - obs_sign*observed + TOLERANCE) / denominator)
    ordered = sorted(roots)
    points = roots | {(a+b)/2 for a, b in zip(ordered, ordered[1:])} | {min(roots)-1, max(roots)+1}
    for theta in points:
        z = [v - d*theta for v, d in zip(outcomes, labels)]
        def difference(case):
            return sum(z[i] for i in case)/2 - sum(z[i] for i in range(5) if i not in case)/3
        extreme = sum(abs(difference(case))+TOLERANCE >= abs(difference((0, 1)))
                      for case in combinations(range(5), 2))
        assert (extreme >= required) == (lower <= theta <= upper)


def test_balanced_two_donors_cannot_exclude_any_shift():
    observed, contrasts = allocation_contrasts([Fraction(4), Fraction(7)], [1, 0], [1, 1])
    assert shift_limits(observed, contrasts) == (None, None, 1)
