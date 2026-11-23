from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from matrisome_omission import omission_summary


def possible_subsets(cases, comparators):
    """Independently score every ordering and every subset that retains both arms."""
    size = cases + comparators
    by_wins = {}
    for case_indices in combinations(range(size), cases):
        labels = [i in case_indices for i in range(size)]
        pairs = [(i, j) for i in range(size) if labels[i]
                 for j in range(size) if not labels[j]]
        higher = sum(i > j for i, j in pairs)
        outcomes = by_wins.setdefault(higher, set())
        for kept_count in range(2, size + 1):
            for kept in combinations(range(size), kept_count):
                retained = [(i, j) for i, j in pairs if i in kept and j in kept]
                if retained:
                    outcomes.add((size - kept_count, sum(i > j for i, j in retained), len(retained)))
    return by_wins


def read_bounds(bounds):
    if bounds is None:
        return None
    return tuple(Fraction(**bounds[key]) for key in ("lower", "upper"))


class OmissionTests(unittest.TestCase):
    def check_enumeration(self, strata, possible):
        pooled = {}
        for outcomes in product(*possible):
            omitted, higher, pairs = map(sum, zip(*outcomes))
            pooled.setdefault(omitted, set()).add(Fraction(higher, pairs))
        baseline, = pooled[0]
        half = Fraction(1, 2)
        if baseline == half:
            threshold = 0
        elif baseline > half:
            threshold = min((k for k, values in pooled.items() if min(values) <= half), default=None)
        else:
            threshold = min((k for k, values in pooled.items() if max(values) >= half), default=None)
        result = omission_summary(strata)
        self.assertEqual(Fraction(**result["baseline"]["fraction"]), baseline)
        self.assertEqual(result["minimum_omissions_to_lose_direction"], threshold)
        for field, k in (("one_omission", 1), ("at_threshold", threshold),
                         ("just_before_threshold", threshold - 1 if threshold else None)):
            expected = (min(pooled[k]), max(pooled[k])) if k in pooled else None
            self.assertEqual(read_bounds(result[field]), expected)
            if expected:
                self.assertEqual(result[field]["omissions"], k)

    def test_all_small_single_stratum_rankings(self):
        for cases, comparators in ((1, 1), (1, 3), (3, 2), (3, 3)):
            for higher, outcomes in possible_subsets(cases, comparators).items():
                with self.subTest(stratum=(cases, comparators, higher)):
                    self.check_enumeration([(cases, comparators, higher)], [outcomes])

    def test_every_pair_of_small_strata_with_different_arm_sizes(self):
        first, second = possible_subsets(2, 2), possible_subsets(1, 3)
        for h1, h2 in product(first, second):
            self.check_enumeration([(2, 2, h1), (1, 3, h2)], [first[h1], second[h2]])

    def test_unequal_phase_weights(self):
        result = omission_summary([(6, 3, 5), (9, 8, 56), (6, 5, 16)])
        self.assertEqual(result["direction"], "above_half")
        self.assertEqual(result["minimum_omissions_to_lose_direction"], 4)
        self.assertEqual(read_bounds(result["one_omission"]), (Fraction(68, 111), Fraction(77, 111)))
        self.assertEqual(read_bounds(result["just_before_threshold"])[0], Fraction(50, 93))
        self.assertEqual(read_bounds(result["at_threshold"])[0], Fraction(41, 84))

    def test_invalid_aggregate_counts(self):
        for strata in ([], [(0, 3, 0)], [(3, -1, 0)], [(2, 3, 7)], [(2, 2, -1)],
                       [(3.0, 2, 1)], [(2, 1, True)]):
            with self.subTest(strata=strata), self.assertRaises(ValueError):
                omission_summary(strata)
