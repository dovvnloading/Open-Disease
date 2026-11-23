"""Check exact same-block concordance bounds against small synthetic orders."""

from collections import defaultdict
from itertools import permutations
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from matrisome_scan_bounds import feasible_within_wins


def enumerate_rankings(blocks):
    observations = [(case, block) for block, sizes in enumerate(blocks)
                    for case, size in zip((True, False), sizes) for _ in range(size)]
    expected = defaultdict(set)
    for order in set(permutations(observations)):
        wins = [(block, other_block) for i, (case, block) in enumerate(order) if case
                for other_case, other_block in order[:i] if not other_case]
        expected[len(wins)].add(sum(block == other for block, other in wins))
    return expected


class TestScanBounds(unittest.TestCase):
    def test_all_possible_totals_against_direct_ordering(self):
        for blocks in [[(1, 2), (2, 1)], [(2, 0), (1, 1), (0, 1)], [(0, 2), (2, 0)]]:
            for total, within in enumerate_rankings(blocks).items():
                with self.subTest(blocks=blocks, total=total):
                    self.assertEqual(feasible_within_wins(blocks, total), sorted(within))

    def test_rome_day_counts_and_total(self):
        self.assertEqual(feasible_within_wins([(3, 0), (3, 3), (2, 3)], 44), [11, 12, 13, 14, 15])

    def test_invalid_counts(self):
        for blocks, total in [([], 0), ([(1, -1)], 0), ([(1.5, 2)], 1),
                              ([(True, 1)], 0), ([(2, 2)], -1), ([(2, 2)], 5),
                              ([(2, 2)], 2.5)]:
            with self.subTest(blocks=blocks, total=total), self.assertRaises(ValueError):
                feasible_within_wins(blocks, total)
