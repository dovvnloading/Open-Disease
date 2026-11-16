"""Check phased-count LD, allele orientation, and the two source receipts."""

import json
import unittest
from copy import deepcopy

from scripts.gwas_3p21_pairwise_ld import INPUT, calculate, verify


class PairwiseLDTests(unittest.TestCase):
    rows = ["C", "T"]
    columns = ["A", "G"]
    counts = {("C", "A"): 345, ("C", "G"): 494, ("T", "A"): 121, ("T", "G"): 46}

    def test_providers_and_signed_direction(self):
        result = verify(json.loads(INPUT.read_text(encoding="utf-8")))
        self.assertEqual(result["counted_alleles"], {"rs1352889": "T", "rs13100173": "A"})
        stats = result["ldpair_count_calculation"]
        self.assertEqual(stats["phased_chromosomes"], 1006)
        self.assertAlmostEqual(stats["signed_r"], 0.233816435489)
        self.assertAlmostEqual(stats["r_squared"], 0.054670125505)
        self.assertAlmostEqual(stats["signed_d_prime"], 0.48684852517)
        self.assertAlmostEqual(stats["pearson_chi_squared_df1_uncorrected"], 54.998146258)
        self.assertLess(stats["pearson_asymptotic_p_value"], 0.0001)
        self.assertEqual(result["ensembl_same_panel_crosscheck"]["r_squared"], 0.053964)

    def test_single_recoding_reverses_signed_statistics(self):
        original = calculate(self.counts, self.rows, self.columns, ("T", "A"))
        flipped = calculate(self.counts, self.rows, self.columns, ("T", "G"))
        for key in ("d", "signed_r", "signed_d_prime"):
            self.assertAlmostEqual(original[key], -flipped[key])
        self.assertEqual(original["r_squared"], flipped["r_squared"])

    def test_exact_independence_and_perfect_negative_ld(self):
        independent = dict.fromkeys(self.counts, 5)
        self.assertEqual(calculate(independent, self.rows, self.columns, ("T", "A"))["signed_r"], 0)
        opposite = {("C", "A"): 4, ("C", "G"): 0, ("T", "A"): 0, ("T", "G"): 3}
        result = calculate(opposite, self.rows, self.columns, ("T", "A"))
        self.assertEqual((result["signed_r"], result["r_squared"], result["signed_d_prime"]), (-1, 1, -1))

    def test_reject_invalid_counts(self):
        for counts in (
            {},
            dict.fromkeys(self.counts, 0),
            {**self.counts, ("C", "A"): -1},
            {**self.counts, ("C", "A"): True},
            {**self.counts, ("C", "A"): 1.5},
            {**self.counts, ("C", "A"): 0, ("C", "G"): 0},
        ):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                calculate(counts, self.rows, self.columns, ("T", "A"))

    def test_reject_changed_receipts_and_marker_order(self):
        original = json.loads(INPUT.read_text(encoding="utf-8"))
        for section in ("ldpair", "ensembl"):
            source = deepcopy(original)
            receipt = source if section == "ldpair" else source["ensembl"]
            receipt["receipt_text"] += " "
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                verify(source)
        original["marker_order"].reverse()
        with self.assertRaisesRegex(ValueError, "order"):
            verify(original)


if __name__ == "__main__":
    unittest.main()
