import unittest
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


path = Path(__file__).resolve().parents[1] / "scripts" / "ovarian_histotype_contrast.py"
spec = spec_from_file_location("ovarian_histotype_contrast", path)
runner = module_from_spec(spec)
spec.loader.exec_module(runner)


class RegistryContrastTest(unittest.TestCase):
    def test_printed_precision_and_nested_denominators(self):
        self.assertEqual(runner.rounding_bounds("1.90"), (1.895, 1.905))
        self.assertEqual(runner.rounding_bounds("1.9"), (1.85, 1.95))
        # Expected counts are 20 in total, 2 clear-cell, 3 endometrioid.
        # Observed counts are 40 in total, 4 clear-cell, 9 endometrioid.
        result = runner.standardized_contrast(
            {
                "country": "fixture",
                "total": "40",
                "total_sir": "2.0",
                "clear_cell": "4",
                "clear_cell_sir": "2.0",
                "endometrioid": "9",
                "endometrioid_sir": "3.0",
            }
        )
        self.assertAlmostEqual(result["ratio_of_sirs"], (13 / 5) / (27 / 15))

    def test_dutch_totals_and_interval_stratification(self):
        rows = runner.read_table("dutch_counts.tsv")
        result = runner.dutch_contrast(rows, "endometriosis")
        self.assertEqual(
            result["rows_exposure_nevus_cols_target_complement"], [[30, 39], [41, 200]]
        )
        self.assertAlmostEqual(
            result["follow_up_stratified_mantel_haenszel_odds_ratio"], 3.690871, places=6
        )
        self.assertAlmostEqual(result["conditional_case_mix_odds_ratio"], 3.732911, places=6)
