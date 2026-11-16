"""Synthetic checks for the two coordinate conventions and missing marker handling."""

import tempfile
import unittest
from pathlib import Path

from scripts.gwas_3p21_rahmioglu_workbook import CANDIDATE, LEAD, analyze_uk, extract


def uk_row(variant, position, gene="ENSG00000186792", slope=0.5, column_l=0.125, flag=0):
    return (gene, 3, 0, 1, "-", 10, 0, variant, 3, position, position, column_l, slope, flag)


class WorkbookAuditTests(unittest.TestCase):
    def test_separate_coordinates_and_absent_is_unobserved(self):
        header = [()] * 3
        rows = header + [
            uk_row("3:49652148", 49614715),
            uk_row("rs1352889", 49614715, gene="ENSG00000126062", slope=-1, column_l=0.5),
            uk_row("3:49652148", 49614714, gene="OTHER"),
            uk_row("3:49200000", 49614715, gene="OTHER2"),
        ]
        result = analyze_uk([(34, rows), (35, header)])
        lead = result["marker_rows"][LEAD]
        self.assertEqual((lead["matched_id"], lead["matched_position"], lead["matched_both"]), (3, 3, 2))
        self.assertEqual(lead["distinct_gene_ids"], 4)
        self.assertEqual(lead["slope_squared_divided_by_L_range"], [2.0, 2.0])
        self.assertEqual(lead["prespecified_genes"]["TMEM115"][0]["row"], 5)
        candidate = result["marker_rows"][CANDIDATE]
        self.assertEqual(candidate["matched_both"], 0)
        self.assertIsNone(candidate["slope_squared_divided_by_L_range"])
        self.assertEqual(candidate["prespecified_genes"]["HYAL3"], [])
        self.assertEqual(sum(result["rows_in_audit_interval_by_numeric_variant_id"].values()), 3)

    def test_second_sheet_and_zero_column_l(self):
        rows = [()] * 3 + [uk_row("3:50332697", 50295266, column_l=0.0)]
        result = analyze_uk([(35, rows)])["marker_rows"][CANDIDATE]
        self.assertEqual(result["matched_both"], 1)
        self.assertEqual(result["prespecified_genes"]["HYAL3"][0]["sheet"], 35)
        self.assertIsNone(result["slope_squared_divided_by_L_range"])

    def test_wrong_workbook_fails_closed(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "different.xlsx"
            path.write_bytes(b"unmatched version")
            with self.assertRaisesRegex(ValueError, "verified 2023"):
                extract(path)

    def test_flags_are_kept_within_gene_across_sheets(self):
        h = [()] * 3
        s34 = h + [uk_row("3:49652148", 49614715, column_l=0.01)]
        s35 = h + [uk_row("3:60000000", 60000000, column_l=0.2, flag=1)]
        result = analyze_uk([(34, s34), (35, s35)])["all_cis_column_L_checks"]["HYAL3"]
        self.assertEqual(result["n_cis_rows"], 2)
        self.assertEqual(result["column_L_min"], 0.01)
        self.assertEqual(result["column_L_max"], 0.2)
        self.assertEqual(result["flagged_best"], [{"sheet": 35, "row": 4, "id": "3:60000000", "column_L": 0.2}])
