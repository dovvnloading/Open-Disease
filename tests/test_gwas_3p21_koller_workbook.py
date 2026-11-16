"""Check table-specific selection against small in-memory worksheets."""

import tempfile
import unittest
from pathlib import Path

from openpyxl import Workbook

from scripts.gwas_3p21_koller_workbook import END, START, audit, extract_table


def sheet(number, rows):
    ws = Workbook().active
    for _ in range(5 if number == 4 else 3):
        ws.append(["header"])
    for row in rows:
        ws.append(row)
    return ws


class KollerWorkbookTests(unittest.TestCase):
    def test_lead_alleles_effect_offsets_and_section_heading(self):
        row = ["rs1", 3, START, "T", "C", "nearest", "known", *range(7, 52)]
        other_chrom = list(row)
        other_chrom[1] = 4
        result = extract_table(sheet(4, [row, ["section title"], other_chrom]), 4)
        self.assertEqual(result["data_rows_scanned"], 2)
        observed, = result["rows"]
        self.assertEqual(observed["row"], 6)
        self.assertEqual(observed["effect_allele"], "T")
        effects = observed["effects"]
        self.assertEqual(effects["cross_ancestry_combined"]["beta"], 7)
        self.assertEqual(effects["eur_combined"]["beta"], 27)
        self.assertEqual(effects["eur_clinical"]["eaf"], 35)

    def test_secondary_signals_fill_down_analysis_and_filter_at_boundaries(self):
        rows = [
            ["EUR clinical", 3, START - 1, "rs_outside", "A", "T", 0, 1, 1, 0, 1, 1],
            [None, 3, END, "rs_edge", "B", "C", 0.3, 1, 1, 0.2, 1, 1],
        ]
        data = extract_table(sheet(5, rows), 5)
        self.assertEqual(data["data_rows_scanned"], 2)
        observed, = data["rows"]
        self.assertEqual(observed["variant"], "rs_edge")
        self.assertEqual(observed["analysis"], "EUR clinical")
        self.assertEqual(observed["joint_beta"], 0.2)

    def test_locus_overlap_even_when_paintor_marker_outside(self):
        rows = [
            ["EUR", "chr3", START - 3, "rs_disjoint", "A", "G", 2, 0.7, f"chr3:1-{START - 1}"],
            ["Cross", "chr3", END + 10, "rs_touch", "C", "T", 4, 0.8, f"chr3:{END}-{END + 15}"],
            [None, "chr3", START, "rs_internal", "G", "A", 5, 0.9, f"chr3:{START}-{START + 1}"],
        ]
        data = extract_table(sheet(7, rows), 7)
        self.assertEqual(data["data_rows_scanned"], 3)
        first, second = data["rows"]
        self.assertTrue(first["locus_overlaps"])
        self.assertFalse(first["marker_in_interval"])
        self.assertTrue(second["marker_in_interval"])
        self.assertEqual(second["analysis"], "Cross")
        rows[0][-1] = "chr4:1-10"
        with self.assertRaisesRegex(ValueError, "chromosome"):
            extract_table(sheet(7, rows), 7)

    def test_coloc_keeps_all_chromosome_three_genes(self):
        row = ["EUR clinical", "Uterus", "ENSG", "other gene", 3, 0.1, 0.7, 0.9, "rs2", 0.6]
        second = list(row)
        second[4] = 4
        data = extract_table(sheet(8, [row, second]), 8)
        self.assertEqual(data["data_rows_scanned"], 2)
        hit, = data["rows"]
        self.assertEqual((hit["gene_symbol"], hit["regional_probability"]), ("other gene", 0.9))
        self.assertEqual(hit["p_gwas_as_labelled"], 0.1)

    def test_workbook_pin_checked_before_parsing(self):
        with tempfile.TemporaryDirectory() as folder:
            file = Path(folder) / "not_the_final_workbook.xlsx"
            file.write_bytes(b"changed")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                audit(file)


if __name__ == "__main__":
    unittest.main()
