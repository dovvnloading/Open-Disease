"""Synthetic checks for allele alignment, exclusions and pre-intersection peak loss."""

from collections import Counter
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from gwas_3p21_sensitivity import START, UNMAPPED, canonical_snv, disease_table, evaluate_gene, joint_row


def disease_row(pos, rsid, effect="G", other="A", p=0.01):
    return dict(chromosome="3", base_pair_location=str(pos), effect_allele=effect,
                other_allele=other, rsid=rsid, beta="0.2", standard_error="0.1",
                effect_allele_frequency="0.4", p_value=str(p))


def qtl_row(pos, rsid, ref="A", alt="G", p=0.01):
    return dict(chromosome="3", position=str(pos), rsid=rsid, ref=ref, alt=alt,
                pvalue=str(p), beta="0.3", se="0.1", maf="0.35", an="258",
                variant=f"chr3_{pos}_{ref}_{alt}", gene_id="ENSG00000186792")


class TestSensitivityGates(unittest.TestCase):
    def test_effect_allele_flipped_to_reference_checked_alt(self):
        rows = [disease_row(START, "one", effect="A", other="G"),
                disease_row(START + 1, "two", effect="T", other="C")]
        disease, markers, counts = disease_table(rows, "AC")
        self.assertEqual(disease[markers["one"]], ("one", -0.2))
        self.assertEqual(disease[markers["two"]], ("two", 0.2))
        self.assertEqual(counts, {"eligible": 2})

    def test_unresolved_and_duplicate_sites(self):
        rows = [disease_row(START, "pal", "A", "T"),
                disease_row(START + 1, "duplicate1"), disease_row(START + 1, "duplicate2")]
        disease, _, counts = disease_table(rows, "AA")
        self.assertFalse(disease)
        self.assertEqual(counts, {"unresolved_palindrome": 1, "duplicate_or_multiallelic_position": 2})
        self.assertEqual(canonical_snv(next(iter(UNMAPPED)), "A", "G", ""), (None, "unmapped"))

    def test_reference_mismatch_and_duplicate_gene_sites(self):
        position = Counter([START])
        r = qtl_row(START, "one", ref="G", alt="A")
        self.assertEqual(joint_row(r, "A", position, {}), ("qtl_ref_mismatch", None))
        r = qtl_row(START, "one")
        self.assertEqual(joint_row(r, "A", Counter({START: 2}), {}), ("duplicate_or_multiallelic_gene_position", None))

    def test_missing_disease_peak_stops_otherwise_complete_gene(self):
        rows = [disease_row(START, "gene-peak"), disease_row(START + 1, "disease-peak")]
        disease, markers, _ = disease_table(rows, "AA")
        match = evaluate_gene([qtl_row(START, "gene-peak")], "AA", disease, markers,
                              ["gene-peak", "disease-peak"], {"gene-peak", "disease-peak"})
        self.assertTrue(match["all_tied_qtl_peaks_retained"])
        self.assertFalse(match["eligible"])
        self.assertEqual(match["credible_set_retained_by_disease_rsid"], 1)
        self.assertEqual(match["missing_required_markers"], ["disease-peak"])

    def test_peak_loss_checked_before_gwas_intersection(self):
        disease, markers, _ = disease_table([disease_row(START, "disease-peak")], "AA")
        qtl = [qtl_row(START, "disease-peak", p=1e-5), qtl_row(START + 1, "lost", p=1e-8)]
        match = evaluate_gene(qtl, "AA", disease, markers, ["disease-peak"], set())
        self.assertFalse(match["eligible"])
        self.assertFalse(match["all_tied_qtl_peaks_retained"])
        self.assertFalse(match["missing_required_markers"])

    def test_complete_joint_coverage_and_alleles(self):
        disease, markers, _ = disease_table([disease_row(START, "peak")], "A")
        match = evaluate_gene([qtl_row(START, "peak")], "A", disease, markers, ["peak"], {"peak"})
        self.assertTrue(match["eligible"])
        self.assertEqual(match["qc_rows"], {"retained": 1})
        self.assertEqual(match["focus_marker_effects"]["peak"]["alt_effect_allele"], "G")
