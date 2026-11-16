"""Offline checks of regional coordinates and published-list bookkeeping."""

import hashlib
import json
import unittest
from copy import deepcopy

from scripts.gwas_3p21_colocalization_audit import GRCH37_CHR3_LENGTH, verify


class LocusAuditTests(unittest.TestCase):
    def setUp(self):
        self.archive = json.dumps({
            "selected": [
                {"symbol": "SHARED", "band": "3p21.31"},
                {"symbol": "TWAS", "band": "3p21.31"},
                {"symbol": "OTHER", "band": "3p21.3"},
            ],
            "measured_bands": ["IRRELEVANT"],
        }).encode()
        self.observations = {
            "disease_lead": {
                "id": "rs1352889", "chromosome": "3", "assembly": "GRCh37", "position": 49652148,
            },
            "coloc_candidate": {
                "id": "rs13100173", "chromosome": "3", "assembly": "GRCh37", "position": 50332697,
            },
            "published_radius_bp": 500000,
            "audit_radius_bp": 1000000,
            "reported_endometrial_eqtl_genes": ["SHARED", "EQTL"],
            "uterus_twas_source_sha256": hashlib.sha256(self.archive).hexdigest(),
        }

    def test_coordinate_coverage(self):
        result = verify(self.observations, self.archive)
        self.assertEqual(result["physical_distance_bp"], 680549)
        self.assertEqual(result["windows"]["published"], {
            "radius_bp": 500000, "start_inclusive": 49152148, "end_inclusive": 50152148,
            "contains_candidate": False, "distance_outside_bp": 180549,
        })
        self.assertEqual(result["windows"]["audit"], {
            "radius_bp": 1000000, "start_inclusive": 48652148, "end_inclusive": 50652148,
            "contains_candidate": True, "distance_outside_bp": 0,
        })

    def test_closed_interval(self):
        for offset, contained, excess in (
            (-500001, False, 1), (-500000, True, 0), (500000, True, 0), (500001, False, 1),
        ):
            with self.subTest(offset=offset):
                self.observations["coloc_candidate"]["position"] = 49652148 + offset
                window = verify(self.observations, self.archive)["windows"]["published"]
                self.assertEqual(window["contains_candidate"], contained)
                self.assertEqual(window["distance_outside_bp"], excess)

    def test_gene_sets_keep_analysis_and_band_separate(self):
        result = verify(self.observations, self.archive)
        self.assertEqual(result["published_gene_annotations"], {
            "reported_mtag_endometrial_eqtl": ["EQTL", "SHARED"],
            "selected_koller_uterus_twas": ["SHARED", "TWAS"],
            "both": ["SHARED"],
            "mtag_endometrial_eqtl_only": ["EQTL"],
            "koller_uterus_twas_only": ["TWAS"],
        })

    def test_changed_archive_is_rejected_even_if_json_equivalent(self):
        with self.assertRaisesRegex(ValueError, "SHA-256"):
            verify(self.observations, self.archive + b" ")

    def test_invalid_marker_or_radius(self):
        for marker in ("disease_lead", "coloc_candidate"):
            for key, value in (
                ("assembly", "GRCh38"), ("chromosome", "4"), ("position", 0),
                ("position", -1), ("position", 1.5), ("position", True),
                ("position", GRCH37_CHR3_LENGTH + 1),
            ):
                with self.subTest(marker=marker, key=key, value=value):
                    data = deepcopy(self.observations)
                    data[marker][key] = value
                    with self.assertRaises(ValueError):
                        verify(data, self.archive)
        for key in ("published_radius_bp", "audit_radius_bp"):
            for value in (0, -1, False, "100", 1.5):
                with self.subTest(key=key, value=value):
                    with self.assertRaises(ValueError):
                        verify({**self.observations, key: value}, self.archive)

    def test_boundaries_clip_to_chromosome(self):
        for pos in (1, GRCH37_CHR3_LENGTH):
            with self.subTest(position=pos):
                self.observations["disease_lead"]["position"] = pos
                self.observations["coloc_candidate"]["position"] = pos
                window = verify(self.observations, self.archive)["windows"]["audit"]
                self.assertGreaterEqual(window["start_inclusive"], 1)
                self.assertLessEqual(window["end_inclusive"], GRCH37_CHR3_LENGTH)
                self.assertTrue(window["contains_candidate"])


if __name__ == "__main__":
    unittest.main()
