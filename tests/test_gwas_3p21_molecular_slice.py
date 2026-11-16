"""Test bounded molecular downloads, digest validation and local re-use."""

import hashlib
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import gwas_3p21_molecular_slice as molecular


class MolecularSliceTests(unittest.TestCase):
    def test_fetch_missing_reference_then_verify_offline(self):
        content = b"ACTG"
        manifest = {"regional_matrices": {}, "small_sources": {
            "GRCh38.fa.txt": {"url": "https://example.org/ref", "accept": "text/plain",
                             "max_bytes": 32, "sha256": hashlib.sha256(content).hexdigest()}
        }}
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with patch.object(molecular, "small_request", return_value=content + b"\n") as fetch:
                molecular.acquire(root, manifest, True)
                molecular.acquire(root, manifest, False)
                self.assertEqual(fetch.call_count, 1)
            (root / "GRCh38.fa.txt").write_bytes(b"AAAA")
            with self.assertRaisesRegex(ValueError, "SHA-256"):
                molecular.acquire(root, manifest, False)

    def test_response_limit(self):
        response = MagicMock()
        response.__enter__.return_value.read.return_value = b"1234"
        with patch.object(molecular, "urlopen", return_value=response):
            with self.assertRaisesRegex(ValueError, "byte limit"):
                molecular.small_request("https://example.org/data", 3, "text/plain")

    def test_no_unverified_partial_file(self):
        idx = b"index"
        source = {"url": "https://example.org/T.tsv.gz", "tbi_bytes": len(idx),
                  "tbi_sha256": hashlib.sha256(idx).hexdigest(), "tsv_sha256": "wrong",
                  "ranges": [], "rows": 1}
        manifest = {"regional_matrices": {"T": source}, "grch38_envelope": [5, 7]}
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(molecular, "small_request", return_value=idx), patch.object(
                molecular, "regional_rows", return_value=(["column"], [["value"]], [])
            ):
                with self.assertRaisesRegex(ValueError, "SHA-256"):
                    molecular.acquire(Path(tmp), manifest, True)
            self.assertEqual(list(Path(tmp).iterdir()), [])
