"""Offline synthetic point-table fixtures for the 3p21.31 regional reader."""

import gzip
import io
from pathlib import Path
import struct
import sys
import tempfile
import unittest
from unittest.mock import patch
import zlib

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from gwas_3p21_disease_slice import extract, verify
from gwas_3p21_tabix import regional_rows


def bgzf_block(raw):
    comp = zlib.compressobj(wbits=-15)
    body = comp.compress(raw) + comp.flush()
    header = bytes.fromhex("1f8b08040000000000ff060042430200")
    size = 18 + len(body) + 8
    return header + struct.pack("<H", size - 1) + body + struct.pack("<II", zlib.crc32(raw), len(raw))


class RangeResponse(io.BytesIO):
    def __init__(self, request, payload):
        assert request.get_header("Accept-encoding") == "identity"
        first, desired_end = map(int, request.get_header("Range")[6:].split("-"))
        last = min(desired_end, len(payload) - 1)
        super().__init__(payload[first:last + 1])
        self.status = 206
        self.headers = {"Content-Range": f"bytes {first}-{last}/{len(payload)}"}


class TestRegionalReader(unittest.TestCase):
    def setUp(self):
        head = b"#chr\tpos\tvalue\n"
        first_data = b"3\t8\tx\n3\t9\ty\n"
        second_data = b"3\t10\tz\n3\t11\tw\n"
        block = bgzf_block(head + first_data)
        self.payload = block + bgzf_block(second_data)
        # The chunk ends inside block two. Include only the complete first row.
        chunk_end = (len(block) << 16) + len(b"3\t10\tz\n")
        self.raw_idx = (
            b"TBI\1" + struct.pack("<8i", 1, 0, 1, 2, 0, 35, 0, 2) + b"3\x00"
            + struct.pack("<iIiQQiQ", 1, 4681, 1, len(head), chunk_end, 1, len(head))
        )

    def query(self, idx=None, first=9, last=10, budget=200_000):
        if idx is None:
            idx = self.raw_idx
        with patch("urllib.request.urlopen", side_effect=lambda request, **_: RangeResponse(request, self.payload)):
            return regional_rows("https://test.invalid/test.gz", gzip.compress(idx), "3", first, last, budget)

    def test_virtual_offsets_and_inclusive_coordinates(self):
        h, rows, receipts = self.query()
        self.assertEqual(h, ["chr", "pos", "value"])
        self.assertEqual(rows, [["3", "9", "y"], ["3", "10", "z"]])
        self.assertEqual(len(receipts), 2)

    def test_zero_based_coordinates(self):
        idx = self.raw_idx[:8] + struct.pack("<i", 0x10000) + self.raw_idx[12:]
        self.assertEqual(self.query(idx, 10, 10)[1], [["3", "9", "y"]])

    def test_missing_bin_has_no_hits(self):
        _, rows, receipts = self.query(first=1_000_000, last=1_000_001)
        self.assertEqual(rows, [])
        self.assertEqual(len(receipts), 1)

    def test_budget_and_bad_index_fail(self):
        with self.assertRaisesRegex(ValueError, "budget"):
            self.query(budget=100)
        with self.assertRaisesRegex(ValueError, "Expected tabix"):
            self.query(b"nope")

    def test_rejects_unbounded_response(self):
        response = io.BytesIO(b"unrelated")
        response.status = 200
        response.headers = {}
        with patch("urllib.request.urlopen", return_value=response), self.assertRaisesRegex(ValueError, "did not honor"):
            regional_rows("https://test.invalid/test.gz", gzip.compress(self.raw_idx), "3", 9, 10, 200_000)


class TestDiseaseTranche(unittest.TestCase):
    def test_wrong_offline_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            directory = Path(name)
            (directory / "gwas_region.tsv").write_text("another release")
            with self.assertRaisesRegex(ValueError, "digest"):
                extract(directory, fetch=False)

    def test_digest_checks_exact_bytes(self):
        self.assertEqual(verify(b"abc", "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"), b"abc")
