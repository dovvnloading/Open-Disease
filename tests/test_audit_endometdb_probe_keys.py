"""Invented assay aliases, probe addresses and opaque quantitative fields."""
from pathlib import Path
import sys
import zlib

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_endometdb_probe_keys import copy_prefix, normalized_summary, raw_summary


def test_prefix_ignores_partial_lines_and_rejects_wrong_table():
    tag = bytes.fromhex("0100da080000")
    encoded = zlib.compress(b"row1\nrow2\npartial")
    block = tag + b"\0" + len(encoded).to_bytes(4, "little") + encoded
    assert copy_prefix(block, 0, tag) == [b"row1", b"row2"]
    with pytest.raises(ValueError, match="header"):
        copy_prefix(block, 0, b"other!")


def test_only_full_unique_manufacturer_address_set_verifies_chip():
    lines = [b"id\tA\t10\t?", b"id\tA\t20\t?",
             b"id\tB\t10\t?", b"id\tB\t20\t?", b"id\tB\t20\t?", b"id\tC\t10\t?"]
    output = raw_summary(lines, {10: "probe1", 20: "probe2"}, {10}, {"A", "B", "C"}, {"B", "C"})
    assert output["complete_wg6_assays"] == 1
    assert output["complete_wg6_eligible_assays"] == 0
    assert output["wg6_addresses_absent_from_ht12"] == 1
    assert output["raw_rows_in_prefix"] == 6
    with pytest.raises(ValueError, match="native assay"):
        raw_summary(lines, {}, set(), {"B", "C"}, set())


def test_duplicate_nonmeasurement_keys_survive_signal_masking():
    rows = [b"a\t?\t0\tGENE\tprobe\t0", b"a\t!?\t0\tGENE\tprobe\t0",
            b"b\t?\t0\tGENE\tprobe\t0", b"a\t?\tkey\tFOSB\tfosb_probe\t42"]
    summary = normalized_summary(rows)
    assert summary["distinct_keys_excluding_signal"] == 3
    assert summary["zero_search_key_and_address_rows"] == 3
    assert summary["first_probe_rows_per_assay_range"] == [1, 2]
    assert summary["fosb_symbol_rows"] == 1
