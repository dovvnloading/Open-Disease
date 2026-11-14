"""Offline synthetic census, with deliberately nonnumeric quantitative fields."""
from pathlib import Path
import sys
import zlib

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_endometdb_chip_census import census, pg_table_blocks, raw_keys


def rows(source):
    header = b"header"
    compressed = zlib.compress(source)
    payload = header + b"\0" + len(compressed).to_bytes(4, "little") + compressed + b"\0" * 5
    frames = pg_table_blocks((payload[n:n + 2] for n in range(0, len(payload), 2)), header)
    return list(raw_keys(frames))


def test_pg_frames_streaming_and_opaque_quantitative_fields():
    data = b"0\t1\t4\t?\t?\t?\t?\n1\t1\t5\tunknown\n\\.\n\n"
    assert rows(data) == [(b"1", b"4"), (b"1", b"5")]
    with pytest.raises(ValueError, match="Incomplete"):
        rows(b"0\t1\t4\t?\n")
    with pytest.raises(ValueError, match="after COPY"):
        rows(data + b"0\t1\t4\t?\n")
    with pytest.raises(ValueError, match="terminator"):
        list(pg_table_blocks([b"head"], b"head"))
    with pytest.raises(ValueError, match="beyond"):
        list(pg_table_blocks([b"head" + b"\0" * 6], b"head"))


def test_repeated_partial_missing_and_unmatched_assays_never_enter_overlap():
    selected = [{"analysis_result_id": str(i), "analysis_run_id": "10", "phase": "secretory", "group": group}
                for i, group in enumerate(["case", "control", "case", "control", "case"])]
    manifest = {"WG6": ({4, 5}, 4), "HT12": ({4, 6}, 6)}
    raw = [(b"0", b"4"), (b"1", b"5"), (b"1", b"4"), (b"0", b"5"),
           (b"2", b"5"), (b"2", b"4"), (b"2", b"5"), (b"3", b"4")]
    ids = {str(i).encode() for i in range(5)}
    status = {}
    report = census(raw, selected, ids, manifest, selected_status=status)
    assert status[b"0"] == "WG6" and status[b"4"] == "missing"
    assert report["selected_by_status"] == {"WG6": 2, "repeated_addresses": 1, "no_exact_manifest": 1, "missing": 1}
    assert report["selected_excess_address_rows"] == 1
    assert report["overlap_by_complete_manifest"]["WG6"] == {"case": 1, "control": 1, "shared_cells": 1}
    with pytest.raises(ValueError, match="outside"):
        census(raw + [(b"X", b"4")], selected, ids, manifest)


def test_identical_signatures_report_ambiguous_chip():
    selected = [{"analysis_result_id": "a", "analysis_run_id": "3", "phase": "p", "group": "case"}]
    report = census([(b"a", b"6")], selected, {b"a"}, {"A": ({6}, 6), "B": ({6}, 6)})
    assert report["selected_by_status"] == {"ambiguous_manifest": 1}


def test_default_source_accepts_documented_single_prefix(tmp_path, monkeypatch):
    import hashlib
    import audit_endometdb_chip_census as reader

    content = b"synthetic archive bytes"
    (tmp_path / "raw-prefix.zip.part").write_bytes(content)
    monkeypatch.setattr(reader, "ZIP_PREFIX_BYTES", len(content))
    monkeypatch.setattr(reader, "PREFIX_SHA256", hashlib.sha256(content).hexdigest())
    assert b"".join(reader.input_bytes(tmp_path, None)) == content
