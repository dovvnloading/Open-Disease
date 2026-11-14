import gzip
from hashlib import sha256
from pathlib import Path
import struct
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import bulk_gse6364 as source


def synthetic_prefix(date, scanner=10000001):
    text = (f"Cols=1164\nRows=1164\nDatHeader=[20..65500] synthetic: "
            f"{date} 10:15:20 {scanner}  M10  \x14 HG-U133_Plus_2.1sq\n").encode("ascii")
    prefix = struct.pack("<6i", 64, 4, 1164, 1164, 1164**2, len(text)) + text
    return gzip.compress(prefix, mtime=0).ljust(2048, b"\0")


def test_scanner_audit_conditions_on_scan_day_and_keeps_ordered_file_pin(tmp_path, monkeypatch):
    groups = [True, False, True, False]
    dates = ["10/01/06", "10/01/06", "10/02/06", "10/03/06"]
    selected = []
    parts = []
    for i, (date, case) in enumerate(zip(dates, groups)):
        part = synthetic_prefix(date)
        (tmp_path / f"synthetic{i}.header.gzpart").write_bytes(part)
        parts.append(part)
        selected.append((i, f"synthetic{i}", case, "unused"))
    monkeypatch.setattr(source, "PREFIX_SHA", sha256(b"".join(parts)).hexdigest())
    summary = source.scan_audit(tmp_path, selected)
    assert summary["native_scanners"] == 1
    assert (summary["scan_days"], summary["mixed_scan_days"], summary["conditional_assignments"]) == (3, 1, 2)
    with pytest.raises(ValueError, match="digest"):
        source.scan_audit(tmp_path, list(reversed(selected)))


def test_bad_cel_geometry_and_truncated_range_fail_closed():
    short = synthetic_prefix("10/01/06")[:1000]
    with pytest.raises(ValueError, match="2,048"):
        source.scanner_fields(short)
    bad = gzip.compress(struct.pack("<6i", 64, 4, 20, 1164, 1164**2, 5) + b"dummy").ljust(2048, b"\0")
    with pytest.raises(ValueError, match="XDA"):
        source.scanner_fields(bad)


def test_native_image_match_requires_dat_and_all_four_corners_across_formats():
    from audit_us_bulk_sources import native_signature

    def binary_header(lines):
        payload = struct.pack("<6i", 64, 4, 1164, 1164, 1164**2, len(lines)) + lines
        return gzip.compress(payload, mtime=0).ljust(2048, b"\0")

    fields = [b"DatHeader=shared synthetic scan", b"GridCornerUL=0 0", b"GridCornerUR=1 0",
              b"GridCornerLR=1 1", b"GridCornerLL=0 1"]
    xda = native_signature(binary_header(b"\n".join(fields)))
    text = b"[CEL]\r\n[HEADER]\r\n" + b"\r\n".join(fields) + b"\r\n[INTENSITY]\r\n"
    other = native_signature(gzip.compress(text, mtime=3).ljust(2048, b"\0"))
    assert xda[0] == "XDA_v4" and other[0] == "text"
    assert xda[1] == other[1]
    fields[-1] = b"GridCornerLL=0 2"
    assert native_signature(binary_header(b"\n".join(fields)))[1] != xda[1]
    with pytest.raises(ValueError, match="Missing"):
        native_signature(binary_header(b"\n".join(fields[:-1])))
