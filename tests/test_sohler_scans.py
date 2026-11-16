from datetime import datetime, timedelta
from pathlib import Path
import struct
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_sohler_scans import HEADER_BYTES, TISSUES, scan, scan_design


def header(chip=b"HG-U133_Plus_2.1sq", header_size=700):
    fixed = struct.pack("<6i", 64, 4, 1164, 1164, 1164**2, header_size)
    dat = b"Cols=1164\nDatHeader=test 01/02/07 10:11:12 1234  M10   " + chip + b" \n"
    return (fixed + dat).ljust(HEADER_BYTES, b" ")


def test_bounded_native_format():
    assert scan(header()) == (datetime(2007, 1, 2, 10, 11, 12), b"1234")
    for wrong in (header()[:-1], header(b"other"), header(header_size=250)):
        with pytest.raises(ValueError):
            scan(wrong)


def test_mixed_blocks_are_distinct_from_within_donor_crossing():
    start = datetime(2007, 1, 1)
    minute = timedelta(minutes=1)
    day = timedelta(days=1)
    lesion, peritoneum, uterus = TISSUES
    records = [("a", lesion, start, "s"), ("a", peritoneum, start+minute, "s"),
               ("b", peritoneum, start+2*minute, "s"), ("b", lesion, start+3*minute, "s"),
               ("c", lesion, start+4*minute, "s"), ("c", peritoneum, start+day, "s"),
               ("d", lesion, start+5*minute, "s"), ("d", peritoneum, start+6*minute, "t"),
               ("e", uterus, start+7*minute, "t")]
    summary = scan_design(records)
    assert summary["paired_donors"] == 4 and summary["donors"] == 5
    assert summary["scanner_days_with_lesion_and_peritoneum"] == 1
    assert summary["same_scanner_pairs"] == 3
    assert summary["same_scanner_day_pairs"] == 2
    assert summary["order_same_scanner_day"] == {"lesion_first": 1, "peritoneum_first": 1, "equal_clock_time": 0}
    with pytest.raises(ValueError, match="Repeated"):
        scan_design(records + [("f", uterus, start, "s")])
    with pytest.raises(ValueError, match="Repeated"):
        scan_design(records + [("d", lesion, start+day, "t")])
