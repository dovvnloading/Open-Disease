from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_heca_cell_libraries import SourceRanges, coverage


def test_sparse_file_reads_adjacent_verified_ranges_and_refuses_gaps():
    stream = SourceRanges(22, {0: b"abcd", 4: b"ef", 20: b"xy"})
    stream.seek(2)
    data = bytearray(4)
    assert stream.readinto(data) == 4 and data == b"cdef"
    with pytest.raises(ValueError, match="outside"):
        stream.readinto(bytearray(1))
    stream.seek(-1, 2)
    assert stream.readinto(data) == 1 and data[0:1] == b"y"


def test_cell_cutoff_applies_to_each_member_before_mixing():
    counts = {("a", "p"): 30, ("a", "q"): 15, ("b", "p"): 21, ("b", "q"): 40, ("c", "r"): 50}
    groups = {"p": True, "q": False, "r": True}
    low = coverage(counts, groups, 1)
    assert (low["mixed_libraries"], low["distinct_cases_in_mixed"], low["gland_cells_in_mixed"]) == (2, 1, 106)
    high = coverage(counts, groups, 20)
    assert (high["mixed_libraries"], high["distinct_controls_in_mixed"], high["gland_cells_in_mixed"]) == (1, 1, 61)
