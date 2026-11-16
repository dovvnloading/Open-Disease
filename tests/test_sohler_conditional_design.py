from datetime import datetime, timedelta
from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from sohler_conditional_design import TISSUES, conditional_design, verify_dat_filename


def test_filename_binding_uses_exact_dat_and_version_delimiters():
    prefix = bytes(24) + b"DatHeader=[1..65535] example:CLS=8100 RWS=8100 \x14 chip\n"
    verify_dat_filename(prefix, "Example.v4.CEL")
    for filename in ("Other.v4.CEL", "Example.CEL", "Example.v3.CEL", "path/Example.v4.CEL"):
        with pytest.raises(ValueError):
            verify_dat_filename(prefix, filename)
    for wrong_header in (prefix.replace(b":CLS", b"-CLS"), prefix.replace(b"\x14", b" ")):
        with pytest.raises(ValueError):
            verify_dat_filename(wrong_header, "Example.v4.CEL")


def test_crossed_days_supply_a_cycle_but_one_way_edge_is_fixed():
    l_tissue, p_tissue, _ = TISSUES
    start = datetime(2007, 1, 1, 9)
    # One within-day pair, a balanced two-day cycle and one unreciprocated edge.
    hours = [(1, 2), (3, 24), (25, 4), (28, 48)]
    rows = [(i, tissue, start+timedelta(hours=t), "scanner")
            for i, (left, right) in enumerate(hours)
            for tissue, t in ((l_tissue, left), (p_tissue, right))]
    design = conditional_design(rows)
    assert design["unconditional_assignments"] == 16
    assert design["scanner_day_margin_preserving_assignments"] == 4
    assert design["movable_pairs"] == 3
    assert design["movable_cross_scanner_day_pairs"] == 2
    assert not design["all_pairs_swapped_preserves_scanner_days"]
    assert design["same_scanner_lesion_scanned_first"] == 3
    assert design["additionally_preserving_same_scanner_lesion_first_count"] == {
        "assignments": 2, "movable_pairs": 2, "movable_cross_scanner_day_pairs": 2,
    }
    assert design["mixed_scanner_days_by_number_of_tissue_switches_among_paired_arrays"] == {1: 1, 3: 1}


def test_separate_scanners_on_one_date_form_distinct_blocks():
    l_tissue, p_tissue, uterus = TISSUES
    time = datetime(2007, 1, 1, 9)
    rows = [("a", l_tissue, time, "s1"), ("a", p_tissue, time, "s2"),
            ("b", p_tissue, time, "s3"), ("b", uterus, time, "s4")]
    design = conditional_design(rows)
    assert design["paired_donors"] == design["different_scanner_order_unassigned"] == 1
    assert design["scanner_day_margin_preserving_assignments"] == 1
    assert design["movable_pairs"] == 0
