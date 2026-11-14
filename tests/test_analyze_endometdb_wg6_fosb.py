"""Exercise the Finnish result with invented donors, address sets and signals."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_endometdb_wg6_fosb as fosb


def source():
    design = {b"case": (("run", "phase"), "case"), b"control": (("run", "phase"), "control")}
    lines = []
    # Target address 11. Tied negative signals give a doubled midrank of 2.
    for donor, values in [(b"case", [-2, -2, 0, 7]), (b"control", [3, 2, -1, 2])]:
        for address, value in zip([10, 11, 12, 13], values):
            lines.append(b"0\t" + donor + b"\t" + str(address).encode() + b"\t"
                         + str(value).encode() + b"\t?\t?\t?")
    lines.insert(2, b"0\texcluded\tthis record has opaque signals")
    return lines, design


def test_tied_negative_and_zero_finite_signals_are_ranked():
    lines, design = source()
    cells = fosb.rank_signals(lines[::-1], design, {10, 11, 12, 13}, 11)
    assert cells == {("run", "phase"): {"case": [2], "control": [4]}}


@pytest.mark.parametrize(("change", "error"), [
    ("null", "Missing or malformed"), ("malformed", "Missing or malformed"),
    ("nan", "Nonfinite"), ("infinity", "Nonfinite"), ("duplicate", "duplicate"),
    ("missing", "Incomplete"), ("extra", "Unexpected"), ("short", "Malformed selected"),
])
def test_compare_requires_finite_signal_at_every_selected_address(change, error):
    lines, design = source()
    if change in ("null", "malformed", "nan", "infinity"):
        bad = {"null": b"\\N", "malformed": b"arbitrary", "nan": b"NaN", "infinity": b"-inf"}[change]
        lines[0] = lines[0].replace(b"\t-2\t", b"\t" + bad + b"\t")
    elif change == "duplicate":
        lines += lines[:1]
    elif change == "missing":
        lines.pop(0)
    elif change == "short":
        lines[0] = b"0\tcase\t10\t22"
    else:
        lines.append(b"0\tcase\t14\t0\t?\t?\t?")
    with pytest.raises(ValueError, match=error):
        fosb.rank_signals(lines, design, {10, 11, 12, 13}, 11)


def test_exact_statistic_keeps_equal_strata_and_all_tied_assignments():
    # First stratum has three possible choices of case: contrasts 1/2, 1/2, -1.
    # Second has two choices: observed +1 and alternative -1.
    # Absolute sum >= observed 3/2 is reached by 2 observed-like assignments + 1 opposite.
    cells = {"a": {"case": [3], "control": [3, 2]}, "b": {"case": [2], "control": [1]}}
    result = fosb.exact_contrast(cells, 4)
    assert result["difference_percentile_points"] == pytest.approx(9.375)
    assert result["total_assignments"] == 6 and result["extreme_assignments"] == 3
    assert result["two_sided_p_fraction"] == "1/2" and result["decision"] == "inconclusive"
    assert result["leave_one_stratum_out_range_percentile_points"] == [6.25, 12.5]


@pytest.mark.parametrize(("a", "b", "decision"), [
    (3, 1, "supports_expected_direction"), (1, 3, "supports_opposite_direction"),
    (3, 3, "inconclusive"),
])
def test_two_sided_decision_and_exact_zero(a, b, decision):
    result = fosb.exact_contrast({i: {"case": [a], "control": [b]} for i in range(6)}, 2)
    assert result["decision"] == decision
    assert result["total_assignments"] == 64
    assert result["extreme_assignments"] == (64 if a == b else 2)


def test_no_shared_design_raises():
    with pytest.raises(ValueError, match="two shared"):
        fosb.exact_contrast({"a": {"case": [3], "control": [2]}}, 4)
    with pytest.raises(ValueError, match="both disease"):
        fosb.exact_contrast({"a": {"case": [3], "control": []}, "b": {"case": [3], "control": [2]}}, 4)


def frozen_donors():
    donors, status = [], {}

    def append(cell, group, chip):
        idx = str(len(donors))
        donors.append({"donor_id": idx, "analysis_result_id": idx, "analysis_run_id": cell[0],
                       "phase": cell[1], "group": group})
        status[idx.encode()] = chip
    for cell, (cases, controls) in fosb.EXPECTED_CELLS.items():
        for group in ["case"] * cases + ["control"] * controls:
            append(cell, group, "WG6")
    for group, run, n in (("case", "20", 3), ("control", "21", 4)):
        for _ in range(n):
            append((run, "other"), group, "WG6")
    for _ in range(6):
        append(("22", "other"), "case", "repeated_addresses")
    return donors, status


def test_design_enforces_unique_donors_verified_chip_and_frozen_overlap():
    donors, status = frozen_donors()
    design = fosb.shared_design(donors, status)
    assert len(design) == 25 and all(status[key] == "WG6" for key in design)
    donors[0]["phase"] = "other"
    with pytest.raises(ValueError, match="Shared run"):
        fosb.shared_design(donors, status)
    donors[0]["donor_id"] = donors[1]["donor_id"]
    with pytest.raises(ValueError, match="donor"):
        fosb.shared_design(donors, status)
