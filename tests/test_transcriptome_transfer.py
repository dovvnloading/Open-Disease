"""Regression tests for the transcriptome-wide cohort and matched-date comparisons."""

from pathlib import Path
import struct
import sys
import zlib

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from rome_diagnostics import blocked_effect
from rome_scan import scanner_fields
from transcriptome_transfer import evaluate_projection


def test_donor_allocations_control_dependence_across_genes():
    case = np.array([True] * 3 + [False] * 3)
    scores = np.array([[0.0, 1.0], [1.0, 3.0], [2.0, 2.0], [5.0, 0.0], [4.0, 2.0], [6.0, 1.0]])
    w = np.array([2.0, -1.0])
    result = evaluate_projection(scores, case, w)
    # Duplicating each gene cannot create more independent units.
    repeat = evaluate_projection(np.tile(scores, (1, 2)), case, np.tile(w, 2))
    assert result["score_difference"] < 0
    assert result["assignments"] == repeat["assignments"] == 20
    assert result["p_two_sided"] == repeat["p_two_sided"] == 0.1
    assert result["score_difference"] == pytest.approx(repeat["score_difference"])


def test_exclude_pure_scan_date_and_equal_weight_other_dates():
    observations = np.array([0.0, 4.0, 7.0, 2.0, 1.0, 90.0, 91.0])
    case = np.array([True, False, True, False, False, True, True])
    dates = np.array([1, 1, 2, 2, 2, 3, 3])
    r = blocked_effect(observations, case, dates)
    assert (r["cases"], r["controls"], r["assignments"]) == (2, 3, 6)
    assert r["within_day_effects"] == [-4.0, 5.5]
    assert r["effect"] == pytest.approx(0.75)
    with pytest.raises(ValueError, match="No mixed"):
        blocked_effect(observations, case, case.astype(int))


def test_scan_parser_recovers_native_date_and_rejects_other_chip():
    def prefix(design):
        metadata = f"DatHeader= 02/09/10 13:00:09 12345678 M10  {design}\n".encode()
        head = struct.pack("<6i", 64, 4, 732, 732, 732**2, len(metadata))
        compressor = zlib.compressobj(wbits=31)
        return (compressor.compress(head + metadata) + compressor.flush()).ljust(2048, b"\0")

    date, instrument = scanner_fields(prefix("HG-U133A_2.1sq"))
    assert date.isoformat() == "2010-02-09" and instrument == "12345678"
    with pytest.raises(ValueError, match="chip design"):
        scanner_fields(prefix("OTHER"))
