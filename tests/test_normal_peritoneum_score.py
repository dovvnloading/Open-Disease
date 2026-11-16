"""Synthetic donors, tissue pairs, and raw signals. No public study values."""
from fractions import Fraction
from itertools import product
from pathlib import Path
import sys

import numpy as np
import pytest
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import normal_peritoneum_score as score
from normal_peritoneum_design import Pair


def synthetic_pairs(n):
    return [Pair(f"toy-{i}", "batch", f"sP{i}", f"sE{i}", (f"P{i}",), (f"E{i}",),
                 "secretory", "stored_false") for i in range(n)]


def test_inverted_rational_limits_are_attainable_breakpoints():
    deltas = [Fraction(x, 7) for x in [-2, 0, 1, 2, 3, 9, 9, 12]]
    lo, hi = score.common_shift_interval(deltas)
    signs = list(product((-1, 1), repeat=len(deltas)))
    def independent(shift):
        centered = [x - shift for x in deltas]
        threshold = abs(sum(centered))
        return sum(abs(sum(x * s for x, s in zip(centered, pattern))) >= threshold
                   for pattern in signs) * 20 > len(signs)
    eps = Fraction(1, 10**6)
    assert [independent(x) for x in (lo-eps, lo, lo+eps, hi-eps, hi, hi+eps)] == [
        False, True, True, True, True, False]
    assert score.common_shift_interval([Fraction(2)] * 6) == (2, 2)
    assert score.tail_count([Fraction(2)] * 6) == 2
    assert score.tail_count([Fraction(0)] * 6) == 64
    assert score.common_shift_interval(deltas[:5]) == (None, None)


def test_underresolved_gate_stops_before_reading_signal_or_annotation():
    def fail_if_accessed():
        raise AssertionError("Signal read before eligibility")
        yield
    with pytest.raises(ValueError, match="Outcome gate"):
        score.score_qualified_pairs(synthetic_pairs(5), fail_if_accessed(), {})


def test_end_to_end_native_order_ties_repeats_and_frozen_addresses(monkeypatch):
    from dataclasses import replace

    pairs = synthetic_pairs(6)
    pairs[0] = replace(pairs[0], peritoneum_assays=("P0", "P0b"))
    pairs[-1] = replace(pairs[-1], medication="stored_true")
    annotation = {addr: gene for gene, (_, addr) in score.TARGETS.items()}
    annotation[3] = 0  # Measured reference-only feature
    monkeypatch.setattr(score, "N_FEATURES", len(annotation))
    addresses = sorted(annotation)
    target_order = [a for _, a in score.TARGETS.values()]
    values = {a: (1 if sign < 0 else 8) for sign, a in score.TARGETS.values()}
    values[3] = 8
    lines = [b"0\tnot-in-lock\tunusable-unselected-record"]
    for i, pair in enumerate(pairs):
        for assay in pair.peritoneum_assays + pair.uterus_assays:
            for a in addresses:
                v = values[a] if assay.startswith("P") and assay != "P0b" else 3
                lines.append(f"0\t{assay}\t{a}\t{v}\t.\t.\t.".encode())
    actual = score.score_qualified_pairs(pairs, reversed(lines), annotation)
    ranks = rankdata([values[a] for a in addresses])
    pp = dict(zip(addresses, 100 * (ranks - 0.5)/len(addresses)))
    signed = [s * (pp[a] - 50) for s, a in score.TARGETS.values()]
    full_delta = sum(signed) / 5
    assert actual["mean_difference_pp"] == pytest.approx(full_delta * 5.5 / 6)
    assert actual["phase_sensitivity"]["mean_difference_pp"] == pytest.approx(full_delta * 4.5 / 5)
    assert actual["two_sided_p"] == 2 / 64
    assert actual["technical_control"]["assay_pairs"] == 1
    assert actual["positive_zero_negative"] == [6, 0, 0]
    assert actual["verified_all_target_midranks"]
    assert actual["assays"] == 13
    assert actual["decision"] == "higher_in_peritoneum"
    assert not any(p.donor in str(actual) for p in pairs)
    bad = dict(annotation)
    bad[target_order[2]] = 0
    with pytest.raises(ValueError, match="physical target"):
        score.score_qualified_pairs(pairs, (), bad)


def test_frozen_plan_deterministic_controls():
    score.synthetic_controls()
    tied = np.ones(7)
    assert score.assay_score(tied, [0, 1, 2, 3, 4]) == 10
