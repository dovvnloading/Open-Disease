from itertools import product
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from adjacent_peritoneum import paired_shift_interval, ranks, summarize


def test_exact_interval_matches_sign_flip_inversion_at_boundaries():
    effects = np.array([-2, -1, -1, 1, 2, 5.0])
    low, high = paired_shift_interval(effects)
    flips = np.array(list(product([-1, 1], repeat=6)))
    for shift, accepted in [
        (low, True),
        (high, True),
        (low - 0.0001, False),
        (high + 0.0001, False),
    ]:
        centered = effects - shift
        pv = np.mean(abs(flips @ centered / 6) >= abs(centered.mean()) - 1e-12)
        assert (pv > 0.05) == accepted
    assert paired_shift_interval(effects[:5]) == [None, None]


def test_percentiles_include_zeros_and_all_measured_genes():
    assert ranks(np.array([0, 0, 2, 7]), [0, 2, 3]) == pytest.approx([25, 62.5, 87.5])
    diffs = np.array([[10, -10], [8, -8], [5, -5], [2, -2], [3, -3], [1, -1.0]])
    out = summarize(diffs, np.array([1, -1]), True)
    assert out["assignments"] == 64 and out["extreme_one"] == 1 and out["extreme_two"] == 2
