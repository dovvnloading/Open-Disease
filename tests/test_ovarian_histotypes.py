from fractions import Fraction
import itertools
import math
from pathlib import Path
import sys

import pytest
from scipy.stats.contingency import odds_ratio

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import ovarian_histotypes as histotypes


def test_nordic_partition_and_ci_orientation():
    rows = [r for r in histotypes.load("nordic_counts.csv") if r["window"] == "10_plus"]
    for row, observed, expected, ratio in zip(rows, [35, 20], [8.9578, 4.7891], [2.7814, 2.6971], strict=True):
        result = histotypes.nordic(row)
        assert result["pair_observed"] == observed
        assert result["expected_pair"] == pytest.approx(expected, abs=0.0001)
        assert result["ratio"] == pytest.approx(ratio, abs=0.0001)
        assert 1 < result["ci_rounding_low"] <= result["low"] < result["ratio"] < result["high"]
        assert result["high"] <= result["ci_rounding_high"]


def test_dutch_conditional_or_and_stratified_tail():
    rows = histotypes.load("dutch_counts.csv")
    single = histotypes.dutch(rows[:2])
    reference = odds_ratio([[13, 19], [27, 88]])
    assert single["conditional_odds_ratio"] == pytest.approx(reference.statistic)
    assert [single["low"], single["high"]] == pytest.approx(reference.confidence_interval())

    result = histotypes.dutch(rows)
    assert result["group_pair_and_total"] == {"endometriosis": [30, 69], "nevus": [41, 241]}
    # Independently enumerate hypergeometric masses using exact integer counts.
    weights = []
    for group_size, selected, total in [(32, 40, 147), (28, 23, 109), (9, 8, 54)]:
        denominator = math.comb(total, group_size)
        weights.append([
            Fraction(math.comb(selected, k) * math.comb(total - selected, group_size - k), denominator)
            for k in range(min(group_size, selected) + 1)
        ])
    tail = sum(
        math.prod(x[1] for x in combination)
        for combination in itertools.product(*(enumerate(w) for w in weights))
        if sum(x[0] for x in combination) >= 30
    )
    assert result["two_sided_p"] == pytest.approx(float(2 * tail), rel=1e-9)
