"""Analyze published ovarian cancer histotype counts. Run from any directory."""
import csv
import itertools
import json
from decimal import Decimal
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.stats import binomtest, hypergeom, nchypergeom_fisher


DATA = Path(__file__).resolve().parents[1] / "evidence" / "long-lag-histotypes"
HISTOTYPES = ("total", "clear_cell", "endometrioid")


def load(name):
    with (DATA / name).open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def nordic(row):
    n, cc, en = [int(row[k]) for k in HISTOTYPES]
    assert 0 < cc + en < n
    sir = [Decimal(row[k + "_sir"]) for k in HISTOTYPES]
    # Treat reported numbers as rounded at their displayed precision.
    half_units = [Decimal("0.5").scaleb(x.as_tuple().exponent) for x in sir]
    k = cc + en
    share_ci = binomtest(k, n).proportion_ci(method="exact")

    def contrast(sirs):
        expected_total, expected_cc, expected_en = [
            obs / float(s) for obs, s in zip((n, cc, en), sirs, strict=True)
        ]
        expected_pair = expected_cc + expected_en
        p0 = expected_pair / expected_total
        assert 0 < p0 < 1
        def scaled_odds(p):
            return (p / (1 - p)) / (p0 / (1 - p0))
        return {
            "expected_pair": expected_pair,
            "expected_other": expected_total - expected_pair,
            "pair_sir": k / expected_pair,
            "remaining_sir": (n - k) / (expected_total - expected_pair),
            "ratio": scaled_odds(k / n),
            "low": scaled_odds(share_ci.low),
            "high": scaled_odds(share_ci.high),
            "two_sided_p": binomtest(k, n, p0).pvalue,
            "excess_share": (k - expected_pair) / (n - expected_total),
        }

    primary = contrast(sir)
    corners = [contrast([v + s * u for v, s, u in zip(sir, signs, half_units, strict=True)])
               for signs in itertools.product((-1, 1), repeat=3)]
    return {"registry": row["registry"], "window": row["window"], "pair_observed": k,
            "remaining_observed": n - k, **primary,
            "ratio_rounding_min": min(x["ratio"] for x in corners),
            "ratio_rounding_max": max(x["ratio"] for x in corners),
            "ci_rounding_low": min(x["low"] for x in corners),
            "ci_rounding_high": max(x["high"] for x in corners)}


def dutch(rows, *, exposure_group="endometriosis", classify_unspecified_controls=False):
    """Exact conditional common case-mix odds ratio stratified by time window."""
    windows = list(dict.fromkeys(r["window"] for r in rows))
    strata, draws = [], 0
    cases_by_group = {exposure_group: [0, 0], "nevus": [0, 0]}
    for window in windows:
        pair = {r["exposure"]: r for r in rows if r["window"] == window}
        counts = []
        for exposure in (exposure_group, "nevus"):
            row = pair[exposure]
            a, total = int(row["clear_cell"]) + int(row["endometrioid"]), int(row["total"])
            if classify_unspecified_controls and exposure == "nevus":
                a += int(row["unspecified_adenocarcinoma"])
            assert 0 <= a <= total
            counts.append((a, total))
            cases_by_group[exposure][0] += a
            cases_by_group[exposure][1] += total
        (a, m), (b, n) = counts
        strata.append((m + n, a + b, m))
        draws += a

    def distribution(log_or, null=False):
        probs = np.array([1.0])
        for M, K, m in strata:
            x = np.arange(m + 1)
            if null:
                mass = hypergeom.pmf(x, M, K, m)
            else:
                mass = nchypergeom_fisher.pmf(x, M, K, m, np.exp(log_or))
            probs = np.convolve(probs, mass)
        return probs / probs.sum()

    def root(target):
        return float(np.exp(brentq(target, -12, 12)))

    def score(log_or):
        p = distribution(log_or)
        return p @ np.arange(len(p)) - draws

    null = distribution(0, null=True)
    return {
        "group_pair_and_total": cases_by_group,
        "conditional_odds_ratio": root(score),
        "low": root(lambda t: distribution(t)[draws:].sum() - 0.025),
        "high": root(lambda t: distribution(t)[:draws + 1].sum() - 0.025),
        "two_sided_p": min(1., 2 * min(null[draws:].sum(), null[:draws + 1].sum())),
    }


if __name__ == "__main__":
    dutch_counts = load("dutch_counts.csv")
    results = {
        "nordic": [nordic(row) for row in load("nordic_counts.csv")],
        "dutch": dutch(dutch_counts),
        "dutch_adenomyosis": dutch(dutch_counts, exposure_group="adenomyosis"),
        "dutch_all_unspecified_control_cancers_reassigned": dutch(
            dutch_counts, classify_unspecified_controls=True
        ),
    }
    print(json.dumps(results, indent=2, allow_nan=False))
