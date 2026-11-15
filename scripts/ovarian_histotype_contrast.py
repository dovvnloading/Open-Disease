"""Recompute long-lag histotype contrasts from published registry aggregates."""

import csv
import json
from decimal import Decimal
from itertools import product
from pathlib import Path

from scipy.stats import binomtest, fisher_exact
from scipy.stats.contingency import odds_ratio


INPUT_DIR = Path(__file__).resolve().parents[1] / "evidence" / "ovarian-histotypes-2026-11-15"


def read_table(filename):
    with (INPUT_DIR / filename).open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def odds(proportion):
    return proportion / (1 - proportion)


def target_fraction(counts, ratios):
    """Expected fraction of cancers that are clear cell or endometrioid."""
    total, clear_cell, endometrioid = (n / sir for n, sir in zip(counts, ratios, strict=True))
    return (clear_cell + endometrioid) / total


def rounding_bounds(printed):
    """Use precision as printed, including trailing zeroes."""
    value = Decimal(printed)
    half_step = Decimal(5).scaleb(value.as_tuple().exponent - 1)
    return (float(value - half_step), float(value + half_step))


def standardized_contrast(row):
    categories = ("total", "clear_cell", "endometrioid")
    counts = [int(row[name]) for name in categories]
    ratio_text = [row[name + "_sir"] for name in categories]
    ratios = [float(x) for x in ratio_text]
    observed = sum(counts[1:])
    total = counts[0]
    reference = target_fraction(counts, ratios)
    result = binomtest(observed, total, reference)
    ci = result.proportion_ci(confidence_level=0.95, method="exact")
    expected_total = total / ratios[0]
    expected_target = reference * expected_total
    rounded_references = [
        target_fraction(counts, corner) for corner in product(*(rounding_bounds(x) for x in ratio_text))
    ]
    # Ratio of SIRs decreases monotonically with the expected target fraction.
    low_ref, high_ref = min(rounded_references), max(rounded_references)
    return {
        "country": row["country"],
        "observed_target": observed,
        "observed_total": total,
        "expected_target": expected_target,
        "expected_total": expected_total,
        "target_sir": observed / expected_target,
        "complement_sir": (total - observed) / (expected_total - expected_target),
        "ratio_of_sirs": odds(observed / total) / odds(reference),
        "ci95": [odds(ci.low) / odds(reference), odds(ci.high) / odds(reference)],
        "p_two_sided": result.pvalue,
        "ratio_rounding_range": [odds(observed / total) / odds(x) for x in (high_ref, low_ref)],
        "ci95_rounding_envelope": [odds(ci.low) / odds(high_ref), odds(ci.high) / odds(low_ref)],
    }


def dutch_contrast(rows, exposure):
    strata = {}
    for row in rows:
        groups = strata.setdefault(row["years"], {})
        group = row["group"]
        if group in groups:
            raise ValueError(f"Duplicate {group} in {row['years']}")
        total = int(row["total"])
        target = int(row["clear_cell"]) + int(row["endometrioid"])
        if not 0 <= target <= total:
            raise ValueError(f"Counts inconsistent in {row}")
        groups[group] = (target, total - target)
    tables = [(groups[exposure], groups["nevus"]) for groups in strata.values()]
    aggregate = [[sum(table[i][j] for table in tables) for j in range(2)] for i in range(2)]
    estimated = odds_ratio(aggregate, kind="conditional")
    ci = estimated.confidence_interval()
    numerator = denominator = 0.0
    for (a, b), (c, d) in tables:
        n = a + b + c + d
        numerator += a * d / n
        denominator += b * c / n
    return {
        "exposure": exposure,
        "rows_exposure_nevus_cols_target_complement": aggregate,
        "conditional_case_mix_odds_ratio": estimated.statistic,
        "ci95": [ci.low, ci.high],
        "p_two_sided": fisher_exact(aggregate).pvalue,
        "follow_up_stratified_mantel_haenszel_odds_ratio": numerator / denominator,
    }


def analyze():
    return {
        "standardized": [standardized_contrast(row) for row in read_table("standardized.tsv")],
        "dutch": [
            dutch_contrast(read_table("dutch_counts.tsv"), exposure)
            for exposure in ("endometriosis", "adenomyosis")
        ],
    }


if __name__ == "__main__":
    print(json.dumps(analyze(), indent=2, allow_nan=False))
