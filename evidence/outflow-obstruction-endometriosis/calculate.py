"""Compute cohort-specific prevalence differences from published group counts."""

import json
from fractions import Fraction
from math import ceil, hypot, sqrt
from pathlib import Path
from statistics import NormalDist


def wilson_interval(events, total):
    if not 0 <= events <= total or total <= 0:
        raise ValueError(f"Invalid binomial counts: {events}/{total}")
    z = NormalDist().inv_cdf(0.975)
    p = events / total
    center = events + z**2 / 2
    half_width = z * sqrt(total * p * (1 - p) + z**2 / 4)
    return (center - half_width) / (total + z**2), (center + half_width) / (total + z**2)


def contrast(events, total, ref_events, ref_total):
    low, high = wilson_interval(events, total)
    ref_low, ref_high = wilson_interval(ref_events, ref_total)
    p, ref_p = events / total, ref_events / ref_total
    difference = p - ref_p
    return {
        "difference_pp": round(100 * difference, 6),
        "newcombe_95_pp": [
            round(100 * (difference - hypot(p - low, ref_high - ref_p)), 6),
            round(100 * (difference + hypot(high - p, ref_p - ref_low)), 6),
        ],
    }


def tipping_point(events, total, ref_events, ref_total):
    p = Fraction(events, total)
    added = max(0, ceil(p * ref_total) - ref_events)
    noncases = ref_total - ref_events
    at_threshold = p - Fraction(ref_events + added, ref_total)
    before = p - Fraction(ref_events + added - 1, ref_total) if added else None
    if added > noncases or at_threshold > 0 or (before is not None and before <= 0):
        raise ArithmeticError("Threshold failed the adjacent-count check")
    return {
        "additional_reference_cases": added,
        "recorded_reference_noncases": noncases,
        "percent_of_reference_noncases": round(100 * added / noncases, 6) if noncases else None,
        "difference_one_case_before_pp": round(100 * float(before), 6) if added else None,
        "difference_at_threshold_pp": round(100 * float(at_threshold), 6),
    }


def calculate(data):
    rows = []
    for cohort in data["cohorts"]:
        total = cohort["higher_retention_group"]["women"]
        ref_total = cohort["reference_group"]["women"]
        for outcome in cohort["outcomes"]:
            events = outcome["higher_retention_events"]
            ref_events = outcome["reference_events"]
            row = {"study_id": cohort["study_id"], "outcome": outcome["outcome"]}
            row.update(contrast(events, total, ref_events, ref_total))
            if outcome["outcome"] == "endometriosis":
                row["tipping_point"] = tipping_point(events, total, ref_events, ref_total)
            rows.append(row)
    return rows


if __name__ == "__main__":
    counts = json.loads(Path(__file__).with_name("counts.json").read_text(encoding="utf-8"))
    print(json.dumps(calculate(counts), indent=2))
