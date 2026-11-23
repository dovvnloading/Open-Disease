"""Reanalyse the published scan-disjoint matrisome rankings using aggregate counts."""

import hashlib
import json
from math import lcm
from pathlib import Path

from matrisome_omission import omission_summary
from matrisome_rank_inference import (
    exact_rank_test, exact_test_from_distribution, holm_three_targets,
    ranked_sum_distribution,
)
from matrisome_scan_bounds import feasible_within_wins

EVIDENCE = Path(__file__).resolve().parents[1] / "evidence"
BENCHMARK = "matrisome-holdout-2026-11-15/transfer-benchmark.json"
SCAN_DESIGN = "genome-transport-2026-11-15/rome_scan_transport.json"
DIGESTS = {
    BENCHMARK: "eca5133ab7fa796492d40d306026bbd20912626c599fae77072e9e21637baaa5",
    SCAN_DESIGN: "46824dbfb93e62ec65e3a18aa5923f1493be174623169d368cd13b6a52355e9a",
}


def read_pinned(relative_path):
    data = (EVIDENCE / relative_path).read_bytes()
    if hashlib.sha256(data).hexdigest() != DIGESTS[relative_path]:
        raise ValueError(f"Aggregate input differs from the audited release: {relative_path}")
    return json.loads(data)


def stratum(row):
    cases, comparators, higher = row["cases"], row["comparators"], row["higher"]
    if row["ties"] != 0:
        raise ValueError("Cross-group ties would require additional rank information")
    if cases * comparators != row["pairs"] or not 0 <= higher <= row["pairs"]:
        raise ValueError("Inconsistent reported pair counts")
    return cases, comparators, higher


def us_phase_mix(phases):
    rows = list(phases.values())
    sizes = [cases * comparators for cases, comparators, _ in rows]
    denominator = lcm(*sizes)
    total_wins = sum(row[2] for row in rows)
    total_pairs = sum(sizes)
    # In the equal-phase average each win is worth 1 / (number of phases * phase pairs).
    weights = [denominator // pairs for pairs in sizes]
    equal = exact_rank_test(rows, weights)
    leave_out = {}
    for name, (_, _, wins), pairs in zip(phases, rows, sizes, strict=True):
        leave_out[name] = {"higher": total_wins - wins, "pairs": total_pairs - pairs}
    return {
        "by_phase": {name: {"cases": cases, "comparators": comparators, "higher": wins,
                            "pairs": cases * comparators}
                     for name, (cases, comparators, wins) in phases.items()},
        "equal_phase_fraction": equal["observed"] / equal["maximum"],
        "equal_phase_test": equal,
        "omit_one_phase_pair_weighted": leave_out,
    }


def rome_scan_sensitivity(rome, design):
    cases, comparators, higher = rome
    mixed = [(day["cases"], day["controls"]) for day in design["date_specific"]]
    mixed_cases = sum(n for n, _ in mixed)
    mixed_controls = sum(m for _, m in mixed)
    if (mixed_cases != design["cases"] or mixed_controls != design["controls"]
            or comparators != mixed_controls or cases - mixed_cases != design["excluded_arrays"]):
        raise ValueError("Rome scanner margins do not reconcile with the ranking cohort")
    # The remaining case-only arrays have no within-day case–comparator pairs.
    blocks = mixed + [(cases - mixed_cases, 0)]
    feasible = feasible_within_wins(blocks, higher)
    null = ranked_sum_distribution(mixed)
    sensitivity = []
    for wins in feasible:
        tail = exact_test_from_distribution(null, wins)
        sensitivity.append({"possible_within_day_wins": wins, "extreme": tail["extreme"],
                            "two_sided_p": tail["two_sided_p"]})
    return {
        "scan_day_margins": [{"cases": n, "comparators": m} for n, m in blocks],
        "within_day_pairs": sum(n * m for n, m in mixed),
        "within_day_assignments": sum(null.values()),
        "cases_scanned_after_comparator_within_day": design["same_date_case_scanned_after_control"],
        "compatible_within_day_tests": sensitivity,
    }


def analyse(benchmark, design):
    fit = benchmark["results"]["scan_disjoint"]
    phases = {name: stratum(row) for name, row in fit["GSE6364_by_phase"].items()}
    rows = list(phases.values())
    pooled = fit["GSE6364_pooled_within_phase"]
    if (sum(n * m for n, m, _ in rows), sum(u for _, _, u in rows), 0) != (
            pooled["pairs"], pooled["higher"], pooled["ties"]):
        raise ValueError("Pooled GSE6364 total does not match the released phases")

    targets = {"GSE6364": rows, "Rome": [stratum(fit["Rome"])], "Mumbai": [stratum(fit["Mumbai"])]}
    tests = {name: exact_rank_test(strata) for name, strata in targets.items()}
    adjusted = holm_three_targets({name: test["two_sided_p"] for name, test in tests.items()})
    for name in tests:
        tests[name]["holm_p"] = adjusted[name]

    return {
        "input_sha256": DIGESTS,
        "fit": "scan_disjoint",
        "array_label_tests": tests,
        "us_phase_sensitivity": us_phase_mix(phases),
        "rome_scan_day_sensitivity": rome_scan_sensitivity(targets["Rome"][0], design),
        "fixed_score_omission_bounds": {
            name: omission_summary(strata) for name, strata in targets.items()
        },
    }


if __name__ == "__main__":
    print(json.dumps(analyse(read_pinned(BENCHMARK), read_pinned(SCAN_DESIGN)),
                     indent=2, sort_keys=True))
