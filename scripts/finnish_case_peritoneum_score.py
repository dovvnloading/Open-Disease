"""Lock or, in a separately authorized task, execute the frozen Finnish case-pair test."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
from statistics import mean

from audit_endometdb_chip_census import copy_lines
from audit_finnish_case_peritoneum_design import audit_with_pairs
from audit_peritoneal_geo import PANEL
from genome_endometdb import load_wg6_annotation, selected_signal_matrix
from normal_peritoneum_preflight import blocks, pair_lock
from normal_peritoneum_score import (
    PANEL_PATH, assay_score, common_shift_interval, frozen_columns, synthetic_controls, tail_count,
)
from verify_normal_peritoneum_score import verified_assay_scores, verify_inference

DESIGN = Path(__file__).resolve().parents[1] / "evidence/peritoneal-followup/finnish-design.json"
RUN_COUNTS = {"9": 1, "10": 5, "12": 1}
LOCK_NAME = "finnish-case-peritoneum-verified-pairs.json"


def locked_assays(pairs):
    if len(pairs) != 7 or len({p.donor for p in pairs}) != 7 or Counter(p.run for p in pairs) != RUN_COUNTS:
        raise ValueError("Expected seven distinct cases in the frozen three-run design")
    ids = [a for p in pairs for a in (*p.lesion, *p.peritoneum)]
    specimens = [s for p in pairs for s in (p.lesion_specimen, p.peritoneum_specimen)]
    if (len(ids) != 14 or len(set(ids)) != 14 or len(set(specimens)) != 14
            or any(not p.lesion or not p.peritoneum for p in pairs)):
        raise ValueError("Expected 14 distinct specimens and WG-6 assays across the seven pairs")
    return sorted(ids)


def preflight(directory, create):
    report, pairs = audit_with_pairs(directory)
    if report != json.loads(DESIGN.read_text()):
        raise ValueError("Address-only audit differs from the frozen case and control design")
    locked_assays(pairs)
    pair_lock(directory / LOCK_NAME, pairs, create=create)
    return pairs


def compare(pairs, lines, annotation):
    """Read the streamed rows twice and release aggregates only after independent agreement."""
    ids = locked_assays(pairs)
    addresses, columns = frozen_columns(annotation, PANEL_PATH)
    synthetic_controls()
    matrix = selected_signal_matrix(lines(), [a.encode() for a in ids], addresses)
    scores = {a: assay_score(row, columns) for a, row in zip(ids, matrix, strict=True)}
    deltas = [mean(scores[a] for a in p.lesion) - mean(scores[a] for a in p.peritoneum) for p in pairs]
    n = len(deltas)
    extreme, assignments = tail_count(deltas), 2**n
    effect = mean(deltas)
    lo, hi = common_shift_interval(deltas)
    deletions = [(sum(deltas) - d) / (n - 1) for d in deltas]
    report = {
        "donors": n, "mean_difference_pp": float(effect),
        "positive_zero_negative": [sum(d > 0 for d in deltas), sum(d == 0 for d in deltas), sum(d < 0 for d in deltas)],
        "extreme_two_sided": extreme, "sign_assignments": assignments, "two_sided_p": extreme / assignments,
        "compatibility_95pct_common_shift_pp": [float(lo), float(hi)],
        "leave_one_donor_out_mean_range_pp": [float(min(deletions)), float(max(deletions))],
    }

    independent = verified_assay_scores(lines(), ids, set(addresses))
    if independent != scores:
        raise ValueError("Independent per-array score disagrees with frozen signed ranks")
    alternate = []
    reduced = []
    for p in pairs:
        delta = sum(independent[a] for a in p.lesion) / len(p.lesion)
        delta -= sum(independent[a] for a in p.peritoneum) / len(p.peritoneum)
        alternate.append(delta)
        if p.run != "10":
            reduced.append(delta)
    if report != verify_inference(alternate):
        raise ValueError("Independent exact tail, interval or aggregate disagrees")
    remaining = [d for p, d in zip(pairs, deltas, strict=True) if p.run != "10"]
    if len(remaining) != 2 or mean(remaining) != sum(reduced, Fraction(0)) / 2:
        raise ValueError("Independent leave-run-10-out estimate disagrees")
    decision = "unresolved"
    if extreme * 20 <= assignments:
        decision = "higher_in_labelled_lesion" if effect > 0 else "lower_in_labelled_lesion"
    report.update(
        assays=len(ids), difference_definition="surgically labelled superficial lesion minus case peritoneum",
        panel_sha256=PANEL, full_feature_denominator=len(addresses),
        leave_run_10_out={"donors": 2, "mean_difference_pp": float(mean(remaining))},
        synthetic_controls_passed=True, independent_verification_passed=True, decision=decision,
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("lock", "execute"),
                        help="lock reads metadata and addresses only; execute reads targeted case signals")
    parser.add_argument("directory", type=Path, help="Ignored cache used for the address-only audit")
    args = parser.parse_args()
    pairs = preflight(args.directory, create=args.mode == "lock")
    if args.mode == "lock":
        report = {"status": "design_locked", "donors": len(pairs), "assays": len(locked_assays(pairs))}
    else:
        report = compare(pairs, lambda: copy_lines(blocks(args.directory)), load_wg6_annotation(args.directory))
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
