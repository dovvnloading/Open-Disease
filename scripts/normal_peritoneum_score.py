"""Execute the frozen noncase tissue comparison after the private physical-chip gate."""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
from statistics import mean

import numpy as np
from scipy.stats import rankdata

from analyze_endometdb_wg6_fosb import N_FEATURES
from audit_endometdb_chip_census import copy_lines
from audit_peritoneal_geo import PANEL
from genome_endometdb import load_wg6_annotation, selected_signal_matrix

ROOT = Path(__file__).resolve().parents[1]
PANEL_PATH = ROOT / "evidence/lesion-transfer/panel.json"
# Entrez -> (signed direction, numeric WG-6 physical address), frozen before access.
TARGETS = {55026: (1, 60064), 9077: (1, 4570725), 3512: (1, 6450672),
           2568: (-1, 3840113), 6299: (-1, 2360678)}


def frozen_columns(annotation, panel_path):
    raw = panel_path.read_bytes()
    if sha256(raw).hexdigest() != PANEL:
        raise ValueError("Frozen panel hash changed")
    rows = json.loads(raw)
    directions = {row["entrez"]: row["direction"] for row in rows}
    if len(rows) != len(TARGETS) or directions != {g: v[0] for g, v in TARGETS.items()}:
        raise ValueError("Frozen gene IDs or directions changed")
    addresses = sorted(annotation)
    observed = {a: annotation[a] for a in addresses if annotation[a] in TARGETS}
    if len(addresses) != N_FEATURES or observed != {a: g for g, (_, a) in TARGETS.items()}:
        raise ValueError("WG-6 denominator or frozen physical target addresses changed")
    return addresses, [addresses.index(address) for _, address in TARGETS.values()]


def assay_score(values, target_columns):
    """Score one assay against every physical feature; also check each midrank."""
    signs = [sign for sign, _ in TARGETS.values()]
    ranked = rankdata(values, method="average")
    signed_sum = 0
    for col, sign in zip(target_columns, signs, strict=True):
        doubled = int(2 * ranked[col] - 1)
        signal = values[col]
        counted = 2 * np.count_nonzero(values < signal) + np.count_nonzero(values == signal)
        if doubled != counted:
            raise ValueError("Independent full-reference midranks disagree")
        signed_sum += sign * doubled
    return Fraction(50 * signed_sum, len(signs) * len(values))


def tail_count(deltas, shift=Fraction(0)):
    centered = tuple(value - shift for value in deltas)
    observed = abs(sum(centered))
    return sum(abs(sum(s * d for s, d in zip(signs, centered))) >= observed
               for signs in product((-1, 1), repeat=len(deltas)))


def common_shift_interval(deltas):
    """Invert exactly at rational breakpoints, counting equality in both tails."""
    n = len(deltas)
    assignments = 2 ** n
    # The two constant assignments enter at all shifts. Every other assignment
    # enters on the closed interval between its two group means. These intervals
    # all contain the overall mean, so their coverage cannot have internal gaps.
    needed = assignments // 20 + 1 - 2
    if needed <= 0:
        return None, None
    total = sum(deltas)
    lowers, uppers = [], []
    for k in range(1, n):
        for subset in combinations(deltas, k):
            subtotal = sum(subset)
            a, b = Fraction(subtotal, k), Fraction(total - subtotal, n - k)
            lowers.append(min(a, b))
            uppers.append(max(a, b))
    lo, hi = sorted(lowers)[needed - 1], sorted(uppers, reverse=True)[needed - 1]
    if not all(20 * tail_count(deltas, x) > assignments for x in (lo, hi)):
        raise ValueError("An endpoint fails exact test inversion")
    return lo, hi


def synthetic_controls():
    cols = list(range(len(TARGETS)))
    p = np.array([6, 4, 3, 1, 2, 1, 3], dtype=float)
    e = np.array([0, 2, 2, 8, 9, 3, 4], dtype=float)
    flat = np.ones_like(p)
    paired = [(flat, flat), (p, p), (p, e), (e, p), (p, e), (p, e)]
    forward = [assay_score(a, cols) - assay_score(b, cols) for a, b in paired]
    reverse = [assay_score(b, cols) - assay_score(a, cols) for a, b in paired]
    if (forward[:2] != [0, 0] or not forward[2] or reverse != [-d for d in forward]
            or tail_count(forward) != tail_count(reverse)):
        raise ValueError("Synthetic paired score or sign-flip control failed")


def paired_results(pairs, scores):
    """Aggregate technical assays within specimen before averaging across donors."""
    deltas, restricted, repeat_differences = [], [], []
    repeat_specimens = 0
    for pair in pairs:
        p = [scores[a] for a in pair.peritoneum_assays]
        e = [scores[a] for a in pair.uterus_assays]
        delta = mean(p) - mean(e)
        deltas.append(delta)
        if pair.sensitivity_phase:
            restricted.append(delta)
        for specimen in (p, e):
            if len(specimen) > 1:
                repeat_specimens += 1
                repeat_differences.extend(abs(a - b) for a, b in combinations(specimen, 2))
    n = len(deltas)
    observed = mean(deltas)
    extreme, total = tail_count(deltas), 2 ** n
    lo, hi = common_shift_interval(deltas)
    if (20 * extreme > total) != (lo <= 0 <= hi):
        raise ValueError("Exact compatibility set and zero-shift test disagree")
    deleted = [(sum(deltas) - d) / (n - 1) for d in deltas]
    decision = "unresolved"
    if 20 * extreme <= total:
        decision = "higher_in_peritoneum" if observed > 0 else "lower_in_peritoneum"
    technical = {"repeated_specimens": repeat_specimens, "assay_pairs": len(repeat_differences)}
    if repeat_specimens >= 2:
        technical["mean_absolute_assay_pair_difference_pp"] = float(mean(repeat_differences))
    return {
        "donors": n, "assays": len(scores), "difference_definition": "peritoneum minus uterus",
        "mean_difference_pp": float(observed),
        "positive_zero_negative": [sum(d > 0 for d in deltas), sum(d == 0 for d in deltas),
                                   sum(d < 0 for d in deltas)],
        "extreme_two_sided": extreme, "sign_assignments": total, "two_sided_p": extreme / total,
        "compatibility_95pct_common_shift_pp": [float(lo), float(hi)],
        "leave_one_donor_out_mean_range_pp": [float(min(deleted)), float(max(deleted))],
        "phase_sensitivity": {
            "donors": len(restricted),
            "phase_counts": dict(Counter(p.sensitivity_phase for p in pairs if p.sensitivity_phase)),
            "mean_difference_pp": float(mean(restricted)) if len(restricted) >= 2 else None,
        },
        "technical_control": technical, "decision": decision,
    }


def score_qualified_pairs(pairs, lines, annotation, panel_path=PANEL_PATH):
    if not 6 <= len(pairs) <= 8 or len({p.donor for p in pairs}) != len(pairs):
        raise ValueError("Outcome gate requires six to eight distinct qualified donor pairs")
    assays = sorted(a for p in pairs for a in (*p.peritoneum_assays, *p.uterus_assays))
    if any(not p.peritoneum_assays or not p.uterus_assays for p in pairs) or len(assays) != len(set(assays)):
        raise ValueError("Frozen pairs must contain unique assays in both tissues")
    addresses, columns = frozen_columns(annotation, panel_path)
    synthetic_controls()
    matrix = selected_signal_matrix(lines, [a.encode() for a in assays], addresses)
    scores = {a: assay_score(values, columns) for a, values in zip(assays, matrix, strict=True)}
    report = paired_results(pairs, scores)
    report.update(panel_sha256=PANEL, full_feature_denominator=len(addresses),
                  synthetic_controls_passed=True, verified_all_target_midranks=True)
    return report


def main():
    from normal_peritoneum_preflight import blocks, preflight

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Ignored directory with the pinned audit inputs")
    args = parser.parse_args()
    cohort, pairs = preflight(args.directory)
    result = {"preflight": cohort}
    if len(pairs) < 6:
        result["status"] = "under_resolved_before_expression"
    else:
        annotation = load_wg6_annotation(args.directory)
        result["outcome"] = score_qualified_pairs(pairs, copy_lines(blocks(args.directory)), annotation)
        result["status"] = "completed"
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
