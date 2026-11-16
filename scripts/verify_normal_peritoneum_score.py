"""Independent arithmetic for the paired noncase experiment, emitting only aggregates."""
from __future__ import annotations

import argparse
from bisect import bisect_left, bisect_right
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
import math
from pathlib import Path
from statistics import mean

from audit_endometdb_chip_census import copy_lines
from audit_endometdb_probe_keys import MANIFESTS, manifest
from normal_peritoneum_preflight import blocks, preflight

SIGNED_ADDRESSES = {60064: 1, 4570725: 1, 6450672: 1, 3840113: -1, 2360678: -1}


def verified_assay_scores(lines, assay_ids, addresses):
    """Compute signed midrank percentiles independently via sorted-value insertion points."""
    assays = {a.encode(): {} for a in assay_ids}
    for line in lines:
        _, assay, rest = line.split(b"\t", 2)
        if assay not in assays:
            continue
        addr_text, signal_text, _ = rest.split(b"\t", 2)
        addr, signal = int(addr_text), float(signal_text)
        values = assays[assay]
        if addr not in addresses or addr in values or not math.isfinite(signal):
            raise ValueError("Invalid address, duplicate or nonfinite signal in verified array")
        values[addr] = signal
    scores = {}
    for assay, values in assays.items():
        if values.keys() != addresses:
            raise ValueError("Incomplete physical address reference in verified assay")
        ordered = sorted(values.values())
        numerator = sum(sign * (bisect_left(ordered, values[addr]) + bisect_right(ordered, values[addr]))
                        for addr, sign in SIGNED_ADDRESSES.items())
        scores[assay.decode()] = Fraction(100 * numerator, 2 * len(addresses) * len(SIGNED_ADDRESSES))
    return scores


def verified_scores(lines, pairs, addresses):
    """Calculate the original noncase tissue comparison from independently ranked assays."""
    if len(pairs) < 6:
        raise ValueError("Too few physically verified donor pairs")
    assay_ids = [a for pair in pairs for a in pair.peritoneum_assays + pair.uterus_assays]
    scores = verified_assay_scores(lines, assay_ids, addresses)
    diffs = []
    sensitivity = []
    rep_gaps = []
    for pair in pairs:
        p, u = ([scores[a] for a in ids] for ids in (pair.peritoneum_assays, pair.uterus_assays))
        delta = mean(p) - mean(u)
        diffs.append(delta)
        if pair.sensitivity_phase:
            sensitivity.append(delta)
        for tissue in (p, u):
            rep_gaps.extend(abs(a - b) for a, b in combinations(tissue, 2))
    return diffs, sensitivity, rep_gaps


def verify_inference(differences):
    """Directly evaluate each distinct exact breakpoint in rational score units."""
    n = len(differences)
    if not 6 <= n <= 8:
        raise ValueError("Outside the frozen verified donor range")
    unit = math.lcm(*(x.denominator for x in differences))
    values = [int(d * unit) for d in differences]
    total = sum(values)
    subsets = [(sum(part), k) for k in range(n + 1) for part in combinations(values, k)]
    assignments = 2**n

    def tail(shift):
        a, b = shift.as_integer_ratio()
        observed = abs(b * total - a * n)
        return sum(abs(b * (total - 2 * s) - a * (n - 2 * k)) >= observed for s, k in subsets)

    knots = sorted({Fraction(s, k) for s, k in subsets if k})
    accepted = [v for v in knots if 20 * tail(v) > assignments]
    if not accepted:
        raise ValueError("No accepted test breakpoint")
    lo, hi = accepted[0], accepted[-1]
    probes = [knots[0] - 1, knots[-1] + 1, *knots, *((a + b) / 2 for a, b in zip(knots, knots[1:]))]
    if any((20 * tail(v) > assignments) != (lo <= v <= hi) for v in probes):
        raise ValueError("Compatibility set failed direct tail verification")
    effect = mean(differences)
    deletions = [(sum(differences) - d) / (n - 1) for d in differences]
    counts = Counter((d > 0) - (d < 0) for d in differences)
    return {
        "donors": n, "mean_difference_pp": float(effect),
        "positive_zero_negative": [counts[1], counts[0], counts[-1]],
        "extreme_two_sided": tail(Fraction(0)), "sign_assignments": assignments,
        "two_sided_p": tail(Fraction(0)) / assignments,
        "compatibility_95pct_common_shift_pp": [float(lo / unit), float(hi / unit)],
        "leave_one_donor_out_mean_range_pp": [float(min(deletions)), float(max(deletions))],
    }


def verify(directory, recorded):
    _, pairs = preflight(directory)
    if len(pairs) < 6:
        raise ValueError("Preflight stopped the outcome experiment")
    filename = "GPL13376-platform.gz"
    addresses = set(manifest(directory / filename, MANIFESTS[filename]))
    if len(addresses) != 48701 or not addresses.issuperset(SIGNED_ADDRESSES):
        raise ValueError("Incorrect WG-6 reference")
    differences, restricted, repeats = verified_scores(copy_lines(blocks(directory)), pairs, addresses)
    checks = verify_inference(differences)
    if recorded["status"] != "completed":
        raise ValueError("Missing primary execution")
    primary = recorded["outcome"]
    if any(primary[key] != value for key, value in checks.items()):
        raise ValueError("Independent primary, interval or deletion calculation differs")
    sensitivity = primary["phase_sensitivity"]
    expected = float(mean(restricted)) if len(restricted) >= 2 else None
    if sensitivity["donors"] != len(restricted) or sensitivity["mean_difference_pp"] != expected:
        raise ValueError("Independent sensitivity calculation differs")
    technical = primary["technical_control"]
    if technical["assay_pairs"] != len(repeats):
        raise ValueError("Independent technical repeat count differs")
    return {"verified": True, "donors": len(pairs), "primary_fields_checked": sorted(checks),
            "checked_phase_sensitivity": True, "checked_technical_pair_count": True}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.directory, json.loads(args.results.read_text())), indent=2))
