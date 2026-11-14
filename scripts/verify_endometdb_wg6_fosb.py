"""Independently verify the frozen Finnish contrast, with a discrete null distribution.

Shares release decoding and clinical source joins with the metadata census.
Uses sorted intensities for ranks and convolves exact per-stratum distributions.
Reads only the eligible donors' avg_signal field. Emits cohort totals only.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
import json
from pathlib import Path
import zlib

import numpy as np

from audit_endometdb import eligible, read_metadata
from audit_endometdb_chip_census import RAW_TAG, input_bytes, pg_table_blocks, raw_dump_bytes
from audit_endometdb_probe_keys import MANIFESTS, manifest

RUN_PHASE = {
    ("10", "secretory"): (2, 6), ("11", "proliferative"): (2, 1),
    ("12", "proliferative"): (1, 1), ("12", "secretory"): (1, 2),
    ("13", "proliferative"): (2, 1), ("13", "secretory"): (1, 1),
    ("14", "proliferative"): (3, 1),
}
ADDRESS = 6040086


def read_selected(directory, prefix):
    clinical = read_metadata(directory / "archive-start.bin", directory / "archive-tail.bin")
    _, selected = eligible(clinical)
    cohort = {r["analysis_result_id"].encode(): (r["analysis_run_id"], r["phase"], r["group"])
              for r in selected if (r["analysis_run_id"], r["phase"]) in RUN_PHASE}
    expected = Counter({(run, phase, group): n for (run, phase), sizes in RUN_PHASE.items()
                        for group, n in zip(("case", "control"), sizes, strict=True)})
    if len(selected) != 38 or Counter(cohort.values()) != expected:
        raise ValueError("Clinical record counts differ from the frozen census")
    manifest_name = "GPL13376-platform.gz"
    mapping = manifest(directory / manifest_name, MANIFESTS[manifest_name])
    if len(mapping) != 48701 or mapping[ADDRESS] != "ILMN_1730706":
        raise ValueError("Unexpected WG-6 manifest")
    per_assay = defaultdict(dict)
    inflater = zlib.decompressobj()
    remainder = b""
    complete = False
    chunks = pg_table_blocks(raw_dump_bytes(input_bytes(directory, prefix)), RAW_TAG)
    for chunk in chunks:
        lines = (remainder + inflater.decompress(chunk)).split(b"\n")
        remainder = lines.pop()
        for line in lines:
            if line == b"\\.":
                if complete:
                    raise ValueError("Repeated COPY terminator")
                complete = True
                continue
            if not line and complete:
                continue
            if complete:
                raise ValueError("Unexpected record after COPY terminator")
            _, assay, suffix = line.split(b"\t", 2)
            if assay in cohort:
                address, signal, _ = suffix.split(b"\t", 2)
                k = int(address)
                if k in per_assay[assay]:
                    raise ValueError("Repeated address in eligible WG-6 assay")
                try:
                    per_assay[assay][k] = float(signal)
                except ValueError:
                    raise ValueError("Missing or malformed selected signal") from None
    if not complete or remainder or not inflater.eof or inflater.unused_data or set(per_assay) != set(cohort):
        raise ValueError("Incomplete input")
    groups = defaultdict(lambda: defaultdict(list))
    for assay, signals in per_assay.items():
        if signals.keys() != mapping.keys():
            raise ValueError("Eligible assay lacks exact WG-6 address set")
        ordered = np.sort(np.fromiter(signals.values(), dtype=float))
        if not np.isfinite(ordered).all():
            raise ValueError("Nonfinite selected signal")
        left = int(np.searchsorted(ordered, signals[ADDRESS], side="left"))
        right = int(np.searchsorted(ordered, signals[ADDRESS], side="right"))
        run, phase, group = cohort[assay]
        groups[run, phase][group].append(left + right)
    return groups


def mean_difference(case_sum, total_sum, ncase, ncontrol):
    return Fraction(case_sum, ncase) - Fraction(total_sum - case_sum, ncontrol)


def calculate(groups):
    distribution = Counter({Fraction(0): 1})
    effects = []
    for counts in groups.values():
        cases = counts["case"]
        combined = cases + counts["control"]
        sizes = len(cases), len(counts["control"])
        total = sum(combined)
        effects.append(mean_difference(sum(cases), total, *sizes))
        permutations = Counter(mean_difference(sum(x), total, *sizes)
                               for x in combinations(combined, len(cases)))
        updated = Counter()
        for left, n in distribution.items():
            for right, m in permutations.items():
                updated[left + right] += n * m
        distribution = updated
    observed = sum(effects)
    extreme = sum(n for val, n in distribution.items() if abs(val) >= abs(observed))
    scale = Fraction(100, 2 * 48701)
    drop = [float((observed - x) * scale / (len(effects) - 1)) for x in effects]
    return {
        "strata": len(effects), "effect_percentile_points": float(observed * scale / len(effects)),
        "extreme": extreme, "assignments": distribution.total(),
        "two_sided_p": extreme / distribution.total(),
        "delete_one_stratum_range": [min(drop), max(drop)],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--prefix", type=Path)
    args = parser.parse_args()
    print(json.dumps(calculate(read_selected(args.directory, args.prefix)), indent=2))
