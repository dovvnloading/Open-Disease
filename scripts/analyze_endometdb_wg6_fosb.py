"""Run the frozen Finnish WG-6 percentile analysis from the pinned offline archive."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
import json
import math
from pathlib import Path

from audit_endometdb import HEAD_SHA256, eligible, pinned, read_metadata
from audit_endometdb_chip_census import (
    RAW_START, RAW_TAG, census, copy_lines, input_bytes, pg_table_blocks, raw_dump_bytes, raw_keys,
)
from audit_endometdb_probe_keys import MANIFESTS, manifest

TARGET = 6040086
N_FEATURES = 48_701
EXPECTED_CELLS = {
    ("10", "secretory"): (2, 6),
    ("11", "proliferative"): (2, 1),
    ("12", "proliferative"): (1, 1),
    ("12", "secretory"): (1, 2),
    ("13", "proliferative"): (2, 1),
    ("13", "secretory"): (1, 1),
    ("14", "proliferative"): (3, 1),
}
RAW_COPY = (b"COPY public.microarray_data (id, sample_analysis_result_id, probeid, avg_signal, "
            b"detection_pval, avg_nbeads, bead_stderr) FROM stdin;")
RAW_FK = (b"ADD CONSTRAINT sample_analysis_result_microarray_data_fk FOREIGN KEY "
          b"(sample_analysis_result_id) REFERENCES public.sample_analysis_result(id)")


def shared_design(selected, status):
    """Recover overlap independently of signal and enforce the frozen census."""
    if (len(selected) != 38 or len({row["donor_id"] for row in selected}) != 38
            or len({row["analysis_result_id"] for row in selected}) != 38
            or Counter(status.values()) != {"WG6": 32, "repeated_addresses": 6}):
        raise ValueError("Clinical donor or physical chip census differs from the plan")
    wg6 = [row for row in selected if status[row["analysis_result_id"].encode()] == "WG6"]
    if Counter(r["group"] for r in wg6) != {"case": 15, "control": 17}:
        raise ValueError("Unexpected clinically eligible WG-6 group counts")
    counts = defaultdict(Counter)
    for row in wg6:
        counts[(row["analysis_run_id"], row["phase"])][row["group"]] += 1
    shared = {key: (c["case"], c["control"]) for key, c in counts.items() if c["case"] and c["control"]}
    if shared != EXPECTED_CELLS:
        raise ValueError("Shared run-by-phase donor cells differ from the plan")
    return {r["analysis_result_id"].encode(): ((r["analysis_run_id"], r["phase"]), r["group"])
            for r in wg6 if (r["analysis_run_id"], r["phase"]) in shared}


def rank_signals(lines, design, addresses, target):
    """Only parse a numeric signal for a donor in the prospectively chosen design."""
    values = {assay: {} for assay in design}
    for line in lines:
        _, assay, rest = line.split(b"\t", 2)
        if assay not in design:
            continue
        tokens = rest.split(b"\t")
        if len(tokens) != 5:
            raise ValueError("Malformed selected assay record")
        if not tokens[0].isdigit():
            raise ValueError("Invalid selected probe address")
        address = int(tokens[0])
        if address not in addresses or address in values[assay]:
            raise ValueError("Unexpected or duplicate probe address in selected assay")
        try:
            signal = float(tokens[1])
        except ValueError:
            raise ValueError("Missing or malformed signal in full WG-6 reference") from None
        if not math.isfinite(signal):
            raise ValueError("Nonfinite signal in full WG-6 reference")
        values[assay][address] = signal
    cells = defaultdict(lambda: {"case": [], "control": []})
    for assay, signals in values.items():
        if signals.keys() != addresses or target not in signals:
            raise ValueError("Incomplete WG-6 reference or missing target")
        fosb = signals[target]
        doubled_rank = sum(2 * (x < fosb) + (x == fosb) for x in signals.values())
        stratum, group = design[assay]
        cells[stratum][group].append(doubled_rank)
    return cells


def exact_contrast(cells, n_features):
    if len(cells) < 2:
        raise ValueError("At least two shared strata are required")
    observed = []
    permutations = []
    for groups in cells.values():
        cases, controls = groups["case"], groups["control"]
        nc, nt = len(cases), len(controls)
        if min(nc, nt) < 1:
            raise ValueError("Each stratum requires both disease groups")
        scores = cases + controls
        total = sum(scores)
        # A case assignment with score sum S has difference S/nc - (total-S)/nt.
        def difference(case_sum):
            return Fraction(case_sum * (nc + nt) - total * nc, nc * nt)

        observed.append(difference(sum(cases)))
        permutations.append([difference(sum(chosen)) for chosen in combinations(scores, nc)])
    # Every label assignment is equally weighted, including tied assignments.
    statistic = sum(observed)
    denominator = math.lcm(*(s.denominator for strata in permutations for s in strata))
    scaled = [[int(v * denominator) for v in strata] for strata in permutations]
    cutoff = int(abs(statistic) * denominator)
    extreme = sum(abs(sum(assignment)) >= cutoff for assignment in product(*scaled))
    assignments = math.prod(map(len, permutations))
    p = Fraction(extreme, assignments)
    points = Fraction(50, n_features)  # doubled midrank to percentile points
    effect = points * statistic / len(cells)
    leave_one_out = [points * (statistic - d) / (len(cells) - 1) for d in observed]
    decision = "inconclusive"
    if p <= Fraction(1, 20) and effect:
        decision = "supports_expected_direction" if effect > 0 else "supports_opposite_direction"
    return {
        "difference_percentile_points": float(effect),
        "extreme_assignments": extreme, "total_assignments": assignments,
        "two_sided_p_fraction": str(p), "two_sided_p": float(p),
        "leave_one_stratum_out_range_percentile_points": [float(min(leave_one_out)), float(max(leave_one_out))],
        "decision": decision,
    }


def analyze(directory, prefix=None):
    from inflate64 import Inflater

    head = directory / "archive-start.bin"
    header = Inflater().inflate(pinned(head, HEAD_SHA256)[48:])[:RAW_START]
    if RAW_COPY not in header or RAW_FK not in header:
        raise ValueError("Native signal column order or raw-to-clinical foreign key differs")
    clinical = read_metadata(head, directory / "archive-tail.bin")
    _, selected = eligible(clinical)
    native_ids = {row["id"].encode() for row in clinical["sample_analysis_result"]}
    manifests = {}
    for (file, digest), (name, address, probe) in zip(
        MANIFESTS.items(), [("WG6", TARGET, "ILMN_1730706"), ("HT12_R1", 7160239, "ILMN_1751607"),
                            ("HT12_R2", 7160239, "ILMN_1751607")], strict=True,
    ):
        mapping = manifest(directory / file, digest)
        if mapping.get(address) != probe:
            raise ValueError("Incorrect physical FOSB probe in manifest")
        manifests[name] = (set(mapping), address)
    wg6_addresses = manifests["WG6"][0]
    if len(wg6_addresses) != N_FEATURES:
        raise ValueError("Unexpected WG-6 manifest size")

    def blocks():
        return pg_table_blocks(raw_dump_bytes(input_bytes(directory, prefix)), RAW_TAG)

    status = {}
    census(raw_keys(blocks()), selected, native_ids, manifests, selected_status=status)
    design = shared_design(selected, status)
    cells = rank_signals(copy_lines(blocks()), design, wg6_addresses, TARGET)
    result = exact_contrast(cells, N_FEATURES)
    if result["total_assignments"] != 12_096:
        raise ValueError("Unexpected frozen label-assignment count")
    return {
        "cohort": "Finnish EndometDB v1", "platform": "HumanWG-6 v2",
        "probe": "ILMN_1730706", "address": TARGET, "reference_features": N_FEATURES,
        "eligible_single_assay_donors": len(selected), "chip_verified_donors": 32,
        "shared_strata": len(cells), "included": {"case": 12, "control": 13}, **result,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--prefix", type=Path, help="Single pinned ZIP prefix, as in the chip census")
    args = parser.parse_args()
    print(json.dumps(analyze(args.directory, args.prefix), indent=2, allow_nan=False))
