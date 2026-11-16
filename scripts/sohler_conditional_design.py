"""Enumerate paired tissue-label assignments that preserve each scanner-day's counts."""

from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
import re

from audit_sohler_peritoneum import TISSUES


def verify_dat_filename(prefix, filename):
    """Match the embedded scanner image name to its versioned SDRF CEL filename."""
    suffix = ".v4.cel"
    if Path(filename).name != filename or not filename.lower().endswith(suffix):
        raise ValueError("Unexpected versioned SDRF CEL basename")
    dat = [line for line in prefix[24:].splitlines() if line.startswith(b"DatHeader=")]
    if len(dat) != 1:
        raise ValueError("Expected one embedded DAT header")
    first, separator, _ = dat[0].removeprefix(b"DatHeader=").partition(b"\x14")
    # The scan image basename ends before the first field of the DAT geometry.
    named = re.match(rb"^\[\d+\.\.\d+\]\s+([^\s:]+):CLS=\d+\s+", first)
    if (not separator or not named or
            named.group(1).decode("ascii").casefold() != filename[:-len(suffix)].casefold()):
        raise ValueError("Embedded DAT image basename disagrees with SDRF CEL filename")


def conditional_design(records):
    """Use scanner times and tissue metadata only. The unmodified assignment is included."""
    lesion, peritoneum, _ = TISSUES
    by_donor = defaultdict(dict)
    for donor, tissue, stamp, scanner in records:
        if tissue in by_donor[donor]:
            raise ValueError("Repeated donor/tissue in conditional design")
        by_donor[donor][tissue] = (stamp, scanner)
    pairs = [(t[lesion], t[peritoneum]) for t in by_donor.values() if lesion in t and peritoneum in t]
    if len(pairs) > 16:
        raise ValueError("Conditional audit is limited to 16 donor pairs")

    blocks = []
    chronological = []
    timeline = defaultdict(list)
    for (l_time, l_scanner), (p_time, p_scanner) in pairs:
        left, right = (l_time.date(), l_scanner), (p_time.date(), p_scanner)
        blocks.append((left, right))
        # Clock times are comparable within a physical scanner only.
        if l_scanner == p_scanner and l_time == p_time:
            raise ValueError("Tied timestamp on one scanner")
        chronological.append((l_time < p_time) if l_scanner == p_scanner else None)
        timeline[left].append((l_time, lesion))
        timeline[right].append((p_time, peritoneum))

    switches = Counter()
    for slots in timeline.values():
        ordered = [tissue for _, tissue in sorted(slots)]
        if len(set(ordered)) == 2:
            switches[sum(a != b for a, b in zip(ordered, ordered[1:]))] += 1

    observed = Counter(left for left, _ in blocks)
    cross_day = {i for i, (left, right) in enumerate(blocks) if left != right}
    movable = set()
    order_movable = set()
    admissible = with_order = 0
    for swaps in product((False, True), repeat=len(pairs)):
        assigned = Counter(right if swap else left for (left, right), swap in zip(blocks, swaps))
        if assigned != observed:
            continue
        admissible += 1
        changed = {i for i, swap in enumerate(swaps) if swap}
        movable.update(changed)
        # Among same-scanner pairs, flipping an early lesion makes it late, and vice versa.
        order_change = sum(1 if chronological[i] else -1 for i in changed if chronological[i] is not None)
        if order_change == 0:
            with_order += 1
            order_movable.update(changed)
    return {
        "paired_donors": len(pairs), "unconditional_assignments": 2**len(pairs),
        "scanner_day_margin_preserving_assignments": admissible,
        "movable_pairs": len(movable), "movable_cross_scanner_day_pairs": len(movable & cross_day),
        "all_pairs_swapped_preserves_scanner_days": Counter(right for _, right in blocks) == observed,
        "same_scanner_lesion_scanned_first": chronological.count(True),
        "same_scanner_peritoneum_scanned_first": chronological.count(False),
        "different_scanner_order_unassigned": chronological.count(None),
        "additionally_preserving_same_scanner_lesion_first_count": {
            "assignments": with_order, "movable_pairs": len(order_movable),
            "movable_cross_scanner_day_pairs": len(order_movable & cross_day),
        },
        "mixed_scanner_days_by_number_of_tissue_switches_among_paired_arrays": dict(sorted(switches.items())),
    }
