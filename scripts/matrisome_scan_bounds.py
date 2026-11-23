"""Identify possible same-block concordance from a published total win count."""

from functools import cache
from numbers import Integral


def feasible_within_wins(blocks, higher):
    """Return every same-block win count compatible with a strict score order.

    Each block supplies (cases, comparators). A win is a case scoring higher
    than a comparator. The total includes pairs from different blocks. Only
    cross-group ties would affect the calculation.
    """
    blocks = list(blocks)
    counts = tuple(n for cases, comparators in blocks for n in (cases, comparators))
    if not blocks or any(not isinstance(n, Integral) or isinstance(n, bool) or n < 0
                         for n in counts):
        raise ValueError("Supply blocks of nonnegative integer group counts")
    counts = tuple(map(int, counts))
    pairs = sum(counts[::2]) * sum(counts[1::2])
    if not isinstance(higher, Integral) or isinstance(higher, bool) or not 0 <= higher <= pairs:
        raise ValueError("Total wins must be an integer between zero and the total pair count")

    @cache
    def attainable(remaining, wins):
        comparators = sum(remaining[1::2])
        if wins < 0 or wins > comparators * sum(remaining[::2]):
            return set()
        if wins == 0:
            return {0}
        totals = set()
        # Choose the highest-ranked remaining observation. A case beats all
        # remaining comparators; a comparator contributes no new wins.
        for i, count in enumerate(remaining):
            if count == 0:
                continue
            is_case = i % 2 == 0
            same_block = remaining[i + 1] if is_case else 0
            lower = remaining[:i] + (count - 1,) + remaining[i + 1:]
            rest = attainable(lower, wins - (comparators if is_case else 0))
            totals.update(same_block + value for value in rest)
        return totals

    return sorted(attainable(counts, int(higher)))
