from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import verify_external_bulk_pathways as independent


def test_exact_transcript_sums_preserve_mathematical_ties():
    data = np.array([[Decimal("0.1")], [Decimal("0.2")], [Decimal("0.3")]], dtype=object)
    primary, alternate = independent.cohort_gene_scores(data, {1: [0, 1], 2: [2]}, {1, 2}, transcripts=True)
    assert primary == [{1: Fraction(50), 2: Fraction(50)}]
    assert alternate == [{1: Fraction(100, 3), 2: Fraction(250, 3)}]


def test_exact_interval_allows_unbounded_design_and_nonzero_tail_counts():
    # In this two-donor study every labeling has the same absolute shift-adjusted contrast.
    result = independent.compatibility_set([Fraction(3), Fraction(1)], [True, False])
    assert result == {"bounds_pp": [None, None], "assignments": 2, "minimum_extreme_for_inclusion": 1,
                      "tail_checks": {"zero": 2}}


def test_scan_unadjusted_loo_and_two_sided_ties_on_synthetic_scores():
    values = [Fraction(x) for x in (1, 1, 2, 2, 3, 3, 4)]
    rows = {"x": values, "GLOBAL": values}
    labels = [True, True, True, True, False, False, False]
    result = independent.summarize(rows, rows, rows, labels, exact=True)
    assert result["GLOBAL"]["extreme"] == 1
    assert result["GLOBAL"]["p"] == 1/35
    assert result["x"]["loo_low"] <= result["x"]["difference"] <= result["x"]["loo_high"]
