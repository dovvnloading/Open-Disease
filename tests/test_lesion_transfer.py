from itertools import product
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import lesion_transfer as transfer


def test_panel_is_discovery_only_and_ties_sort_by_entrez():
    ids = np.array([50, 40, 30, 20, 10, 60])
    diffs = np.array([[4, 4], [-3, -3], [4, 4], [-4, -4], [2, 2], [1, 1]])
    locked = transfer.frozen_panel(ids, diffs)
    assert [r["entrez"] for r in locked] == [30, 50, 20, 40, 10]
    assert [r["direction"] for r in locked] == [1, 1, -1, -1, 1]


def test_global_projection_preserves_donor_correlation():
    weights = np.array([2., -1., 1.])
    diffs = np.array([[2., -1., 4.], [-2., 4., -3.], [1., 3., 2.]])
    projected = transfer.projection(weights, diffs)
    np.testing.assert_allclose(projected, [1.75, -0.75, 3.25])
    tails = transfer.summarize(projected)
    expected = [sum(s*d for s, d in zip(sign, projected)) for sign in product([-1, 1], repeat=3)]
    assert tails["one_sided_p"] == pytest.approx(sum(x >= sum(projected) for x in expected)/8)


def test_gene_conjunction_rejects_opposite_signs_and_preserves_tail_resolution():
    ids = np.array([10, 20, 30])
    first = np.repeat([[1.], [1.], [0.]], 5, axis=1)
    second = np.repeat([[1.], [-1.], [0.]], 4, axis=1)
    rows = transfer.gene_results(ids, first, second)
    assert [r[-2] for r in rows] == [2/16, 1., 1.]
    assert [r[-1] for r in rows] == [6/16, 1., 1.]
