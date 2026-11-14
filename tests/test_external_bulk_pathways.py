from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from external_bulk_pathways import evaluate, project


def test_common_membership_excludes_unmeasured_rows_with_original_denominator():
    full = {"p": list(range(20))}
    swiss = {"p": {0, 1, 2, 18, 19}}
    assert project(full, swiss, set(range(16))) == ({"p": list(range(16))}, {"p": [0, 1, 2]})
    with pytest.raises(ValueError, match="coverage|Insufficient"):
        project(full, swiss, set(range(15)))


def test_allocation_resolution_and_equal_path_weighting():
    case = np.array([True] * 4 + [False] * 3)
    x = np.column_stack([case * 40., case * 20., case * 0.])
    primary = {"a": [1], "b": [2, 3]}
    dataset = (np.array([1, 2, 3]), x, case, x)
    rows, interval = evaluate(dataset, primary, primary, inference=True)
    assert rows["GLOBAL"]["difference"] == 25
    assert rows["GLOBAL"]["extreme"] == 1 and rows["GLOBAL"]["p"] == 1 / 35
    assert interval["assignments"] == 35 and interval["tail_checks"]["zero"] == 1
    assert rows["a"]["holm"] == 2 / 35
    rows, interval = evaluate(dataset, primary, primary, inference=False)
    assert "p" not in rows["GLOBAL"] and interval is None
