from copy import deepcopy
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_swiss_glandular as study


def test_batch_permutation_ties_weighting_and_omission():
    case = np.array([True, False, True, False, False])
    batch = np.array([1, 1, 2, 2, 2])
    # Within-batch case-minus-control effects are +4 and -2. Equal weighting
    # gives +1. The six null values are +1,+2.5,+2.5,-3,-1.5,-1.5.
    x = np.tile(np.array([7., 3., 0., 2., 2.])[:, None], (1, 6))
    result = study.summarize(case, batch, x, 2*x)
    row = result["endpoints"]["GLOBAL"]
    assert result["assignments"] == 6
    assert row["difference"] == pytest.approx(1)
    assert row["extreme"] == 6 and row["p"] == 1
    assert row["cell_first"] == pytest.approx(2)
    assert (row["batch_omit_low"], row["batch_omit_high"]) == (-2, 4)
    assert not any(item["supported"] for item in result["endpoints"].values())


def test_gate_blocks_even_a_strong_component():
    case = np.array([True, False] * 9)
    batch = np.repeat(np.arange(9), 2)
    contrast = 25 + 20*case
    x = np.column_stack([np.full(18, 30), contrast, 75-contrast, *([np.full(18, 25)]*3)])
    result = study.summarize(case, batch, x, x)
    row = result["endpoints"][study.ENDPOINTS[1]]
    assert row["holm"] < .05 and not row["supported"]
    assert result["endpoints"]["GLOBAL"]["p"] == 1
    result["cells"] = 1000
    other = deepcopy(result)
    study.compare_independent(result, other)
    other["endpoints"]["GLOBAL"]["batch_omit_high"] += 0.001
    with pytest.raises(ValueError, match="batch_omit_high"):
        study.compare_independent(result, other)


@pytest.mark.parametrize("fault", ["missing", "duplicate", "label", "global", "nan"])
def test_private_donor_scores_reject_incomplete_or_inconsistent_input(tmp_path, fault):
    path = tmp_path / "scores.tsv"
    rows = [["1", "a", "0", key, "35", "35"] for key in study.ENDPOINTS]
    if fault == "missing":
        rows.pop()
    elif fault == "duplicate":
        rows.append(rows[0])
    elif fault == "label":
        rows[1][2] = "1"
    elif fault == "global":
        rows[0][-1] = "36"
    else:
        rows[0][-1] = "nan"
    header = ["donor_index", "capture", "case", "endpoint", "primary", "cell_first"]
    path.write_text("\n".join("\t".join(row) for row in [header, *rows]) + "\n")
    with pytest.raises(ValueError):
        study.read_scores(path)
