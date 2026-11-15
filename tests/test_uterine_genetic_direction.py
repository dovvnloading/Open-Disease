import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from uterine_genetic_direction import score_cohort


def test_equal_band_weights_and_signs():
    genes = np.array([10, 20, 30])
    cases = np.array([True, False, True, False])
    groups = np.array(["a", "a", "b", "b"])
    panel = [{"entrez": gene, "symbol": str(gene), "band": band, "direction": sign}
             for gene, band, sign in [(10, "A", 1), (20, "A", -1), (30, "B", 1)]]
    matrix = np.array([[23, 18, 31], [20, 20, 30], [23, 18, 31], [20, 20, 30]])
    result = score_cohort(panel, genes, matrix, cases, groups)
    assert result["case_minus_control_pp"] == pytest.approx(1.75)
    assert result["omit_one_band_pp"] == {"A": 1., "B": 2.5}
    assert (result["positive_tail"], result["two_sided_tail"], result["assignments"]) == (1, 2, 4)
