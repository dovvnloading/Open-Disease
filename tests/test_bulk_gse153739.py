from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from bulk_gse153739 import fixed_value, gene_scores


def test_two_stage_ranking_handles_summed_decimal_ties_and_unmapped_transcripts():
    ids = np.array([7, 7, 91, 0])
    expression = np.array([[fixed_value(x)] for x in ("0.1", "0.2", "0.3", "0.0000000001")])
    genes, primary, alternate = gene_scores(ids, expression)
    assert genes.tolist() == [7, 91]
    assert primary.tolist() == [[50, 50]]  # exact equality of 0.1+0.2 and 0.3
    assert alternate.tolist() == [[50, 87.5]]
    assert fixed_value("1e-10") == fixed_value("0.0000000001") == 1


@pytest.mark.parametrize("value", ["nan", "Infinity", "-0.01", "1.1e-10", "undetected"])
def test_invalid_transcript_values_fail(value):
    with pytest.raises(ValueError, match="FPKM"):
        fixed_value(value)
