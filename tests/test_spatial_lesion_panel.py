from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from spatial_lesion_panel import CELL_TYPES, TISSUES, summarize_regions


def test_region_replicates_and_five_genes_do_not_inflate_pair_count():
    records = []
    for index in range(5):
        for cell in CELL_TYPES:
            for side, tissue in enumerate(TISSUES):
                for technical in [0.0, 2.0]:
                    vector = np.repeat(technical + side * (index + 1), 5)
                    passed = not (index == 4 and side == 1)
                    records.append(
                        ((str(index), tissue, cell), passed, vector, np.ones(5, dtype=bool))
                    )
    symbols = list("ABCDE")
    for selected, pairs, effect in [(False, 5, 3.0), (True, 4, 2.5)]:
        summary = summarize_regions(records, symbols, [1] * 5, selected)
        for result in summary.values():
            assert result["pairs"] == pairs and result["positive_pairs"] == pairs
            assert result["effect_pp"] == pytest.approx(effect)
            assert result["one_sided_p"] == 1 / 2**pairs
            assert result["holm_p"] == 3 / 2**pairs
