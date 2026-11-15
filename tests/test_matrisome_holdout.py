from pathlib import Path

import numpy as np
import pytest
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from audit_matrisome_holdout import read_keyed
from matrisome_transfer_benchmark import pair_comparison, training_reference


def test_count_ties_as_half_pairs():
    row = pair_comparison(np.array([1, 2, 2, 0]), np.array([True, True, False, False]))
    assert row == {"cases": 2, "comparators": 2, "higher": 2, "ties": 1, "pairs": 4, "fraction": 0.625}


def test_reference_averages_phases_equally():
    phases = np.array(["p", "p", "p", "s", "s", "m", "m"])
    cases = np.array([True, True, False, True, False, True, False])
    data = np.array([[2, 1], [4, 1], [0, 1], [0, 3], [0, 1], [0, 2], [0, 0]])
    assert np.allclose(training_reference(data, cases, phases), [0.6, 0.8])
    with pytest.raises(ValueError, match="phase overlap"):
        training_reference(data[:-1], cases[:-1], phases[:-1])


def test_duplicate_author_keys_cannot_silently_overwrite(tmp_path: Path):
    file = tmp_path / "example.tsv"
    file.write_text("sample_name\tcondition\nx\t0\nx\t1\n")
    with pytest.raises(ValueError, match="Repeated"):
        read_keyed(file)
