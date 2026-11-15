"""Synthetic regression cases for two easily confused genetic-panel tests."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from check_genetic_anchor_transfer import panel_check, signed_check
from restore_mumbai_ena_order import recover


def test_perfect_label_separation_can_be_nonspecific():
    genes = np.array([1, 2, 3, 10180])
    values = np.array([0., 1., 2., 3., 4., 5.])
    # The background carries exactly the same contrast as the candidate panel.
    matrix = np.column_stack([values, -values, values + 7, 2 * values])
    cases = np.array([False, False, False, True, True, True])
    result = panel_check(genes, matrix, cases, genes, [1, 2])
    assert result["self_contained"]["extreme"] == 2
    assert result["self_contained"]["assignments"] == 20
    assert result["competitive"]["value"] == pytest.approx(0, abs=1e-12)
    assert result["competitive"]["extreme"] == 20


def test_genetic_direction_and_deletion_keep_the_predefined_sign():
    genes = np.array([10, 10180])
    cases = np.array([True, True, False, False])
    matrix = np.array([[1., 10.], [2., 11.], [5., 2.], [6., 3.]])
    source = [{"symbol": "X", "gene_id": "10", "twas_z": "-5"},
              {"symbol": "RBM6", "gene_id": "10180", "twas_z": "6"}]
    result = signed_check(genes, matrix, cases, source)
    assert result["extreme"] == result["rbm6_only"]["extreme"] == result["excluding_rbm6"]["extreme"] == 1
    assert result["assignments"] == 6
    assert result["per_gene_effects_pp"] == {"X": -4., "RBM6": 8.}


def test_missing_or_uninformative_data_do_not_become_signed_evidence():
    case = np.array([True, False, True, False])
    source = [{"symbol": "X", "gene_id": "1", "twas_z": 4},
              {"symbol": "RBM6", "gene_id": "10180", "twas_z": -4}]
    with pytest.raises(ValueError, match="Missing"):
        signed_check(np.array([1]), np.zeros((4, 1)), case, source)
    with pytest.raises(ValueError, match="Constant"):
        signed_check(np.array([1, 10180]), np.ones((4, 2)), case, source)
    with pytest.raises(ValueError, match="Unexpected"):
        recover(b"run_accession\n")
