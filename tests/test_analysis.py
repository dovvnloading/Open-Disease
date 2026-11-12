import numpy as np
import pandas as pd
import pytest

from open_disease.analysis import analyze_pseudobulk, bh_adjust
from open_disease.synthetic import generate_synthetic


def test_bh_is_order_independent_and_bounded():
    p = np.array([0.03, 0.001, 0.8, 0.02, 1.0])
    q = bh_adjust(p)
    assert q[1] <= q[3] <= q[0] <= q[2] <= q[4]
    assert np.all((q >= p) & (q <= 1))
    assert np.allclose(q, bh_adjust(p[::-1])[::-1])


def test_phase_separation_and_patient_counts(tmp_path):
    pb, _ = generate_synthetic(tmp_path)
    result = analyze_pseudobulk(pd.read_csv(pb))
    assert len(result) == 2 * 2 * 50
    assert set(result.n_case) == {6}
    assert set(result.n_control) == {6}
    assert set(result.phase) == {"proliferative", "secretory"}
    g = result[(result.gene_id == "SYN000") & (result.phase == "proliferative")]
    assert (g.effect_log2_cpm > 1).all()
    assert (g.q_value < 0.10).all()


def test_missing_gene_is_not_mistaken_for_zero(tmp_path):
    pb, _ = generate_synthetic(tmp_path)
    data = pd.read_csv(pb).iloc[1:]
    with pytest.raises(ValueError, match="Incomplete gene matrix"):
        analyze_pseudobulk(data)


def test_insufficient_donors_produce_no_statistical_claim(tmp_path):
    pb, _ = generate_synthetic(tmp_path)
    result = analyze_pseudobulk(pd.read_csv(pb), min_donors_per_group=7)
    assert set(result.status) == {"insufficient_donors"}
    assert (result.q_value == 1).all()
