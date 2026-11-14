"""Synthetic source and physical-address alignment checks, without study records."""
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import genome_endometdb as fin


def test_annotation_preserves_reference_only_rows_and_multiprobe_genes():
    keys = ("Array_Address_Id", "Species", "Source", "Entrez_Gene_ID")
    raw = [
        ["001", "Homo sapiens", "RefSeq", "22"], ["002", "Homo sapiens", "RefSeq", "22"],
        ["003", "Homo sapiens", "RefSeq", "2 /// 3"], ["004", "Homo sapiens", "Unigene", "5"],
        ["005", "Mus musculus", "RefSeq", "22"], ["006", "Homo sapiens", "RefSeq", "0"],
        ["007", "Homo sapiens", "ILMN_Controls", ""], ["[Controls]", "", "", ""],
    ]
    rows = [dict(zip(keys, row)) for row in raw]
    assert fin.annotation_entrez(rows) == {1: 22, 2: 22, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0}
    with pytest.raises(ValueError, match="Duplicate"):
        fin.annotation_entrez(rows + [rows[0]])


def raw_inputs():
    design = {b"toy_b": (("2", "phase_b"), "case"), b"toy_a": (("1", "phase_a"), "control")}
    mapping = {20: 51, 10: 0, 30: 51}
    lines = [b"0\tignored\tmalformed opaque signal"] + [
        b"0\t" + sample + b"\t" + str(address).encode() + b"\t" + str(signal).encode() + b"\t.\t.\t."
        for sample, values in [(b"toy_b", [-2, -2, 9]), (b"toy_a", [5, 6, 5])]
        for address, signal in zip([30, 20, 10], values)
    ]
    return lines, design, mapping


def test_private_row_order_and_full_physical_universe_are_aligned():
    lines, design, mapping = raw_inputs()
    matrix, genes, cases, strata = fin.selected_signals(reversed(lines), design, mapping)
    np.testing.assert_array_equal(matrix, [[5, 6, 5], [9, -2, -2]])
    np.testing.assert_array_equal(genes, [0, 51, 51])
    np.testing.assert_array_equal(cases, [False, True])
    np.testing.assert_array_equal(strata, ["1-phase_a", "2-phase_b"])


@pytest.mark.parametrize(("change", "error"), [
    ("missing", "Incomplete"), ("duplicate", "duplicated"), ("address", "physical address"),
    ("malformed", "Malformed"), ("bad_value", "malformed"), ("infinite", "Nonfinite"),
])
def test_selected_assays_require_unique_complete_finite_physical_values(change, error):
    lines, design, mapping = raw_inputs()
    if change == "missing":
        lines.pop(1)
    elif change == "duplicate":
        lines.append(lines[1])
    else:
        lines[1] = {
            "address": b"0\ttoy_b\t99\t2\t.\t.\t.",
            "malformed": b"0\ttoy_b\t30\t2",
            "bad_value": b"0\ttoy_b\t30\tNA\t.\t.\t.",
            "infinite": b"0\ttoy_b\t30\tinf\t.\t.\t.",
        }[change]
    with pytest.raises(ValueError, match=error):
        fin.selected_signals(lines, design, mapping)


def test_load_finland_returns_aligned_gene_percentiles(monkeypatch):
    lines, design, mapping = raw_inputs()
    monkeypatch.setattr(fin, "load_wg6_annotation", lambda _: mapping)
    monkeypatch.setattr(fin, "verified_wg6_source", lambda _: (design, mapping.keys(), lines))
    genes, scores, cases, strata = fin.load_finland(Path("."))
    np.testing.assert_array_equal(genes, [51])
    np.testing.assert_allclose(scores[:, 0], [(5 / 6 + 1 / 3) * 50, 100 / 3])
    np.testing.assert_array_equal(cases, [False, True])
    np.testing.assert_array_equal(strata, ["1-phase_a", "2-phase_b"])


def test_alternative_null_variance_and_maxT_preserve_synchronized_labels():
    from verify_finnish_genome_screen import recompute, sorted_gene_scores
    genes, ranks = sorted_gene_scores(np.array([[1., -2., 1., 0.]]), np.array([1, 0, 3, 1]))
    np.testing.assert_array_equal(genes, [1, 3])
    np.testing.assert_allclose(ranks, [[(75 + 37.5) / 2, 75]])
    # All four within-pair permutations. Two orthogonal effects share max |z|=sqrt(2).
    x = np.array([[1, 0, 8], [0, 1, 8], [1, 1, 3], [0, 0, 3]], dtype=float)
    cases = np.array([True, False, True, False])
    strata = np.array(["first", "first", "second", "second"])
    effects, z, raw_count, max_count, n = recompute(x, cases, strata, block=1)
    np.testing.assert_allclose(effects, [1, 0, 0])
    np.testing.assert_allclose(z, [2 ** .5, 0, 0])
    np.testing.assert_array_equal(raw_count, [2, 4, 4])
    np.testing.assert_array_equal(max_count, [4, 4, 4])
    assert n == 4


def test_phase_diagnostic_weights_runs_equally_in_each_phase():
    from describe_cross_cohort_phase import phase_effects
    values = np.array([[3.], [1.], [4.], [0.], [0.], [8.], [12.]])
    cases = np.array([True, False, True, False, False, True, False])
    strata = np.array(["1-proliferative"] * 2 + ["2-proliferative"] * 3 + ["1-secretory"] * 2)
    result = phase_effects(values, cases, strata)
    assert result["proliferative"][0] == 2
    assert result["secretory"][0] == 1
    np.testing.assert_allclose(result["proliferative"][1], [3.])
    np.testing.assert_allclose(result["secretory"][1], [-4.])
    cases[1] = True
    with pytest.raises(ValueError, match="overlap"):
        phase_effects(values, cases, strata)
