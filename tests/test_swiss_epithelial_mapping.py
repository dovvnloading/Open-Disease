from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import swiss_epithelial_mapping as hg
from swiss_prior_cohorts import pathway_effects


def entry(hgnc, symbol, gene, ensembl="", **kw):
    return dict(hgnc_id=hgnc, symbol=symbol, entrez_id=gene, ensembl_gene_id=ensembl,
                status="Approved", **kw)


def test_only_exact_unambiguous_matches_survive_in_original_row_order():
    records = [
        entry("HGNC:1", "FIRST", "10", "ENSG00000000001", alias_symbol="ALIAS"),
        entry("HGNC:2", "SECOND", "20", "ENSG00000000002"),
        entry("HGNC:3", "THIRD", "30", "ENSG00000000003"),
        entry("HGNC:4", "MULTI", "40"), entry("HGNC:5", "ANOTHER", "40"),
        entry("HGNC:6", "AMBIG", "60"), entry("HGNC:7", "AMBIG", "70"),
        entry("HGNC:8", "MIXED", "80|90"), entry("HGNC:9", "SPLIT", "90", "ENSG00000000009|X"),
        {**entry("HGNC:10", "OLD", "100"), "status": "Entry Withdrawn"},
    ]
    labels = ["FIRST", "ENSG00000000002.7", "THIRD", "ENSG00000000003", "MULTI", "AMBIG",
              "ALIAS", "first", "MIXED", "ENSG00000000009.2", "OLD", "ENSG00000000002.extra"]
    mapped, qc = hg.crosswalk(records, labels)
    assert mapped == {"10": ("FIRST", 1), "20": ("ENSG00000000002.7", 2)}
    assert qc == {"source_rows": 12, "no_match": 6, "ambiguous_label": 1,
                  "multiple_hgnc_entries": 1, "multiple_feature_rows": 2,
                  "matched_symbol": 1, "matched_ensembl": 1}
    assert hg.crosswalk(list(reversed(records)), labels) == (mapped, qc)


@pytest.mark.parametrize("labels", [["A", "A"], ["A", ""]])
def test_bad_axis_is_rejected(labels):
    with pytest.raises(ValueError, match="nonempty and unique"):
        hg.crosswalk([], labels)


def test_membership_projection_keeps_missing_and_cross_path_shared_members(tmp_path):
    mapping = {"11": ("X", 2)}
    g = {"family": "gene", "id": "11", "members": 1}
    p = {"family": "pathway", "id": "R-HSA-1", "entrez_ids": [11, 12], "members": 2}
    rows = hg.candidate_projection([g, p], mapping)
    table = tmp_path / "map.tsv"
    hg.write_projection(table, rows)
    assert table.read_text().splitlines() == [
        "family\tid\tentrez_id\tfeature_label\trow_index",
        "gene\t11\t11\tX\t2", "pathway\tR-HSA-1\t11\tX\t2", "pathway\tR-HSA-1\t12\t\t",
    ]
    assert hg.source_receipt({}, rows, table)["distinct_mapped_genes"] == 1
    with pytest.raises(ValueError, match="Duplicate"):
        hg.candidate_projection([p, p], mapping)
    with pytest.raises(ValueError, match="membership"):
        hg.candidate_projection([{**p, "entrez_ids": [11, 11]}], mapping)


def test_verified_mapping_refuses_wrong_reference_axis_or_lock(tmp_path, monkeypatch):
    (tmp_path / hg.HGNC_FILE).write_text("unexpected")
    with pytest.raises(ValueError, match="HGNC"):
        hg.load_projection(tmp_path)
    monkeypatch.setattr(hg, "HGNC_SIZE", (tmp_path / hg.HGNC_FILE).stat().st_size)
    monkeypatch.setattr(hg, "HGNC_SHA256", hg.sha256_file(tmp_path / hg.HGNC_FILE))
    labels = tmp_path / "external-feature-labels.txt"
    labels.write_text("X\n")
    with pytest.raises(ValueError, match="feature label"):
        hg.load_projection(tmp_path, labels_path=labels)
    monkeypatch.setattr(hg, "LABEL_SHA256", hg.sha256_file(labels))
    with pytest.raises(ValueError, match="nominations"):
        hg.load_projection(tmp_path, labels, labels_path=labels)


def test_prior_membership_comparison_preserves_equal_stratum_weight_and_sign():
    # Different group proportions ensure that pooling across strata would fail this test.
    genes = np.array([11, 12])
    scores = np.array([[7, 3], [3, 7], [0, 3], [0, 3], [2, 2]], dtype=float)
    cases = np.array([1, 0, 1, 1, 0], dtype=bool)
    strata = np.array(["x", "x", "y", "y", "y"])
    effects = pathway_effects((genes, scores, cases, strata), [[11, 12], [11]])
    np.testing.assert_allclose(effects, [-.25, 1])
    with pytest.raises(ValueError, match="absent"):
        pathway_effects((genes, scores, cases, strata), [[13]])
