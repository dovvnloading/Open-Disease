"""Synthetic coverage for physical-array ranking and the Rome validation cohort."""
import csv
import gzip
from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import genome_rome
from open_disease.microarray_ranks import gene_percentiles


def write_table(path, rows):
    with gzip.open(path, "wt", newline="") as fh:
        csv.writer(fh, delimiter="\t", lineterminator="\n").writerows(rows)
    return path


def test_ranks_include_unmapped_and_collapse_redundant_probes():
    # Feature one is unmatched but ties with a biological measurement.
    data = [[3, 3, 9, -1], [2, 8, 0, 8]]
    genes, scores = gene_percentiles(data, [0, 12, 5, 12])
    assert list(genes) == [5, 12]
    np.testing.assert_array_equal(scores, [[87.5, 31.25], [12.5, 75.0]])
    _, scaled = gene_percentiles(np.exp(data), [0, 12, 5, 12])
    np.testing.assert_array_equal(scores, scaled)


@pytest.mark.parametrize("signals,ids", [
    ([[np.nan, 0]], [1, 2]), ([[1, 2]], [1]), ([[1, 2]], [-1, 0]),
    ([[1, 2]], [0, 0]), ([[1, 2]], [1.0, 2.0]),
])
def test_rejects_incomplete_or_invalid_gene_matrix(signals, ids):
    with pytest.raises(ValueError):
        gene_percentiles(signals, ids)


def test_annotation_excludes_controls_and_ambiguous_entrez(tmp_path):
    rows = [["!platform_table_begin"], ["ID", "Gene ID", "Gene symbol"],
            ["p1", "11", "A"], ["p2", "12 /// 13", "B /// C"], ["p3", "", ""],
            ["AFFX-p4", "14", "D"], ["p5", "15", "F /// G"], ["!platform_table_end"]]
    p = write_table(tmp_path / "ann.gz", rows)
    assert genome_rome.read_annotation(p) == {"p1": 11, "p2": 0, "p3": 0, "AFFX-p4": 0, "p5": 0}
    rows.insert(-1, rows[2])
    with pytest.raises(ValueError, match="repeated"):
        genome_rome.read_annotation(write_table(p, rows))


def matrix_rows():
    groups = ["control"] * 6 + ["lesion", "case"] * 8
    ids = [f"GSM9900{i:03}" for i in range(22)]
    descriptions = {"control": "Normal (control)", "case": "Normal (eutopic)",
                    "lesion": "Pathological (ectopic)"}
    return [
        ["!Series_geo_accession", "GSE25628"], ["!Series_platform_id", "GPL571"],
        ["!Sample_geo_accession", *ids], ["!Sample_platform_id", *(["GPL571"] * 22)],
        ["!Sample_data_processing", *(["RMA + quantile normalization"] * 22)],
        ["!Sample_characteristics_ch1", *["donor: "+("healthy donor" if g == "control" else "patient")
                                         for g in groups]],
        ["!Sample_characteristics_ch1", *(["tissue: endometrium"] * 22)],
        ["!Sample_characteristics_ch1", *["disease state: " + descriptions[g] for g in groups]],
        ["!series_matrix_table_begin"], ["ID_REF", *ids],
        ["partial", "null", *([1] * 21)],
        ["retained", *["ignored" if g == "lesion" else 2 for g in groups]],
        ["reference", *["ignored" if g == "lesion" else 5 for g in groups]],
        ["!series_matrix_table_end"],
    ]


def test_matrix_uses_uterine_columns_and_common_feature_reference(tmp_path):
    rows = matrix_rows()
    p = write_table(tmp_path / "matrix.gz", rows)
    mapping = {"partial": 0, "retained": 45, "reference": 0}
    x, genes, cases = genome_rome.read_matrix(p, mapping, min_complete=2)
    assert x.shape == (14, 2)
    np.testing.assert_array_equal(x, np.tile([2, 5], (14, 1)))
    assert list(genes) == [45, 0] and list(cases) == [False] * 6 + [True] * 8
    with pytest.raises(ValueError, match="insufficient"):
        genome_rome.read_matrix(p, mapping, min_complete=3)
    rows[-3][1] = "corrupt"
    with pytest.raises(ValueError, match="Malformed"):
        genome_rome.read_matrix(write_table(p, rows), mapping, min_complete=2)


def test_matrix_and_annotation_require_one_to_one_feature_rows(tmp_path):
    p = write_table(tmp_path / "matrix.gz", matrix_rows())
    with pytest.raises(ValueError, match="Incomplete"):
        genome_rome.read_matrix(p, {"partial": 0, "retained": 2, "reference": 0, "absent": 3}, 2)


def test_incomplete_eligible_probe_cannot_change_gene_composition(tmp_path):
    p = write_table(tmp_path / "matrix.gz", matrix_rows())
    mapping = {"partial": 45, "retained": 45, "reference": 0}
    with pytest.raises(ValueError, match="Incomplete eligible"):
        genome_rome.read_matrix(p, mapping, min_complete=2)


def test_validation_cross_checker_uses_only_locked_genes_and_projected_members():
    from verify_rome_pathway_screen import membership_rows

    members = list(range(31, 46))
    candidates = [{"family": "gene", "id": "33", "direction": "lower"},
                  {"family": "pathway", "id": "R-HSA-42", "direction": "higher",
                   "members": len(members), "entrez_ids": members}]
    expected = [("gene", "33", "lower", 33)]
    expected.extend(("pathway", "R-HSA-42", "higher", gene) for gene in members)
    assert list(membership_rows(candidates)) == expected
    with pytest.raises(ValueError, match="Duplicate"):
        list(membership_rows(candidates * 2))
