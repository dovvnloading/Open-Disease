import pandas as pd
import pytest

from open_disease.pseudobulk import aggregate_cells
from open_disease.schema import validate_pseudobulk


def toy_cells():
    return pd.DataFrame([
        dict(study_id="A", donor_id=donor, diagnosis=diagnosis,
             phase="secretory", hormone_treatment="none", cell_type="epithelial",
             cell_id=f"{donor}-c1", gene_id=gene, count=value)
        for donor, diagnosis, amounts in [("d1", "case", (5, 0)), ("d2", "control", (1, 3))]
        for gene, value in zip(("GENE1", "GENE2"), amounts)
        if value != 0
    ])


def test_pseudobulk_fills_missing_zero_counts_without_inventing_donors():
    pb = aggregate_cells(toy_cells())
    zero = pb[(pb.donor_id == "d1") & (pb.gene_id == "GENE2")].iloc[0]
    assert zero.counts == 0
    assert zero.library_size == 5
    assert len(pb) == 4
    assert pb.groupby("donor_id").n_cells.first().to_dict() == {"d1": 1, "d2": 1}


def test_donor_label_conflicts_are_rejected():
    cells = toy_cells()
    cells.loc[1, "diagnosis"] = "case"
    with pytest.raises(ValueError, match="Donor metadata conflict"):
        aggregate_cells(cells)


def test_cell_gene_duplicates_are_rejected():
    cells = pd.concat([toy_cells(), toy_cells().iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="Duplicate gene"):
        aggregate_cells(cells)


def test_pseudobulk_inconsistent_library_is_rejected():
    pb = aggregate_cells(toy_cells())
    pb.loc[0, "library_size"] = 100
    with pytest.raises(ValueError, match="Conflicting library_size"):
        validate_pseudobulk(pb)


def test_treatment_outside_baseline_requires_new_protocol():
    pb = aggregate_cells(toy_cells())
    pb.loc[0, "hormone_treatment"] = "yes"
    with pytest.raises(ValueError, match="hormone_treatment"):
        validate_pseudobulk(pb)
