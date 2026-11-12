"""Aggregate cell counts at donor level, never treat cells as patients."""
from __future__ import annotations

import pandas as pd
from .schema import CELL, META, DONOR, require_columns, validate_metadata, validate_nonnegative_int, validate_pseudobulk


def aggregate_cells(cells: pd.DataFrame, *, max_dense_rows: int = 5_000_000) -> pd.DataFrame:
    """Small-data adapter: cell-level long counts -> donor pseudobulk.

    Expands explicit zeros within each (study, cell_type). For huge atlases,
    aggregate in Seurat/Scanpy and import donor-level counts instead.
    """
    require_columns(cells, CELL)
    validate_metadata(cells)
    cells = cells.copy()
    cells["count"] = validate_nonnegative_int(cells["count"], "count")
    for column in ("cell_id", "gene_id"):
        if cells[column].isna().any() or (cells[column].astype(str).str.strip() == "").any():
            raise ValueError(f"{column} must be present")
    cell_key = DONOR + ["cell_id"]
    if cells.duplicated(cell_key + ["gene_id"]).any():
        raise ValueError("Duplicate gene within the same cell")
    if (cells.groupby(cell_key)["cell_type"].nunique() > 1).any():
        raise ValueError("cell_id changes cell_type within a donor")
    sample_cols = META + ["cell_type"]
    dense = []
    for _, group in cells.groupby(["study_id", "cell_type"], sort=True):
        genes = sorted(group["gene_id"].unique())
        samples = group[sample_cols].drop_duplicates().reset_index(drop=True)
        if len(genes) * len(samples) > max_dense_rows:
            raise ValueError("Expansion too large: aggregate with Seurat/Scanpy, then export donor-level counts")
        grid = samples.merge(pd.DataFrame({"gene_id": genes}), how="cross")
        counts = group.groupby(sample_cols + ["gene_id"], as_index=False)["count"].sum()
        counts = counts.rename(columns={"count": "counts"})
        grid = grid.merge(counts, on=sample_cols + ["gene_id"], how="left")
        grid["counts"] = grid["counts"].fillna(0).astype("int64")
        cell_counts = group.groupby(sample_cols)["cell_id"].nunique().rename("n_cells").reset_index()
        grid = grid.merge(cell_counts, on=sample_cols, how="left")
        grid["library_size"] = grid.groupby(sample_cols)["counts"].transform("sum")
        dense.append(grid)
    return validate_pseudobulk(pd.concat(dense, ignore_index=True).sort_values(
        ["study_id", "donor_id", "cell_type", "gene_id"]).reset_index(drop=True))
