"""Donor-separated validation. A donor is (study_id, donor_id)."""
from __future__ import annotations

import numpy as np
import pandas as pd

DONOR = ["study_id", "donor_id"]
META = DONOR + ["diagnosis", "phase", "hormone_treatment"]
CELL = META + ["cell_type", "cell_id", "gene_id", "count"]
PB = META + ["cell_type", "gene_id", "counts", "library_size", "n_cells"]


def require_columns(df: pd.DataFrame, required: list[str]) -> None:
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))
    if df.empty:
        raise ValueError("Input table is empty")


def validate_metadata(df: pd.DataFrame) -> None:
    require_columns(df, META)
    for column in META + (["cell_type"] if "cell_type" in df else []):
        if df[column].isna().any() or (df[column].astype(str).str.strip() == "").any():
            raise ValueError(f"Null or blank {column} is not permitted")
    if not set(df["diagnosis"].unique()) <= {"case", "control"}:
        raise ValueError("diagnosis must be case or control")
    # Baseline restricted to the untreated atlas population.
    if not set(df["hormone_treatment"].unique()) <= {"none"}:
        raise ValueError("Baseline requires hormone_treatment=none; pre-specify other cohorts separately")
    varying = df.groupby(DONOR)[["diagnosis", "phase", "hormone_treatment"]].nunique(dropna=False)
    if (varying > 1).any().any():
        raise ValueError("Donor metadata conflict: diagnosis/phase/treatment changes within a donor")


def validate_nonnegative_int(series: pd.Series, name: str) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if numeric.isna().any() or not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError(f"{name} must contain finite numbers")
    if (numeric < 0).any() or (numeric != np.floor(numeric)).any():
        raise ValueError(f"{name} must contain nonnegative integers")
    return numeric.astype("int64")


def validate_pseudobulk(df: pd.DataFrame) -> pd.DataFrame:
    require_columns(df, PB)
    validate_metadata(df)
    out = df.copy()
    for name in ("counts", "library_size", "n_cells"):
        out[name] = validate_nonnegative_int(out[name], name)
    if (out["library_size"] == 0).any() or (out["n_cells"] == 0).any():
        raise ValueError("library_size and n_cells must be positive")
    if out["gene_id"].isna().any() or (out["gene_id"].astype(str).str.strip() == "").any():
        raise ValueError("gene_id must be present")
    if out.duplicated(DONOR + ["cell_type", "gene_id"]).any():
        raise ValueError("Duplicate donor/cell_type/gene record")
    sample = DONOR + ["cell_type"]
    if (out.groupby(sample)[["library_size", "n_cells"]].nunique() > 1).any().any():
        raise ValueError("Conflicting library_size or n_cells within pseudobulk sample")
    if (out.groupby(sample)["counts"].sum() >
            out.groupby(sample)["library_size"].first()).any():
        raise ValueError("Sum of counts exceeds library_size; use full-transcriptome library sizes")
    return out
