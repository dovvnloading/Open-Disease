"""Audit donor-level expression cohort without exposing patient records.

Example: python scripts/audit_cohort.py --pseudobulk data/processed/study.csv --out results/audit.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
from open_disease.qc import cohort_audit


def main() -> None:
    parser = argparse.ArgumentParser(description="Structural research eligibility audit, NOT a clinical quality approval")
    parser.add_argument("--pseudobulk", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-donors-per-group", type=int, default=3)
    parser.add_argument("--min-cells-per-sample", type=int, default=20)
    parser.add_argument("--min-genes", type=int, default=20)
    args = parser.parse_args()
    report = cohort_audit(pd.read_csv(args.pseudobulk),
                          min_donors_per_group=args.min_donors_per_group,
                          min_cells_per_sample=args.min_cells_per_sample,
                          min_genes=args.min_genes)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out, index=False)
    print(report.to_string(index=False))
    print("Review fields remain REQUIRED even when baseline_eligible=true.")
    if not report["baseline_eligible"].all():
        raise SystemExit("One or more cohort strata failed the structural baseline audit")


if __name__ == "__main__":
    main()
