"""Conditional Rome comparison and four technical-control contrasts."""

import argparse
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
import re

import numpy as np

from analyze_gse25628_fosb import checked_groups
from cross_cohort_discovery import assignment_design, holm
from rome_scan import audit
from transcriptome_transfer import check_file, study_inputs, unique_gene_symbols
from verify_melbourne_s20 import HEADER, WORKBOOK


def controls(directory):
    matrix = directory / "GSE25628_series_matrix.txt.gz"
    header = defaultdict(list)
    probes = {}
    with gzip.open(matrix, "rt", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if row == ["!series_matrix_table_begin"]:
                break
            if row and row[0].startswith(("!Series_", "!Sample_")):
                header[row[0]].append(row[1:])
        ids, groups = checked_groups(header)
        if next(reader) != ["ID_REF", *ids]:
            raise ValueError("Misaligned Rome matrix")
        selected = np.array(groups) != "lesion"
        for row in reader:
            if row == ["!series_matrix_table_end"]:
                break
            if row[0].startswith("AFFX-"):
                probes[row[0]] = np.array(row[1:], dtype=float)[selected]
    results = {}
    for key, label in [("ACTB", "AFFX-HSAC07/X00351"), ("GAPDH", "AFFX-HUMGAPDH/M33197")]:
        results[key + "_log2_3minus5"] = probes[label + "_3_at"] - probes[label + "_5_at"]
    for name, species in [
        ("hybridization", ["Ec-bioB", "Ec-bioC", "Ec-bioD", "P1-cre"]),
        ("RNA_preparation", ["Bs-dap", "Bs-lys", "Bs-phe", "Bs-thr"]),
    ]:
        means = []
        for tag in species:
            subset = [
                val
                for k, val in probes.items()
                if re.fullmatch(r"AFFX-r2-" + tag + r"-(3|5|M)(_s)?_at", k)
            ]
            if len(subset) not in (2, 3):
                raise ValueError("Missing expected exogenous probe sets")
            means.append(np.mean(subset, axis=0))
        results[name + "_spike_mean_log2"] = np.mean(means, axis=0)
    if not np.isfinite(np.stack(list(results.values()))).all():
        raise ValueError("Nonfinite technical-control values")
    return results


def blocked_effect(scores, case, days):
    mixed = [d for d in set(days) if len(np.unique(case[days == d])) == 2]
    use = np.isin(days, mixed)
    if not np.any(use):
        raise ValueError("No mixed scan dates")
    obs, weights, _ = assignment_design(case[use], days[use])
    score = scores[use]
    observed = float(obs @ score)
    extreme = int((np.abs(weights @ score) >= abs(observed) - 1e-11).sum())
    return dict(
        effect=observed,
        extremes=extreme,
        assignments=len(weights),
        p_two_sided=extreme / len(weights),
        cases=int(case[use].sum()),
        controls=int((~case[use]).sum()),
        within_day_effects=sorted(
            float(scores[case & (days == d)].mean() - scores[~case & (days == d)].mean())
            for d in mixed
        ),
    )


def cycle_genes(private):
    from openpyxl import load_workbook

    path = private / "fung-2018/Fung_supplementary_tables.xlsx"
    check_file(path, WORKBOOK)
    symbols = unique_gene_symbols(private / "gse25628")
    book = load_workbook(path, read_only=True, data_only=True)
    genes = set()
    try:
        for table in ("TableS2", "TableS3", "TableS4", "TableS5"):
            rows = iter(book[table].values)
            title = next(rows)[0]
            if not title.startswith("Table S" + table.removeprefix("TableS") + "."):
                raise ValueError("Unexpected cycle contrast table")
            if list(next(rows)[:10]) != HEADER:
                raise ValueError("Unexpected cycle table header")
            for row in rows:
                symbol = row[8].strip() if isinstance(row[8], str) else ""
                if symbol in symbols:
                    genes.add(symbols[symbol])
    finally:
        book.close()
    return genes


def run(root):
    directory = root / "data/private/gse25628"
    days, receipt = audit(directory)
    if (
        receipt["ordered_sha256"]
        != "360626a37e13529bb37331a8f389a9ebe0e597ba1da43cbd1500e0843d439825"
    ):
        raise ValueError("Native scan metadata changed")
    days = np.array(days)
    fin, _, genes, _, xr, y, _, _ = study_inputs(root)
    weight = np.array([fin[g] for g in genes])
    score = ((xr - xr.mean(axis=0)) / xr.std(axis=0, ddof=1)) @ (weight / np.abs(weight).sum())
    endpoint = blocked_effect(score, y, days)
    cycling = cycle_genes(root / "data/private")
    flagged = np.isin(genes, list(cycling))
    phase = {}
    scaled = (xr - xr.mean(axis=0)) / xr.std(axis=0, ddof=1)
    for label, mask in [
        ("excluding_phase_table_genes", ~flagged),
        ("phase_table_genes_only", flagged),
    ]:
        w = weight[mask]
        result = blocked_effect(scaled[:, mask] @ (w / np.abs(w).sum()), y, days)
        result["genes"] = int(mask.sum())
        phase[label] = result
    assays = {name: blocked_effect(values, y, days) for name, values in controls(directory).items()}
    for item, adj in zip(
        assays.values(), holm([a["p_two_sided"] for a in assays.values()]), strict=True
    ):
        item["p_holm"] = float(adj)
    return dict(
        native_scan_metadata=receipt,
        Finnish_weighted_score=endpoint,
        technical_controls=assays,
        cycle_specificity=phase,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path("."))
    args = p.parse_args()
    print(json.dumps(run(args.repo), indent=2, allow_nan=False))
