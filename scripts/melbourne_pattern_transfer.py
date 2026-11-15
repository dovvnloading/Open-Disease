"""Exact Rome donor-label tests of published independent Melbourne expression patterns."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from openpyxl import load_workbook

from analyze_melbourne_genes import SOURCE_HASH
from cross_cohort_discovery import assignment_design, holm
from genome_rome import checked_annotation, load_rome
from open_disease.geo_microarray import sha256_file

FINNISH_TABLE_HASH = "9df65b29615ec56767e9abaaf3844f8ebde5e21baf8cdab193ac24857a170da5"
SHEETS = {"TableS20": 15262, "TableS3": 1555}


def symbol_bridge(directory):
    """Use only one-to-one symbol/Entrez pairs in the valid GPL571 features."""
    mapped = checked_annotation(directory)
    by_symbol = defaultdict(set)
    by_id = defaultdict(set)
    with gzip.open(directory / "GPL571.annot.gz", "rt", newline="") as f:
        for line in f:
            if line.rstrip() == "!platform_table_begin":
                break
        for row in csv.DictReader(f, delimiter="\t"):
            if row["ID"] == "!platform_table_end":
                break
            gene = mapped[row["ID"]]
            if gene:
                symbol = row["Gene symbol"]
                by_symbol[symbol].add(gene)
                by_id[gene].add(symbol)
    return {symbol: next(iter(ids)) for symbol, ids in by_symbol.items()
            if len(ids) == 1 and len(by_id[next(iter(ids))]) == 1}


def signature(sheet, bridge, shared):
    rows = iter(sheet.values)
    title, header = next(rows), next(rows)
    if (not title[0].startswith(f"Table {sheet.title[5:]}.") or
            tuple(header[:5]) != ("PROBE_ID", "logFC", "AveExpr", "t", "P.Value")
            or header[8] != "ILMN_GENE"):
        raise ValueError("Unexpected Melbourne table layout")
    observed = set()
    effects = defaultdict(list)
    for row in rows:
        if row[0] in observed or not str(row[0]).startswith("ILMN_") or not np.isfinite(row[3]):
            raise ValueError("Invalid Melbourne probe identifier or statistic")
        observed.add(row[0])
        if sheet.title == "TableS3":
            if (row[10] not in {"Up in ES", "Down in ES"} or
                    (row[3] > 0) != (row[10] == "Up in ES")):
                raise ValueError("Published phase-direction annotation disagrees with statistic")
        symbol = row[8].strip() if isinstance(row[8], str) else ""
        gene = bridge.get(symbol)
        if gene in shared:
            effects[gene].append(float(row[3]))
    if len(observed) != SHEETS[sheet.title] or not effects:
        raise ValueError("Melbourne probe universe changed or no genes mapped")
    ids = sorted(effects)
    return ids, np.array([np.mean(effects[gene]) for gene in ids]), sum(map(len, effects.values()))


def standardize(scores):
    scale = scores.std(axis=0, ddof=1)
    if np.any(scale <= 0):
        raise ValueError("A reference gene has zero variation")
    return (scores - scores.mean(axis=0)) / scale


def exact_correlations(z, signature_t, weights):
    """Compute Pearson r across genes for each whole-donor assignment via a Gram matrix."""
    centered = z - z.mean(axis=1, keepdims=True)
    ref = signature_t - signature_t.mean()
    norm = np.linalg.norm(ref)
    if not norm:
        raise ValueError("Constant external reference")
    numerator = weights @ (centered @ ref)
    denominator = np.sqrt(np.einsum("ij,ij->i", weights @ (centered @ centered.T), weights))
    if np.any(denominator <= 0):
        raise ValueError("Constant permuted gene contrast")
    return numerator / (denominator * norm)


def calculate(scores, case, signature_t, fin_stat):
    observed, weights, _ = assignment_design(case, np.zeros(len(case), dtype=int))
    perm = exact_correlations(standardize(scores), signature_t, weights)
    index, = np.where(np.all(np.isclose(weights, observed, rtol=0, atol=1e-14), axis=1))[0]
    rho = float(perm[index])
    rome_t = observed @ standardize(scores)
    direct_r = np.corrcoef(rome_t, signature_t)[0, 1]
    if not np.isclose(rho, direct_r, atol=1e-12, rtol=0):
        raise ValueError("Direct and projected Pearson correlations disagree")
    loo = []
    for i in range(len(case)):
        keep = np.arange(len(case)) != i
        labels = case[keep]
        z = standardize(scores[keep])
        contrast = z[labels].mean(axis=0) - z[~labels].mean(axis=0)
        loo.append(float(np.corrcoef(contrast, signature_t)[0, 1]))
    extreme = int((np.abs(perm) >= abs(rho) - 1e-12).sum())
    return dict(rome_r=rho, rome_extreme=extreme, assignments=len(perm), p_two_sided=extreme / len(perm),
                null_central_95=np.quantile(perm, [0.025, 0.975]).tolist(),
                rome_loo_r_range=[min(loo), max(loo)],
                finland_r=float(np.corrcoef(fin_stat, signature_t)[0, 1]),
                rome_abs_magnitude_r=float(np.corrcoef(np.abs(rome_t), np.abs(signature_t))[0, 1]),
                finland_abs_magnitude_r=float(np.corrcoef(np.abs(fin_stat), np.abs(signature_t))[0, 1]))


def scan_day_sensitivity(scores, case, dates, ref):
    """Keep the 14-array gene scale fixed and compare groups within mixed dates."""
    mixed = [d for d in np.unique(dates) if len(set(case[dates == d])) == 2]
    keep = np.isin(dates, mixed)
    c = case[keep]
    z = standardize(scores)[keep]
    observed, weights, _ = assignment_design(c, dates[keep])
    null = exact_correlations(z, ref, weights)
    r = float(np.corrcoef(observed @ z, ref)[0, 1])
    direct_null = np.array([np.corrcoef(row, ref)[0, 1] for row in weights @ z])
    if not np.allclose(null, direct_null, rtol=0, atol=1e-11):
        raise ValueError("Direct and Gram-matrix date-conditioned correlations disagree")
    extreme = int((np.abs(null) >= abs(r) - 1e-12).sum())
    return dict(cases=int(c.sum()), controls=int((~c).sum()), n_dates=len(mixed),
                excluded=int((~keep).sum()), r=r, extreme=extreme, assignments=len(null),
                p_two_sided=extreme / len(null))


def run(directory, workbook, finnish, scan_audit=None):
    if sha256_file(workbook) != SOURCE_HASH or sha256_file(finnish) != FINNISH_TABLE_HASH:
        raise ValueError("Unexpected Melbourne workbook or Finnish screen")
    with finnish.open(newline="") as f:
        fin = {int(row["id"]): float(row["statistic"]) for row in csv.DictReader(f, delimiter="\t")}
    if len(fin) != 12044:
        raise ValueError("Unexpected number of Finnish genes")
    genes, scores, case, _ = load_rome(directory)
    if len(case) != 14 or sum(case) != 8:
        raise ValueError("Expected 8 Rome cases and 6 controls")
    dates = None
    if scan_audit is not None:
        from rome_scan_dates import scan_records
        records = scan_records(scan_audit, directory / "GSE25628_series_matrix.txt.gz")
        uterine = [(group, stamp) for group, stamp, _ in records if group != "lesion"]
        if not np.array_equal(case, [group == "case" for group, _ in uterine]):
            raise ValueError("Rome native scan and expression group order disagree")
        dates = np.array([stamp.date() for _, stamp in uterine])
    bridge = symbol_bridge(directory)
    book = load_workbook(workbook, read_only=True, data_only=True)
    result = []
    try:
        for name, n_source in SHEETS.items():
            ids, reference, probes = signature(book[name], bridge, set(fin))
            cols = np.searchsorted(genes, ids)
            if not np.array_equal(genes[cols], ids):
                raise ValueError("Rome is missing a reference gene")
            row = calculate(scores[:, cols], case, reference, np.array([fin[i] for i in ids]))
            if dates is not None:
                row["mixed_scan_dates"] = scan_day_sensitivity(scores[:, cols], case, dates, reference)
            result.append(dict(sheet=name, published_probes=n_source, mapped_probes=probes, genes=len(ids), **row))
    finally:
        book.close()
    for row, adjusted in zip(result, holm([r["p_two_sided"] for r in result]), strict=True):
        row["p_holm"] = float(adjusted)
    if dates is not None:
        blocks = [r["mixed_scan_dates"] for r in result]
        for block, p in zip(blocks, holm([b["p_two_sided"] for b in blocks]), strict=True):
            block["p_holm"] = float(p)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rome", type=Path, required=True)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--finnish", type=Path, required=True)
    parser.add_argument("--scan-audit", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(args.rome, args.workbook, args.finnish, args.scan_audit), indent=2))
