"""Test external transfer of existing endometriosis transcriptome-wide contrasts.

Run in the Open-Disease repository root with PYTHONPATH=src:scripts.
Only aggregate outputs are returned. No patient predictions or fitted classifiers.
"""

from collections import defaultdict
import argparse
import csv
import gzip
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from bulk_gse153739 import load_mumbai
from genome_rome import checked_annotation, load_rome
from cross_cohort_discovery import holm
from verify_melbourne_s20 import HEADER, TITLE, WORKBOOK, sheet_rows

FINNISH_HASH = "9df65b29615ec56767e9abaaf3844f8ebde5e21baf8cdab193ac24857a170da5"
NOMINATION_HASH = "ba140f2bb026b0416e1c705ab28f616f430cc4295f94d174feaf5d15947c1e31"


def check_file(path, digest):
    if sha256(path.read_bytes()).hexdigest() != digest:
        raise ValueError(f"Input differs from the pinned snapshot: {path}")


def finnish(exposure):
    path = exposure / "discovery_gene.tsv"
    check_file(path, FINNISH_HASH)
    with path.open() as f:
        rows = {int(r["id"]): float(r["statistic"]) for r in csv.DictReader(f, delimiter="\t")}
    if len(rows) != 12044 or 2354 in rows or not np.isfinite(list(rows.values())).all():
        raise ValueError("Unexpected Finnish discovery universe")
    nominations = exposure / "nominations.json"
    check_file(nominations, NOMINATION_HASH)
    excluded = set()
    for item in json.loads(nominations.read_text())["candidates"]:
        if item["family"] == "gene":
            excluded.add(int(item["id"]))
        else:
            excluded.update(item["entrez_ids"])
    return rows, excluded


def unique_gene_symbols(directory):
    # Require bidirectional, unique symbol-to-gene and gene-to-symbol relationships
    # in the same audited probe annotation used by the Rome loader.
    checked = checked_annotation(directory)
    symbol_ids, id_symbols = defaultdict(set), defaultdict(set)
    with gzip.open(directory / "GPL571.annot.gz", "rt") as f:
        for line in f:
            if line.strip() == "!platform_table_begin":
                break
        for row in csv.DictReader(f, delimiter="\t"):
            probe = row["ID"]
            if probe == "!platform_table_end":
                break
            entrez = checked[probe]
            if entrez:
                symbol_ids[row["Gene symbol"]].add(entrez)
                id_symbols[entrez].add(row["Gene symbol"])
    unique = {}
    for symbol, ids in symbol_ids.items():
        if len(ids) == 1:
            (entrez,) = ids
            if id_symbols[entrez] == {symbol}:
                unique[symbol] = entrez
    return unique


def melbourne(workbook, directory):
    check_file(workbook, WORKBOOK)
    symbols = unique_gene_symbols(directory)
    rows = sheet_rows(workbook)
    if next(rows)[1]["A"] != TITLE:
        raise ValueError("Unexpected source table title")
    header = next(rows)[1]
    if [header[k] for k in "ABCDEFGHIJ"] != HEADER:
        raise ValueError("Unexpected source table header")
    seen, by_gene = set(), defaultdict(list)
    for _, row in rows:
        gene = row.get("I", "").strip()
        if row["A"] in seen:
            raise ValueError("Repeated Melbourne probe")
        seen.add(row["A"])
        if gene in symbols:
            t = float(row["D"])
            if not np.isfinite(t):
                raise ValueError("Invalid source statistic")
            by_gene[symbols[gene]].append(t)
    if len(seen) != 15262:
        raise ValueError("Incomplete Melbourne table")
    return {gene: ts[0] for gene, ts in by_gene.items() if len(ts) == 1}


def contrast(scores, case):
    return scores[case].mean(axis=0) - scores[~case].mean(axis=0)


def evaluate_projection(scores, case, weights):
    if (
        scores.shape != (len(case), len(weights))
        or case.dtype != np.bool_
        or min(case.sum(), (~case).sum()) < 2
    ):
        raise ValueError("Invalid score or donor design")
    if np.any(np.ptp(scores, axis=0) == 0) or not np.any(weights):
        raise ValueError("Constant gene or all-zero weight vector")
    if not np.isfinite(scores).all() or not np.isfinite(weights).all():
        raise ValueError("Nonfinite input")
    weights = weights / np.abs(weights).sum()
    scaled = (scores - scores.mean(axis=0)) / scores.std(axis=0, ddof=1)
    projected = scaled @ weights
    diff = float(contrast(projected, case))
    counts = 0
    alternative_counts = 0
    omissions = []
    # Equivalent calculation: sum pre-calculated within-donor gene projections
    # with a pooled-total identity and Python fsum rather than two group means.
    gene_projection = np.sum(scaled * weights[None, :], axis=1)
    exact_total = math.fsum(gene_projection)
    n, cases = len(case), int(case.sum())
    for chosen in combinations(range(n), cases):
        selected = np.zeros(n, dtype=bool)
        selected[list(chosen)] = True
        delta = float(contrast(projected, selected))
        sump = math.fsum(gene_projection[list(chosen)])
        alternative = sump / cases - (exact_total - sump) / (n - cases)
        if abs(delta - alternative) > 2e-13:
            raise ValueError("Permutation calculations disagree")
        counts += abs(delta) >= abs(diff) - 1e-12
        alternative_counts += abs(alternative) >= abs(diff) - 1e-12
    if counts != alternative_counts or not np.isclose(
        diff, contrast(scaled, case) @ weights, atol=1e-13
    ):
        raise ValueError("Projection and direct contrast disagree")
    for i in range(n):
        retained = np.arange(n) != i
        omissions.append(contrast(projected[retained], case[retained]))
    external = contrast(scaled, case)
    informative = (weights != 0) & (np.abs(external) > 1e-12)
    denom = math.comb(n, cases)
    return dict(
        genes=len(weights),
        cases=cases,
        controls=n - cases,
        score_difference=diff,
        extreme=counts,
        assignments=denom,
        p_two_sided=counts / denom,
        omit_one_range=[min(omissions), max(omissions)],
        unscaled_difference_pp=float(contrast(scores, case) @ weights),
        spearman=float(spearmanr(weights, external).statistic),
        genes_nonzero=int(informative.sum()),
        same_sign=int(np.sum(np.sign(weights[informative]) == np.sign(external[informative]))),
    )


def study_inputs(root):
    private = root / "data/private"
    exposure = root / "evidence/cross-cohort-2026-11-14"
    fin, removed = finnish(exposure)
    rome, xr, yr, _ = load_rome(private / "gse25628")
    mumbai, xm, ym, _ = load_mumbai(
        private / "gse153739", private / "gse266265/hgnc_complete_set_2026-10-06.tsv"
    )
    genes = np.array(sorted(fin.keys() & set(rome) & set(mumbai)))
    xr, xm = xr[:, np.searchsorted(rome, genes)], xm[:, np.searchsorted(mumbai, genes)]
    eligible = (np.ptp(xr, axis=0) > 0) & (np.ptp(xm, axis=0) > 0)
    initially = len(genes)
    genes, xr, xm = genes[eligible], xr[:, eligible], xm[:, eligible]
    return fin, removed, genes, initially, xr, yr, xm, ym


def run(root):
    fin, removed, genes, initially, xr, yr, xm, ym = study_inputs(root)
    private = root / "data/private"
    mweights = melbourne(private / "fung-2018/Fung_supplementary_tables.xlsx", private / "gse25628")
    checks = {}
    for name, source, excluded in [
        ("Finland", fin, set()),
        ("Finland_without_nominees", fin, removed),
        ("Melbourne_unoriented", mweights, set()),
    ]:
        keep = np.array([g in source and g not in excluded for g in genes])
        w = np.array([source[g] for g in genes[keep]])
        checks[name] = {
            study: evaluate_projection(x[:, keep], y, w)
            for study, x, y in [("Rome", xr, yr), ("Mumbai", xm, ym)]
        }
        if name != "Finland_without_nominees":
            results = checks[name].values()
            for result, adjusted in zip(
                results, holm([r["p_two_sided"] for r in results]), strict=True
            ):
                result["p_holm"] = float(adjusted)
    shared = [g for g in genes if g in mweights]
    return dict(
        initially_shared=initially,
        variable_shared=len(genes),
        melbourne_mapped_single=len(mweights),
        finland_melbourne_spearman_unoriented=float(
            spearmanr([fin[g] for g in shared], [mweights[g] for g in shared]).statistic
        ),
        families=checks,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    args = parser.parse_args()
    print(json.dumps(run(args.repo.resolve()), indent=2, allow_nan=False))
