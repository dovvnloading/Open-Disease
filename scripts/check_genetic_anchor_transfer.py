"""Check genetic-panel specificity in Rome and transfer fixed panels into Mumbai."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from bulk_gse153739 import load_mumbai
from cross_cohort_discovery import TOL, assignment_design, holm
from evaluate_genetic_anchors import FINLAND, read_table
from genome_rome import load_rome
from open_disease.geo_microarray import sha256_file


def standardized(matrix):
    centered = matrix - matrix.mean(axis=0)
    sd = centered.std(axis=0, ddof=1)
    variable = sd > 1e-12
    return centered[:, variable] / sd[variable], variable


def upper_tail(value, reference):
    extreme = int(np.count_nonzero(reference >= value - TOL))
    return {"value":float(value), "extreme": extreme, "assignments":len(reference),
            "p":extreme / len(reference)}


def panel_check(genes, matrix, labels, common, panel):
    measured = np.intersect1d(genes, common)
    x, variable = standardized(matrix[:,np.searchsorted(genes, measured)])
    included = measured[variable]
    mask = np.isin(included, panel)
    obs, w, _ = assignment_design(labels, np.zeros(len(labels), dtype=int))
    # Every assignment has the same case-control counts and null variance.
    scale = np.sqrt(1 / labels.sum() + 1 / (~labels).sum())
    square = np.square((obs @ x) / scale)
    in_null, out_null = np.zeros(len(w)), np.zeros(len(w))
    for start in range(0, len(included), 128):
        sl = slice(start, start+128)
        z2 = np.square((w @ x[:,sl]) / scale)
        in_null += z2[:,mask[sl]].sum(axis=1)
        out_null += z2[:,~mask[sl]].sum(axis=1)
    in_null /= mask.sum()
    out_null /= (~mask).sum()
    observed = square[mask].mean()
    result = {
        "panel_genes":int(mask.sum()), "background_genes":int((~mask).sum()),
        "missing_panel_ids":sorted(set(map(int,panel))-set(map(int,included[mask]))),
        "inactive_common_genes":int((~variable).sum()), "self_contained":upper_tail(observed, in_null),
        "background_mean_squared_z":float(square[~mask].mean()),
        "competitive":upper_tail(observed-square[~mask].mean(), in_null-out_null),
    }
    rbm = np.flatnonzero(included == 10180)
    if len(rbm):
        result["rbm6_abs_z_rank"] = 1 + int((square > square[rbm[0]]+TOL).sum())
    return result


def signed_check(genes, matrix, labels, twas):
    ids = np.array([int(r["gene_id"]) for r in twas])
    if not np.isin(ids, genes).all():
        raise ValueError("Missing Mumbai TWAS gene")
    x = matrix[:,np.searchsorted(genes, ids)]
    scaled, variable = standardized(x)
    if not variable.all():
        raise ValueError("Constant TWAS component")
    aligned = scaled * np.sign([float(r["twas_z"]) for r in twas])
    obs, w, _ = assignment_design(labels, np.zeros(len(labels), dtype=int))
    score = aligned.mean(axis=1)
    result = upper_tail(obs @ score, w @ score)
    for mask, key in [(ids == 10180, "rbm6_only"), (ids != 10180, "excluding_rbm6")]:
        sensitivity = aligned[:,mask].mean(axis=1)
        result[key] = upper_tail(obs @ sensitivity, w @ sensitivity)
    result["per_gene_effects_pp"] = {r["symbol"]:float(v) for r,v in zip(twas, obs @ x, strict=True)}
    return result


def main(a):
    if sha256_file(a.finland) != FINLAND:
        raise ValueError("Unexpected Finnish screen checksum")
    common = np.array([int(r["id"]) for r in read_table(a.finland)])
    registry = read_table(a.panel / "rahmioglu_table1.tsv")
    panel = np.intersect1d(common,[int(r["gene_id"]) for r in registry])
    twas = read_table(a.panel / "koller_uterus.tsv")
    if len(panel) != 15 or {r["symbol"] for r in twas} != {"GMPPB","CRISPLD1","RBM6","CISD2"}:
        raise ValueError("Expected frozen source panels")
    genes, x, c, _ = load_rome(a.rome)
    results = {"rome_smr":panel_check(genes, x, c, common, panel)}
    genes, x, c, alt = load_mumbai(a.mumbai, a.hgnc)
    if (len(c),c.sum()) != (7,4):
        raise ValueError("Unexpected Mumbai design")
    results["mumbai_smr"] = panel_check(genes, x, c, common, panel)
    results["mumbai_twas"] = signed_check(genes, x, c, twas)
    results["mumbai_twas_alternative"] = signed_check(genes, alt, c, twas)
    primary = [results["mumbai_smr"]["self_contained"],results["mumbai_twas"]]
    for item,p in zip(primary,holm([r["p"] for r in primary]),strict=True):
        item["holm_2"] = float(p)
    a.out.write_text(json.dumps(results,indent=2,allow_nan=False)+"\n")
    print(json.dumps(results,indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["finland", "rome", "mumbai", "hgnc", "panel", "out"]:
        parser.add_argument("--"+name,type=Path,required=True)
    main(parser.parse_args())
