"""Repeat Rome genetic-panel tests within native scan dates shared by both diagnoses."""
import argparse
import json
from pathlib import Path

import numpy as np
from check_genetic_anchor_transfer import upper_tail
from cross_cohort_discovery import assignment_design, holm
from evaluate_genetic_anchors import FINLAND, read_table
from genome_rome import load_rome
from open_disease.geo_microarray import sha256_file
from rome_scan_dates import uterine_scan_dates


def analyze(rome, scans, finland, panel):
    if sha256_file(finland) != FINLAND:
        raise ValueError("Unexpected Finnish summary")
    genes, x, c, _ = load_rome(rome)
    dates = np.array(uterine_scan_dates(scans, rome / "GSE25628_series_matrix.txt.gz"))
    mixed = [day for day in np.unique(dates) if len(set(c[dates == day]))==2]
    keep = np.isin(dates, mixed)
    obs, w, _ = assignment_design(c[keep],dates[keep])
    if (len(w),len(mixed),sum(c[keep]),sum(~c[keep])) != (200,2,5,6):
        raise ValueError("Unexpected shared-day design")
    common = np.array([int(r["id"]) for r in read_table(finland)])
    if not np.isin(common, genes).all():
        raise ValueError("Changed Rome background")
    x = x[keep][:,np.searchsorted(genes,common)]
    # Subtract per-date means before permutation to eliminate constant offsets.
    for day in mixed:
        group = dates[keep] == day
        x[group] -= x[group].mean(axis=0)
    null = w @ x
    sd = null.std(axis=0)
    if not (sd > 1e-12).all():
        raise ValueError("Constant gene in the shared universe")
    z, null_z = (obs @ x)/sd, null/sd
    ref = read_table(panel/"rahmioglu_table1.tsv")
    symbols = {int(r["gene_id"]):r["paper_symbol"] for r in ref}
    panel_mask = np.isin(common,list(symbols))
    if panel_mask.sum() != 15:
        raise ValueError("Unexpected conditional genetic panel")
    z2, null2 = z**2, null_z**2
    primary = upper_tail(z2[panel_mask].mean(),null2[:,panel_mask].mean(axis=1))
    primary["background_mean_squared_z"] = float(z2[~panel_mask].mean())
    competitive = upper_tail(z2[panel_mask].mean()-z2[~panel_mask].mean(),
                             null2[:,panel_mask].mean(axis=1)-null2[:,~panel_mask].mean(axis=1))
    k = list(common).index(10180)
    rbm6 = upper_tail(-z[k],-null_z[:,k])
    rbm6["effect_pp"] = float(obs @ x[:,k])
    for r,p in zip([primary,rbm6],holm([primary["p"],rbm6["p"]]),strict=True):
        r["holm_2"] = float(p)
    individual = {}
    effects = obs @ x
    for i in np.flatnonzero(panel_mask):
        individual[symbols[common[i]]] = dict(effect_pp=float(effects[i]), **upper_tail(abs(z[i]),abs(null_z[:,i])))
    return {"shared_dates":len(mixed), "cases":int(sum(c[keep])), "controls":int(sum(~c[keep])),
            "smr":primary,"smr_competitive":competitive,"rbm6":rbm6,
            "smr_genes_descriptive":individual}


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ["rome","scans","finland","panel","out"]:
        parser.add_argument("--"+name,type=Path,required=True)
    a=parser.parse_args()
    result=analyze(a.rome,a.scans,a.finland,a.panel)
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="smr_genes_descriptive"},indent=2))
