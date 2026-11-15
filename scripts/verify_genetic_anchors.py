"""Re-enumerate primary Rome and Mumbai genetic-panel tests by cohort mean differences."""
from itertools import combinations, product
import argparse
import json
from pathlib import Path

import numpy as np
from bulk_gse153739 import load_mumbai
from evaluate_genetic_anchors import read_table
from genome_rome import load_rome
from rome_scan_dates import uterine_scan_dates


def verify_design(genes, x, case, dates, scores, source, reference):
    groups = [np.flatnonzero(dates == d) for d in np.unique(dates)]
    groups = [idx for idx in groups if case[idx].any() and not case[idx].all()]
    x = x[:,np.searchsorted(genes,scores)]
    options = []
    for idx in groups:
        count = int(case[idx].sum())
        options.append([x[list(sub)].mean(axis=0)-x[list(set(idx)-set(sub))].mean(axis=0)
                        for sub in combinations(idx,count)])
    delta = np.array([np.mean(p,axis=0) for p in product(*options)])
    obs = np.array([x[idx[case[idx]]].mean(axis=0)-x[idx[~case[idx]]].mean(axis=0)
                    for idx in groups]).mean(axis=0)
    z = delta/delta.std(axis=0)
    zobs = obs/delta.std(axis=0)
    if source == "smr":
        null = np.square(z).mean(axis=1)
        observed = np.square(zobs).mean()
    else:
        signs = np.array([source[g] for g in scores])
        null = (z*signs).mean(axis=1)
        observed = (zobs*signs).mean()
    total = len(delta)
    extreme = int((null >= observed-1e-10).sum())
    # A global fixed change of scale cancels in the signed-score comparisons.
    if extreme != reference["extreme"] or abs(extreme/total-reference["p"]) > 1e-12:
        raise ValueError("Independent enumeration disagrees")
    return extreme,total


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("rome","scans","mumbai","hgnc","panel"):
        parser.add_argument("--"+name,type=Path,required=True)
    a=parser.parse_args()
    scan=json.loads((a.panel/"scan_day.json").read_text())
    transfer=json.loads((a.panel/"transfer.json").read_text())
    primary=json.loads((a.panel/"omnibus.json").read_text())
    smr=[int(r["gene_id"]) for r in read_table(a.panel/"rahmioglu_expression.tsv") if r["measured_both"]=="True"]
    twas={int(r["gene_id"]):np.sign(float(r["twas_z"])) for r in read_table(a.panel/"koller_uterus.tsv")}
    genes,x,case,_=load_rome(a.rome)
    days=np.array(uterine_scan_dates(a.scans,a.rome/"GSE25628_series_matrix.txt.gz"))
    unblocked=np.zeros(len(case))
    for days_ref,ref1,ref2 in [(unblocked,transfer["rome_smr"]["self_contained"],
                                dict(extreme=primary["uterus_twas"]["extreme_one_sided"],p=primary["uterus_twas"]["p"])),
                               (days,scan["smr"],scan["rbm6"])]:
        print(verify_design(genes,x,case,days_ref,smr,"smr",ref1))
        print(verify_design(genes,x,case,days_ref,[10180],twas,ref2))
    genes,x,case,_=load_mumbai(a.mumbai,a.hgnc)
    days=np.zeros(len(case))
    print(verify_design(genes,x,case,days,np.intersect1d(smr,genes),"smr",transfer["mumbai_smr"]["self_contained"]))
    print(verify_design(genes,x,case,days,list(twas),twas,transfer["mumbai_twas"]))
