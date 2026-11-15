"""Describe fixed uterine TWAS genes using the audited eight-biopsy RNA cohort."""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from analyze_gse315857_fosb import MATRIX_HASH, N_NAMED, checked_cohort
from bulk_gse153739 import ATTRS, DIGESTS, mapping
from evaluate_genetic_anchors import read_table
from open_disease.geo_microarray import sha256_file
from swiss_epithelial_mapping import HGNC_SHA256


def describe(directory, gencode, hgnc, panel):
    matrix = directory / "counts.download"
    for file,pin in [(matrix,MATRIX_HASH),(gencode/ATTRS,DIGESTS[ATTRS]),(hgnc,HGNC_SHA256)]:
        if sha256_file(file) != pin:
            raise ValueError("Unverified input")
    groups = checked_cohort(directory / "childgpl.download",directory / "ena.download")
    with gzip.open(matrix,"rt",encoding="utf-8") as stream:
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        if header[0] != "EnsembleID" or len(set(header[1:])) != 8 or set(header[1:]) != groups.keys():
            raise ValueError("Group labels do not align")
        case = np.array([groups[label] == "Eosis" for label in header[1:]])
        rows = [row for row in reader if row[0]]
    tx = [r[0] for r in rows]
    if len(tx) != N_NAMED or len(set(tx)) != N_NAMED:
        raise ValueError("Unexpected transcript universe")
    values = np.asarray([r[1:] for r in rows],dtype=float)
    pr = 100*(rankdata(values,method="average",axis=0)-1)/(N_NAMED-1)
    ids,_ = mapping(gencode,hgnc,tx)
    panel_rows = read_table(panel/"koller_uterus.tsv")
    signed, per_gene = [], {}
    for r in panel_rows:
        mask = ids == int(r["gene_id"])
        stats = {"transcripts":int(mask.sum())}
        per_gene[r["symbol"]] = stats
        if not mask.any():
            continue
        v = pr[mask].mean(axis=0)
        if v.std(ddof=1) == 0:
            raise ValueError("Constant genetic-anchor measure")
        stats["effect_pp"] = float(v[case].mean()-v[~case].mean())
        signed.append(np.sign(float(r["twas_z"])) * (v-v.mean())/v.std(ddof=1))
    if not signed:
        raise ValueError("No assayed genes")
    score = np.mean(signed,axis=0)
    leave_one = []
    for i in range(8):
        keep = np.arange(8) != i
        leave_one.append(score[case & keep].mean()-score[~case & keep].mean())
    return {"cases":int(case.sum()),"controls":int((~case).sum()),"gene_effects":per_gene,
            "n_genes":len(signed), "score_difference":float(score[case].mean()-score[~case].mean()),
            "leave_one_donor_fixed_scale":[float(min(leave_one)),float(max(leave_one))]}


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ("directory","gencode","hgnc","panel","out"):
        parser.add_argument("--"+name,type=Path,required=True)
    args=parser.parse_args()
    result=describe(args.directory,args.gencode,args.hgnc,args.panel)
    args.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2))
