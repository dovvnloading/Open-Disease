"""Compare published endometriosis genetic anchors with case-control uterine ranks."""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np
from openpyxl import load_workbook
from scipy.stats import rankdata

from cross_cohort_discovery import TOL, assignment_design, holm, write_tsv
from genome_rome import load_rome
from open_disease.geo_microarray import sha256_file
from verify_melbourne_s20 import HEADER, TITLE, WORKBOOK, sheet_rows

FINLAND = "9df65b29615ec56767e9abaaf3844f8ebde5e21baf8cdab193ac24857a170da5"
KOLLER = "c021b6c20c668644939fd83381307a2c8cf08f18d95156b33e65fb0301997815"
THRESHOLD = 3.72e-6


def read_table(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def twas_panel(book_path, crosswalk):
    mapping = {r["symbol"]: r["gene_id"] for r in read_table(crosswalk)}
    book = load_workbook(book_path, data_only=True, read_only=True)
    selected = []
    n_total = n_uterus = 0
    try:
        rows = iter(book["Supplementary Table 13"].values)
        if "EUR combined definition" not in next(rows)[0]:
            raise ValueError("Unexpected genetic phenotype")
        next(rows)
        header = next(rows)
        if tuple(header[6:8]) != ("Uterus z-score", "Uterus p-value"):
            raise ValueError("Unexpected uterine TWAS columns")
        for row in rows:
            transcript, symbol = row[:2]
            z, p = row[6:8]
            n_total += 1
            if (z, p) == ("NA", "NA"):
                continue
            if not np.isfinite(z) or not 0 <= p <= 1:
                raise ValueError("Invalid uterine TWAS result")
            n_uterus += 1
            if p < THRESHOLD:
                selected.append({"symbol": symbol, "gene_id": mapping[symbol],
                                 "ensembl": transcript, "twas_z": z, "twas_p": p})
    finally:
        book.close()
    if (n_total, n_uterus) != (13457, 2533) or set(mapping) != {r["symbol"] for r in selected}:
        raise ValueError("TWAS coverage or eligible set changed")
    return selected


def melbourne_rows(path, registry):
    by_name = {symbol: r["gene_id"] for r in registry for symbol in (r["paper_symbol"], r["current_symbol"])}
    probes = defaultdict(list)
    rows = sheet_rows(path)
    if next(rows)[1]["A"] != TITLE or [value for value in next(rows)[1].values()][:10] != HEADER:
        raise ValueError("Unexpected Melbourne Table S20")
    seen = set()
    for _, row in rows:
        if row["A"] in seen:
            raise ValueError("Duplicate Melbourne probe")
        seen.add(row["A"])
        gene = by_name.get(row["I"].strip())
        if gene is not None:
            probes[gene].append((row["A"], float(row["E"])))
    if len(seen) != 15262:
        raise ValueError("Unexpected Melbourne Table S20 size")
    return probes


def mean_square_test(z_null, z_obs):
    if not z_obs.size:
        raise ValueError("Empty panel")
    observed = np.square(z_obs).mean()
    null = np.square(z_null).mean(axis=1)
    extreme = int((null >= observed - TOL).sum())
    return {"genes": len(z_obs), "mean_squared_z": float(observed),
            "extreme": extreme, "p": extreme / len(null)}


def main(args):
    for path, pin in [(args.finland, FINLAND), (args.koller, KOLLER), (args.melbourne, WORKBOOK)]:
        if sha256_file(path) != pin:
            raise ValueError(f"Source checksum failed for {path.name}")
    args.out.mkdir(parents=True, exist_ok=True)
    registry = read_table(args.definitions / "rahmioglu_table1.tsv")
    if len(registry) != 23 or len({r["gene_id"] for r in registry}) != 23:
        raise ValueError("Incomplete or duplicated genetic panel")
    twas = twas_panel(args.koller, args.definitions / "koller_crosswalk.tsv")
    mel = melbourne_rows(args.melbourne, registry + [
        {"gene_id": r["gene_id"], "paper_symbol": r["symbol"], "current_symbol": r["symbol"]} for r in twas
    ])
    finnish_rows = read_table(args.finland)
    fin = {r["id"]: r for r in finnish_rows}
    f_ranks = rankdata([-abs(float(r["statistic"])) for r in finnish_rows], method="min")
    ranks = {r["id"]: int(rank) for r, rank in zip(finnish_rows, f_ranks, strict=True)}
    genes, scores, case, strata = load_rome(args.rome)
    obs, weights, _ = assignment_design(case, strata)
    if (len(weights), case.sum(), len(case), len(fin)) != (3003, 8, 14, 12044):
        raise ValueError("Unexpected cohort design")
    shared = (set(fin) & {str(g) for g in genes})
    if shared != set(fin):
        raise ValueError("Incomplete Rome coverage of the frozen common universe")
    measured = [r for r in registry if r["gene_id"] in shared]
    ids = np.array([int(r["gene_id"]) for r in measured])
    x = scores[:, np.searchsorted(genes, ids)]
    x = x - x.mean(axis=0)
    sd_null = np.sqrt(1 / case.sum() + 1 / (~case).sum()) * x.std(axis=0, ddof=1)
    if np.any(sd_null <= 0):
        raise ValueError("Uninformative gene")
    z_obs, z_null = (obs @ x) / sd_null, (weights @ x) / sd_null
    effect = obs @ x
    extreme = (np.abs(z_null) >= np.abs(z_obs) - TOL).sum(axis=0)
    index = {r["gene_id"]: i for i, r in enumerate(measured)}
    results = []
    test_inputs = []
    for r in registry:
        g = r["gene_id"]
        probes = mel[g]
        entry = {"symbol":r["paper_symbol"], "gene_id":g, "measured_both": g in index,
                 "finland_effect_pp":"", "finland_rank_abs_z":"", "finland_p":"",
                 "rome_effect_pp":"", "rome_extreme_two_sided":"", "rome_p":"",
                 "same_sign":"", "conjunction_p":"", "holm_23":"",
                 "melbourne_probes":len(probes),
                 "melbourne_within_gene_p":min(1., len(probes)*min(p for _,p in probes)) if probes else ""}
        if g in index:
            i = index[g]
            f_effect, f_p = float(fin[g]["effect_pp"]), float(fin[g]["p_raw"])
            r_p = extreme[i] / len(weights)
            agree = bool(f_effect * effect[i] > 0)
            conj = max(f_p, r_p) if agree else 1.
            entry.update(finland_effect_pp=f_effect, finland_rank_abs_z=ranks[g], finland_p=f_p,
                         rome_effect_pp=effect[i], rome_extreme_two_sided=int(extreme[i]),
                         rome_p=r_p, same_sign=agree, conjunction_p=conj)
        results.append(entry)
        test_inputs.append(entry["conjunction_p"] if g in index else 1.)
    for r, p in zip(results, holm(test_inputs), strict=True):
        if r["measured_both"]:
            r["holm_23"] = p
    write_tsv(args.out/"rahmioglu_expression.tsv", list(results[0]), results)
    summary = {"rome_permutations": len(weights), "table1_all": mean_square_test(z_null, z_obs)}
    for key, mask in [
        ("table1_omit_abo", [r["paper_symbol"] != "ABO" for r in measured]),
        ("table1_endometrium_eqtl", ["endometrium_eQTL" in r["molecular_evidence"].split(",") for r in measured]),
    ]:
        summary[key] = mean_square_test(z_null[:,mask], z_obs[mask])
    t_measured = [r for r in twas if r["gene_id"] in shared]
    t_ids = [int(r["gene_id"]) for r in t_measured]
    t_mat = scores[:, np.searchsorted(genes, t_ids)]
    sd = t_mat.std(axis=0, ddof=1)
    if not t_ids or np.any(sd <= 0):
        raise ValueError("Inadequate TWAS panel measurements")
    aligned = ((t_mat - t_mat.mean(axis=0)) / sd) * np.sign([r["twas_z"] for r in t_measured])
    score = aligned.mean(axis=1)
    observed, null = float(obs @ score), weights @ score
    t_extreme = int((null >= observed - TOL).sum())
    loo = []
    for i in range(len(case)):
        keep = np.arange(len(case)) != i
        loo.append(score[case & keep].mean() - score[~case & keep].mean())
    gene_loo = [float(obs @ np.delete(aligned, i, axis=1).mean(axis=1)) for i in range(len(t_ids)) if len(t_ids) > 1]
    summary["uterus_twas"] = {
        "source_models":2533, "eligible":len(twas), "measured":len(t_ids),
        "aligned_mean_sd_difference":observed, "extreme_one_sided":t_extreme,
        "p": t_extreme / len(weights),
        "range_omit_donor_fixed_scale":[float(min(loo)), float(max(loo))],
        "range_omit_gene":[min(gene_loo), max(gene_loo)] if gene_loo else None,
        "finland_mean_aligned_z":float(np.mean([
            np.sign(r["twas_z"])*float(fin[r["gene_id"]]["statistic"]) for r in t_measured])),
    }
    for r, p in zip([summary["table1_all"], summary["uterus_twas"]],
                    holm([summary["table1_all"]["p"], summary["uterus_twas"]["p"]]), strict=True):
        r["holm_2"] = float(p)
    t_index = {r["gene_id"]: i for i, r in enumerate(t_measured)}
    for r in twas:
        g = r["gene_id"]
        r["measured_both"] = g in t_index
        r["finland_effect_pp"] = fin[g]["effect_pp"] if g in t_index else ""
        r["rome_effect_pp"] = float(obs @ t_mat[:,t_index[g]]) if g in t_index else ""
        probes = mel[g]
        r["melbourne_probes"] = len(probes)
        r["melbourne_within_gene_p"] = min(1.,len(probes)*min(p for _,p in probes)) if probes else ""
    mel_inputs = [r["melbourne_within_gene_p"] if r["melbourne_probes"] else 1. for r in twas]
    for r,p in zip(twas,holm(mel_inputs),strict=True):
        r["melbourne_holm_4"] = p if r["melbourne_probes"] else ""
    write_tsv(args.out/"koller_uterus.tsv", list(twas[0]), twas)
    (args.out/"omnibus.json").write_text(json.dumps(summary, indent=2, allow_nan=False)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ["finland", "rome", "koller", "melbourne", "definitions", "out"]:
        parser.add_argument("--"+flag, type=Path, required=True)
    main(parser.parse_args())
