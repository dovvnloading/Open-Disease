"""Measure coverage and targeted evidence in Fung et al. Table S20.

Only aggregate probe and gene statistics are written. Install openpyxl to
read the cached source workbook. See the frozen Melbourne plan.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

from cross_cohort_discovery import holm, write_tsv
from open_disease.geo_microarray import sha256_file

SOURCE_HASH = "b83a3d32a5882f750d978f4de1a6473fdd5d7627dbccf7d7c1774162e104e8e3"
TARGETS = {"ABO": (28, "higher"), "KIF3C": (3797, "lower"), "ZNF835": (90485, "lower"),
           "GLTP": (51228, "higher"), "F11": (2160, "lower")}
TITLE = ("Table S20. Differential expression of genes expressed in >90% of samples. "
         "Differentially expressed genes in Endometriosis Cases and Controls, ranked by p-value")
HEADER = ("PROBE_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "FC",
          "ILMN_GENE", "CHROMOSOME")
PROBE_FIELDS = ("symbol", "entrez_id", "probe_id", "source_logFC", "two_sided_p")
GENE_FIELDS = ("symbol", "entrez_id", "finnish_direction", "s20_probes", "measurement", "gene_p",
               "holm_input_p", "holm_p", "two_sided_result", "direction_agreement")
PROBE_ID = re.compile(r"ILMN_[1-9][0-9]*")


def check_nominations(lock):
    genes = [(r["id"], r["direction"]) for r in lock["candidates"] if r["family"] == "gene"]
    if genes != [(str(gene), direction) for gene, direction in TARGETS.values()]:
        raise ValueError("Source candidate IDs or directions disagree with the locked family")


def candidate_probes(rows, expected_count=15262):
    title = next(rows, ())
    header = next(rows, ())
    if (not title or title[0] != TITLE or tuple(header[:len(HEADER)]) != HEADER
            or any(value is not None for value in header[len(HEADER):])):
        raise ValueError("Unexpected Table S20 title or header")
    seen, selected = set(), []
    for row in rows:
        if len(row) < len(HEADER):
            raise ValueError("Unexpected Table S20 row width")
        probe = row[0]
        if not isinstance(probe, str) or not PROBE_ID.fullmatch(probe) or probe in seen:
            raise ValueError("Invalid or repeated probe in Table S20")
        seen.add(probe)
        symbol = row[8].strip() if isinstance(row[8], str) else ""
        if symbol not in TARGETS:
            continue
        coef, pv = row[1], row[4]
        if (type(coef) not in (float, int) or not math.isfinite(coef) or
                type(pv) not in (float, int) or not 0 <= pv <= 1):
            raise ValueError("Invalid coefficient or nominal probability for a matching probe")
        selected.append(dict(symbol=symbol, entrez_id=TARGETS[symbol][0], probe_id=probe,
                             source_logFC=coef, two_sided_p=pv))
    if len(seen) != expected_count:
        raise ValueError("Table S20 does not contain the frozen high-detection probe universe")
    return selected


def summarize(probes):
    genes = []
    for symbol, (gene, direction) in TARGETS.items():
        observed = [p["two_sided_p"] for p in probes if p["symbol"] == symbol]
        adjusted = min(1., len(observed) * min(observed)) if observed else None
        genes.append(dict(symbol=symbol, entrez_id=gene, finnish_direction=direction,
                          s20_probes=len(observed), measurement="measured" if observed else "unmeasured",
                          gene_p=adjusted, holm_input_p=adjusted if observed else 1.))
    for row, adjusted in zip(genes, holm([g["holm_input_p"] for g in genes]), strict=True):
        row["holm_p"] = float(adjusted) if row["s20_probes"] else None
        row["direction_agreement"] = "source_orientation_unresolved" if row["s20_probes"] else "unmeasured"
        if row["s20_probes"] == 0:
            row["two_sided_result"] = "untested"
        else:
            row["two_sided_result"] = "supported" if adjusted <= .05 else "did_not_pass"
    return genes


def run(workbook, lock_path, output):
    from openpyxl import load_workbook

    check_nominations(json.loads(lock_path.read_text(encoding="utf-8")))
    if sha256_file(workbook) != SOURCE_HASH:
        raise ValueError("Unexpected Melbourne workbook checksum")
    book = load_workbook(workbook, data_only=True, read_only=True)
    try:
        sheet = book["TableS20"]
        if sheet.max_row != 15264 or sheet.max_column != 16:
            raise ValueError("Unexpected TableS20 dimensions")
        probes = candidate_probes(iter(sheet.values))
    finally:
        book.close()
    order = {name: i for i, name in enumerate(TARGETS)}
    probes.sort(key=lambda r: (order[r["symbol"]], int(r["probe_id"].split("_")[-1])))
    genes = summarize(probes)
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "melbourne_genes.tsv", GENE_FIELDS, genes)
    write_tsv(output / "melbourne_probes.tsv", PROBE_FIELDS, probes)
    return genes, probes


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--nominations", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    genes, probes = run(args.workbook, args.nominations, args.out)
    print(f"Validated 15,262 probes; {len(probes)} matching probes, "
          f"{sum(g['s20_probes'] > 0 for g in genes)} measured genes in a family of {len(genes)}.")
