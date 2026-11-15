"""Replay the five-gene lesion panel on the deposited GSE282532 transcript FPKMs."""
import argparse
import csv
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
from openpyxl import load_workbook
from scipy.stats import rankdata

from genetic_anchor_tissue import ROOT, holm, one
from lesion_transfer import summarize
from open_disease.geo_microarray import sha256_file
from swiss_epithelial_mapping import HGNC_FILE, HGNC_SHA256


BOOK = "GSE282532_mRNA_FPKM.xlsx"
MATRIX = "GSE282532_series_matrix.txt.gz"
PANEL = ROOT / "evidence/lesion-transfer/panel.json"
PINS = {
    BOOK: "fd2be96816176f62d3ad2632b7bb81f8ab383a9663a3259c8222e1336c532911",
    MATRIX: "5fe3ea41dcc26c85a01c749e2c2a29a82aff6d87ef6ea4b2c4241d28dd2b5564",
    HGNC_FILE: HGNC_SHA256,
    PANEL.name: "4c634a65d648a83af8418e7ad1bf50e03d1a83982d0265d47d9ccbd06aba0114",
}


def checked(path):
    if sha256_file(path) != PINS[path.name]:
        raise ValueError(f"Input digest changed: {path.name}")
    return path


def parse_pairs(path, columns):
    header = defaultdict(list)
    with gzip.open(path, "rt", newline="") as f:
        for row in csv.reader(f, delimiter="\t"):
            if row:
                header[row[0]].append(row[1:])
    title = one(header, "!Sample_title")
    tissue = one(header, "!Sample_source_name_ch1")
    library, alias = header["!Sample_description"]
    if len(set(columns)) != 10 or set(columns) != set(alias):
        raise ValueError("GEO sample descriptions and workbook columns disagree")
    pairs = defaultdict(dict)
    for t, site, lib, a in zip(title, tissue, library, alias, strict=True):
        title_match = re.fullmatch(r"(EM|OMA)_(\d+)", t)
        alias_match = re.fullmatch(r"(EM|OSA)_([A-Z]+)", a)
        if not title_match or not alias_match or lib != f"Library name: {t}":
            raise ValueError("Unrecognized source-defined library key")
        uterus = site == "eutopic endometrium"
        if site not in {"ectopic endometrium", "eutopic endometrium"} or not (
            uterus == (title_match[1] == "EM") == (alias_match[1] == "EM")
        ):
            raise ValueError("Conflicting tissue labels")
        key = title_match[2], alias_match[2]
        if site in pairs[key]:
            raise ValueError("Duplicate tissue within a donor")
        pairs[key][site] = columns.index(a)
    if len(pairs) != 5 or len({k[0] for k in pairs}) != 5 or len({k[1] for k in pairs}) != 5:
        raise ValueError("The two pairing keys disagree")
    return ([p["ectopic endometrium"] for p in pairs.values()],
            [p["eutopic endometrium"] for p in pairs.values()])


def panel_ids(path, panel):
    ids = {str(r["entrez"]): r["direction"] for r in panel}
    eligible = {}
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["entrez_id"] in ids and r["status"] == "Approved":
                gene = r["ensembl_gene_id"]
                if gene in eligible:
                    raise ValueError("Repeated annotation for panel gene")
                eligible[gene] = (r["symbol"], ids[r["entrez_id"]])
    if len(eligible) != len(panel):
        raise ValueError("A locked member lacked a unique HGNC mapping")
    return eligible


def analyze(directory, hgnc):
    targets = panel_ids(checked(hgnc), json.loads(checked(PANEL).read_text()))
    book = load_workbook(checked(directory / BOOK), read_only=True, data_only=True)
    try:
        if book.sheetnames != ["mRNA_FPKM"]:
            raise ValueError("Unexpected worksheets")
        rows = iter(book.active.values)
        header = next(rows)
        if header[:4] != ("Transcript_id", "Gene_id", "Official_Symbol", "Gene_type"):
            raise ValueError("Unexpected transcript annotation")
        lesion, uterus = parse_pairs(checked(directory / MATRIX), list(header[4:]))
        genes, seen, values = [], set(), []
        for transcript, gene, symbol, _, *expression in rows:
            if not re.fullmatch(r"ENSG\d+", str(gene)):
                continue
            if transcript in seen or not re.fullmatch(r"ENST\d+", str(transcript)):
                raise ValueError("Invalid or repeated transcript")
            if gene in targets and targets[gene][0] != symbol:
                raise ValueError("HGNC Ensembl ID and deposited symbol disagree")
            seen.add(transcript)
            genes.append(gene)
            values.append(expression)
    finally:
        book.close()
    values = np.asarray(values, dtype=float)
    if values.shape != (83064, 10) or not np.isfinite(values).all() or (values < 0).any():
        raise ValueError("Invalid or missing transcript FPKM")
    ids, group = np.unique(genes, return_inverse=True)
    totals = np.zeros((len(ids), len(header) - 4))
    np.add.at(totals, group, values)
    percentile = 100 * (rankdata(totals, axis=0) - 0.5) / len(ids)
    transcript_pct = 100 * (rankdata(values, axis=0) - 0.5) / len(values)
    output, primary, alternative = {}, [], []
    gene_array = np.array(genes)
    for gene, (symbol, direction) in targets.items():
        index = np.flatnonzero(ids == gene).item()
        selected = gene_array == gene
        ranked = percentile[index]
        delta = ranked[lesion] - ranked[uterus]
        oriented = delta * direction
        primary.append(oriented)
        alternate = transcript_pct[selected].mean(axis=0)
        alternative.append((alternate[lesion] - alternate[uterus]) * direction)
        result = summarize(oriented)
        output[symbol] = dict(ensembl=gene, transcripts=int(selected.sum()), direction=direction,
                              difference_pp=float(delta.mean()), concordant=result["higher"],
                              one_sided_p=result["one_sided_p"])
    for symbol, adj in zip(output, holm([r["one_sided_p"] for r in output.values()]), strict=True):
        output[symbol]["holm_p"] = float(adj)
    return {
        "study": "GSE282532", "pairs": len(lesion), "transcripts": len(values), "genes": len(ids),
        "panel_genes": output, "primary_pp": summarize(np.mean(primary, axis=0)),
        "transcript_first_pp": summarize(np.mean(alternative, axis=0)),
        "input_sha256": PINS,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("hgnc", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze(args.data, args.hgnc), indent=2, allow_nan=False))
