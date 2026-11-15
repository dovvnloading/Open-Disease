"""Evaluate a fixed, scan-disjoint matrisome ranking rule across uterine cohorts."""
import argparse
from collections import Counter
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from audit_matrisome_holdout import audit, read_keyed
from audit_us_bulk_sources import native_signature
from bulk_gse153739 import load_mumbai
from bulk_gse6364 import GPL_SHA, MATRIX_SHA, checked_platform
from cross_cohort_discovery import align
from genome_rome import load_rome
from open_disease.geo_microarray import sha256_file
from reconcile_gse51981_sources import CHAR, DIAGNOSIS, PATHOLOGY, GROUP_CODES, read_quick_view

HASHES = {"GPL570.annot.gz": GPL_SHA, "GSE6364.gz": MATRIX_SHA,
          "GSE51981.gz": "810e560c2b10dabb2f92dde03f64a2f9af4d72d51360f29c102569942dd84211",
          "all_matrisome_genes.json": "22a228b4a8802a3137fa7e77db35f165e7eb730f127d261a373c412c49065b35"}


def mapped_genes(annotation, symbols, mapping):
    by_symbol = {s: set() for s in symbols}
    with gzip.open(annotation, "rt") as stream:
        for line in stream:
            if line.strip() == "!platform_table_begin":
                break
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["ID"] == "!platform_table_end":
                break
            symbol, gene = row["Gene symbol"], mapping[row["ID"]]
            if symbol in by_symbol and gene:
                by_symbol[symbol].add(gene)
    return sorted({next(iter(ids)) for ids in by_symbol.values() if len(ids) == 1})


def array_ranks(matrix, mapping, genes, selected):
    # Rank against the entire audited platform before selecting matrisome features.
    with gzip.open(matrix, "rt") as stream:
        for line in stream:
            if line.strip() == "!series_matrix_table_begin":
                break
        reader = csv.reader(stream, delimiter="\t")
        header = next(reader)
        keys = [key for key in header[1:] if key in selected]
        if len(keys) != len(selected) or set(keys) != selected:
            raise ValueError("Matrix columns and selected arrays disagree")
        indices = [header.index(k) for k in keys]
        profiles = np.empty((len(keys), len(mapping)))
        probes = []
        for row in reader:
            if row == ["!series_matrix_table_end"]:
                break
            if len(row) != len(header) or len(probes) == len(mapping):
                raise ValueError("Malformed or excess platform feature")
            profiles[:, len(probes)] = [float(row[i]) for i in indices]
            probes.append(row[0])
    if set(probes) != mapping.keys() or len(probes) != len(mapping) or not np.isfinite(profiles).all():
        raise ValueError("Incomplete platform or nonfinite array signals")
    probe_genes = np.array([mapping[p] for p in probes])
    ranks = (rankdata(profiles, axis=1) - 0.5) * (100 / len(probes))
    return keys, np.column_stack([ranks[:, probe_genes == g].mean(axis=1) for g in genes])


def training_reference(profiles, cases, phases):
    contrasts = []
    for phase in sorted(set(phases)):
        cases_in_phase = (phases == phase) & cases
        controls_in_phase = (phases == phase) & ~cases
        if cases_in_phase.any() and controls_in_phase.any():
            contrasts.append(profiles[cases_in_phase].mean(axis=0) - profiles[controls_in_phase].mean(axis=0))
    if len(contrasts) != 3:
        raise ValueError("Insufficient phase overlap in reference")
    vector = np.mean(contrasts, axis=0)
    if np.linalg.norm(vector) == 0:
        raise ValueError("Zero reference vector")
    return vector / np.linalg.norm(vector)


def pair_comparison(scores, cases):
    comparison = scores[cases, None] - scores[~cases]
    higher, ties, pairs = int((comparison > 0).sum()), int((comparison == 0).sum()), comparison.size
    return {"cases": int(cases.sum()), "comparators": int((~cases).sum()), "higher": higher,
            "ties": ties, "pairs": pairs, "fraction": (higher + 0.5*ties)/pairs}


def run(directory, scans, rome, mumbai, hgnc):
    for file, checksum in HASHES.items():
        if sha256_file(directory / file) != checksum:
            raise ValueError(f"Unexpected benchmark input {file}")
    audit(scans)
    membership = read_keyed(scans / "coldata.tsv")
    clinical = read_quick_view(scans / "GSE51981.soft")
    groups = {}
    for key, record in clinical.items():
        f = dict(field.split(": ", 1) for field in record[CHAR])
        groups[key] = "case" if f[DIAGNOSIS] == "Endometriosis" else GROUP_CODES[f[DIAGNOSIS], f[PATHOLOGY]]
    train = {k for k, row in membership.items() if row["batch"] == "3" and groups[k] != "other"}
    test = {k for k, row in membership.items() if row["batch"] == "2"}
    test_images = {native_signature((scans / "GSE6364" / (key + ".part")).read_bytes())[1] for key in test}
    clean = {key for key in train if native_signature((scans / "GSE51981" / (key + ".part")).read_bytes())[1] not in test_images}
    mapping = checked_platform(directory / "GPL570.annot.gz")
    symbols = json.loads((directory / "all_matrisome_genes.json").read_text())["symbols"]
    affy_genes = mapped_genes(directory / "GPL570.annot.gz", symbols, mapping)
    rg, rp, rc, _ = load_rome(rome)
    mg, mp, mc, _ = load_mumbai(mumbai, hgnc)
    genes = np.array(sorted(set(affy_genes) & set(rg) & set(mg)), dtype=int)
    train_keys, xp = array_ranks(directory / "GSE51981.gz", mapping, genes, train)
    test_keys, xt = array_ranks(directory / "GSE6364.gz", mapping, genes, test)
    train_cases = np.array([groups[k] == "case" for k in train_keys])
    test_cases = np.array([membership[k]["condition"] == "endometriosis" for k in test_keys])
    train_phase, test_phase = [np.array([membership[k]["phase"] for k in keys]) for keys in (train_keys, test_keys)]
    results = {}
    for name, keep in (("scan_disjoint", clean), ("with_reused_scans", train)):
        mask = np.array([k in keep for k in train_keys])
        weights = training_reference(xp[mask], train_cases[mask], train_phase[mask])
        score = xt @ weights
        rows = {phase: pair_comparison(score[test_phase == phase], test_cases[test_phase == phase])
                for phase in sorted(set(test_phase))}
        pairs, higher, ties = [sum(row[k] for row in rows.values()) for k in ("pairs", "higher", "ties")]
        results[name] = {"reference_cases": int(train_cases[mask].sum()),
                         "reference_comparators": int((~train_cases[mask]).sum()),
                         "reference_counts_by_phase": {phase: dict(Counter(groups[k] for k in keep if membership[k]["phase"] == phase))
                                                       for phase in sorted(set(train_phase))},
                         "GSE6364_by_phase": rows, "GSE6364_pooled_within_phase": {"pairs": pairs,
                             "higher": higher, "ties": ties, "fraction": (higher + 0.5*ties)/pairs},
                         "Rome": pair_comparison(align(rg, rp, genes) @ weights, rc),
                         "Mumbai": pair_comparison(align(mg, mp, genes) @ weights, mc)}
    return {"reference_symbols": len(symbols), "GPL570_mapped_genes": len(affy_genes),
            "common_entrez_ids": genes.tolist(), "results": results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("directory", "scans", "rome", "mumbai", "hgnc"):
        parser.add_argument(key, type=Path)
    args = parser.parse_args()
    print(json.dumps(run(**vars(args)), sort_keys=True, indent=2))
