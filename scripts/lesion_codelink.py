"""Evaluate the frozen lesion panel in a third cohort on a CodeLink microarray."""

import argparse
import csv
import gzip
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from cross_cohort_discovery import holm
from genetic_anchor_tissue import paired_summary
from genome_rome import read_annotation
from open_disease.geo_microarray import sha256_file

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "evidence/lesion-transfer/panel.json"
PINS = {
    "panel.json": "4c634a65d648a83af8418e7ad1bf50e03d1a83982d0265d47d9ccbd06aba0114",
    "GPL2895.annot.gz": "91e4f679c2a99e6274583b98906c28d4702146e85c63f26b927dc183d0041ce9",
    "GSE5108_series_matrix.txt.gz": "d76b8615d51ceec2a9a7009821b10c3ea95d72cc1cf86bc7a34b8941c8646f11",
    "GDS3092.soft.gz": "7e612e32e2a56a1b0eef5076ce3734792e1b754917d02af8d239bddd98c1f1a2",
}


def read_curated(path):
    groups = []
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("^SUBSET = "):
                groups.append({})
            elif line.startswith("!subset_"):
                key, value = line.strip().split(" = ", 1)
                groups[-1][key.removeprefix("!subset_")] = value
            elif line.startswith("!dataset_table_begin"):
                break
    pairs = [r["sample_id"].split(",") for r in groups if r["type"] == "individual"]
    if len(pairs) != 11 or any(len(pair) != 2 for pair in pairs):
        raise ValueError("Expected eleven two-specimen GEO individual subsets")
    return pairs


def load_design(header, paired_ids):
    cols = header["!Sample_geo_accession"]
    if len(set(cols)) != 22 or set(cols) != {i for pair in paired_ids for i in pair}:
        raise ValueError("Curated and submitted sample inventories differ")
    pairs, sites, mismatch = [], [], []
    for group in paired_ids:
        indices = {
            header["!Sample_source_name_ch1"][cols.index(x)]: cols.index(x)
            for x in group
        }
        if set(indices) != {"endometriosis", "endometrium"}:
            raise ValueError("Curated individual lacks one array from each tissue")
        i, j = indices["endometriosis"], indices["endometrium"]
        ids = []
        for k, tissue in [(i, "Endometriosis"), (j, "Endometrium")]:
            m = re.fullmatch(tissue + r" (\d+)(B?) WG", header["!Sample_title"][k])
            if not m:
                raise ValueError("Unrecognized submitted title")
            ids.append(m.groups())
        if ids[0][0] != ids[1][0]:
            raise ValueError("Donor numbers differ within curated pair")
        mismatch.append(ids[0] != ids[1])
        diagnosis = header["!Sample_characteristics_ch1"][i]
        sites.append(
            "ovarian"
            if "endometrioma" in diagnosis
            else (
                "peritoneal"
                if "peritoneal endometriosis" in diagnosis
                else "unspecified"
            )
        )
        pairs.append((i, j))
    if sum(mismatch) != 1 or Counter(sites) != {
        "ovarian": 6,
        "peritoneal": 3,
        "unspecified": 2,
    }:
        raise ValueError("Site mix or pairing audit changed")
    return np.array(pairs), np.array(sites), np.array(mismatch)


def analyze(directory):
    paths = {
        name: (PANEL if name == "panel.json" else directory / name) for name in PINS
    }
    if any(sha256_file(paths[k]) != d for k, d in PINS.items()):
        raise ValueError("Changed source or panel lock")
    panel = json.loads(PANEL.read_text())
    mapping = read_annotation(paths["GPL2895.annot.gz"])
    curated = read_curated(paths["GDS3092.soft.gz"])
    with gzip.open(paths["GSE5108_series_matrix.txt.gz"], "rt") as f:
        rows = csv.reader(f, delimiter="\t")
        header = {}
        for row in rows:
            if row == ["!series_matrix_table_begin"]:
                break
            if row:
                header[row[0]] = row[1:]
        if next(rows) != ["ID_REF", *header["!Sample_geo_accession"]]:
            raise ValueError("Source matrix columns differ from header")
        pairs, sites, mismatch = load_design(header, curated)
        features, values = [], []
        for row in rows:
            if row == ["!series_matrix_table_end"]:
                break
            features.append(row[0])
            values.append([float(x) for x in row[1:]])
        else:
            raise ValueError("Missing matrix terminator")
    v = np.array(values)
    if (
        v.shape != (54359, 22)
        or len(set(features)) != len(v)
        or not np.isfinite(v).all()
    ):
        raise ValueError("Unverified reference rows")
    ids = np.array([mapping[p] for p in features])
    ranks = (rankdata(v, axis=0) - 0.5) * (100 / len(v))
    deltas, results, ps = [], {}, []
    for gene in panel:
        key, direction = gene["entrez"], gene["direction"]
        selected = ids == key
        if not selected.any():
            raise ValueError(f"Frozen panel gene {key} lacks an eligible probe")
        scores = ranks[selected].mean(axis=0)
        diff = scores[pairs[:, 0]] - scores[pairs[:, 1]]
        summary = paired_summary(diff)
        tail = summary["upper"] if direction == 1 else summary["lower_tail"]
        ps.append(tail / summary["assignments"])
        results[key] = dict(
            effect_pp=summary["effect_pp"],
            eligible_probes=int(selected.sum()),
            p_directional=ps[-1],
            concordant_pairs=int((diff * direction > 0).sum()),
        )
        deltas.append(diff * direction)
    for row, adj in zip(results.values(), holm(ps), strict=True):
        row["holm"] = float(adj)
    score = np.mean(deltas, axis=0)
    summary = paired_summary(score)
    summary["p_directional"] = summary["upper"] / summary["assignments"]
    summary["omit_suffix_mismatch"] = paired_summary(score[~mismatch])
    summary["by_site"] = {
        site: dict(
            n=int((sites == site).sum()),
            effect_pp=float(score[sites == site].mean()),
            positive=int((score[sites == site] > 0).sum()),
        )
        for site in sorted(set(sites))
    }
    return dict(
        pairs=11,
        reference_features=len(ids),
        annotation_features=len(mapping),
        panel=summary,
        genes=results,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    print(json.dumps(analyze(p.parse_args().directory), indent=2, allow_nan=False))
