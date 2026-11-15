"""Independently verify the complete CodeLink panel output from pinned GEO sources."""

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import numpy as np

from verify_genetic_anchor_pairs import compare, exact_pair_test

PINS = {
    "panel.json": "4c634a65d648a83af8418e7ad1bf50e03d1a83982d0265d47d9ccbd06aba0114",
    "GPL2895.annot.gz": "91e4f679c2a99e6274583b98906c28d4702146e85c63f26b927dc183d0041ce9",
    "GSE5108_series_matrix.txt.gz": "d76b8615d51ceec2a9a7009821b10c3ea95d72cc1cf86bc7a34b8941c8646f11",
    "GDS3092.soft.gz": "7e612e32e2a56a1b0eef5076ce3734792e1b754917d02af8d239bddd98c1f1a2",
}


def study_design(path, columns, meta):
    subsets = []
    study_fields = {}
    with gzip.open(path, "rt") as source:
        for line in source:
            if line.startswith("!dataset_table_begin"):
                if Counter(next(source).rstrip().split("\t")[2:]) != Counter(columns):
                    raise ValueError(
                        "Curated expression columns disagree with submitted accessions"
                    )
                break
            key, _, value = line.rstrip().partition(" = ")
            if key == "^SUBSET":
                subsets.append({})
            elif key.startswith("!subset_"):
                subsets[-1][key.removeprefix("!subset_")] = value
            elif key.startswith("!dataset_"):
                study_fields[key] = value
    if not (
        study_fields["!dataset_platform"] == "GPL2895"
        and study_fields["!dataset_reference_series"] == "GSE5108"
        and study_fields["!dataset_sample_count"] == "22"
    ):
        raise ValueError("Curated dataset identifiers disagree")
    tissues = {
        x["description"]: set(x["sample_id"].split(","))
        for x in subsets
        if x["type"] == "disease state"
    }
    people = [x["sample_id"].split(",") for x in subsets if x["type"] == "individual"]
    if (
        len(people) != 11
        or set(tissues) != {"endometriosis", "control"}
        or Counter(s for pair in people for s in pair) != Counter(columns)
    ):
        raise ValueError("GDS donor or tissue inventory differs")
    pairs, sites, mismatch = [], [], []
    for pair in people:
        (lesion,) = set(pair) & tissues["endometriosis"]
        (uterus,) = set(pair) & tissues["control"]
        i, j = columns.index(lesion), columns.index(uterus)
        if (meta["!Sample_source_name_ch1"][i], meta["!Sample_source_name_ch1"][j]) != (
            "endometriosis",
            "endometrium",
        ):
            raise ValueError("Curated and submitted tissue assignments differ")
        a = re.fullmatch(r"Endometriosis (\d+)(B?) WG", meta["!Sample_title"][i])
        b = re.fullmatch(r"Endometrium (\d+)(B?) WG", meta["!Sample_title"][j])
        if not a or not b or a[1] != b[1]:
            raise ValueError("Numeric donor identity check failed")
        mismatch.append(a[2] != b[2])
        diagnosis = meta["!Sample_characteristics_ch1"][i]
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
        raise ValueError("Unexpected lesion-site or suffix-mismatch audit")
    return pairs, sites, mismatch


def verify(directory, lock, result_path):
    for name, digest in PINS.items():
        with (lock if name == "panel.json" else directory / name).open("rb") as stream:
            if hashlib.file_digest(stream, "sha256").hexdigest() != digest:
                raise ValueError("Changed pinned file")
    panel = json.loads(lock.read_text())
    with gzip.open(directory / "GPL2895.annot.gz", "rt", newline="") as stream:
        for line in stream:
            if line.rstrip() == "!platform_table_begin":
                break
        annotation = csv.DictReader(stream, delimiter="\t")
        mapped, all_probes = defaultdict(set), set()
        for row in annotation:
            probe = row["ID"]
            if probe == "!platform_table_end":
                break
            if probe in all_probes:
                raise ValueError("Duplicate CodeLink annotation")
            all_probes.add(probe)
            if not probe.startswith("AFFX") and re.fullmatch(
                r"[a-zA-Z0-9][a-zA-Z0-9_.-]*", row["Gene symbol"]
            ):
                mapped[row["Gene ID"]].add(probe)
    with gzip.open(
        directory / "GSE5108_series_matrix.txt.gz", "rt", newline=""
    ) as stream:
        fields = csv.reader(stream, delimiter="\t")
        meta = {}
        for row in fields:
            if row == ["!series_matrix_table_begin"]:
                break
            if row and row[0].startswith("!Sample_"):
                meta[row[0]] = row[1:]
        header = next(fields)
        columns = header[1:]
        if columns != meta["!Sample_geo_accession"] or len(set(columns)) != 22:
            raise ValueError("Series column mismatch")
        values = {}
        for row in fields:
            if row == ["!series_matrix_table_end"]:
                break
            if row[0] in values:
                raise ValueError("Duplicate matrix probe")
            values[row[0]] = [float(v) for v in row[1:]]
    if (
        len(values) != 54359
        or len(all_probes) != 56448
        or not set(values) <= all_probes
    ):
        raise ValueError("Incorrect feature universe")
    pairs, sites, mismatch = study_design(directory / "GDS3092.soft.gz", columns, meta)
    matrix = np.array(list(values.values()))
    if matrix.shape != (54359, 22) or not np.isfinite(matrix).all():
        raise ValueError("Unaudited signal matrix")
    ordered = np.sort(matrix, axis=0)
    genes, oriented, counts = {}, {}, {}
    for gene in panel:
        key = str(gene["entrez"])
        selected = mapped[key] & values.keys()
        counts[key] = len(selected)
        nums = [0] * len(pairs)
        for probe in selected:
            r = [
                int(
                    np.searchsorted(col, val, "left")
                    + np.searchsorted(col, val, "right")
                )
                for col, val in zip(ordered.T, values[probe], strict=True)
            ]
            for k, (i, j) in enumerate(pairs):
                nums[k] += r[i] - r[j]
        scale = 50 / (len(values) * len(selected))
        raw = exact_pair_test(nums, scale)
        direction = gene["direction"]
        oriented[key] = [d * direction for d in nums]
        tail = raw["upper"] if direction == 1 else raw["lower_tail"]
        genes[key] = dict(
            effect_pp=raw["effect_pp"],
            eligible_probes=len(selected),
            p_directional=Fraction(tail, raw["assignments"]),
            concordant_pairs=sum(d > 0 for d in oriented[key]),
        )
    order = sorted(genes, key=lambda gene: genes[gene]["p_directional"])
    for i, key in enumerate(order):
        genes[key]["holm"] = min(
            1,
            max(
                (len(order) - j) * genes[order[j]]["p_directional"]
                for j in range(i + 1)
            ),
        )
    common = math.lcm(*counts.values())
    total = [
        sum(oriented[key][i] * (common // counts[key]) for key in genes)
        for i in range(len(pairs))
    ]
    scale = 50 / (len(values) * len(genes) * common)
    summary = exact_pair_test(total, scale)
    summary["p_directional"] = Fraction(summary["upper"], summary["assignments"])
    summary["omit_suffix_mismatch"] = exact_pair_test(
        [d for d, flag in zip(total, mismatch, strict=True) if not flag], scale
    )
    by_site = {}
    for site in set(sites):
        group = [d for s, d in zip(sites, total, strict=True) if s == site]
        by_site[site] = dict(
            n=len(group),
            effect_pp=scale * sum(group) / len(group),
            positive=sum(d > 0 for d in group),
        )
    summary["by_site"] = by_site
    observed = dict(
        pairs=len(pairs),
        reference_features=len(values),
        annotation_features=len(all_probes),
        panel=summary,
        genes=genes,
    )
    compare(observed, json.loads(result_path.read_text()))
    return dict(
        paired_score_estimate_verified=True,
        genes_verified=len(genes),
        all_reported_fields_verified=True,
        same_direction_pairs=summary["higher"],
        total_pairs=len(pairs),
        upper=summary["upper"],
        assignments=summary["assignments"],
        curated_tissue_subsets_checked=True,
        curated_individual_subsets_checked=len(pairs),
        suffix_mismatch_omission_checked=True,
        subgroup_aggregates_checked=len(by_site),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ("directory", "lock", "results"):
        parser.add_argument(arg, type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.directory, args.lock, args.results), indent=2))
