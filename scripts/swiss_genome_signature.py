"""Test the Finnish signed expression profile in proliferative Swiss glandular donors."""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.stats import rankdata

from audit_gse266265_epithelial import DIGEST, WORKBOOK_DIGEST, verify, write_sample_map
from cross_cohort_discovery import assignment_design
from genome_signature_transfer import DISCOVERY_SHA, score_profiles, weights_and_candidates
from swiss_epithelial_mapping import (
    HGNC_FILE, HGNC_SHA256, HGNC_SIZE, LABEL_COUNT, LABEL_SHA256,
    LOCK_FILE, LOCK_SHA256, candidate_projection, crosswalk, write_projection,
)
from swiss_global_interval import allocation_contrasts, shift_limits
from open_disease.geo_microarray import sha256_file
from verify_swiss_glandular import load_pathways, read_integers


def glandular_ranks(stream):
    """Sum each native sparse column into its donor. Retain all zero features."""
    with stream.open("rb") as handle:
        features, donors, cells = map(int, read_integers(handle, 3))
        if (features, donors, cells) != (29034, 13, 3564):
            raise ValueError("Unexpected Swiss glandular data dimensions")
        labels = read_integers(handle, donors)
        batches = read_integers(handle, donors)
        sizes = read_integers(handle, donors)
        if sorted(set(labels)) != [0, 1] or labels.sum() != 6 or min(sizes) < 20:
            raise ValueError("Unexpected donor allocation")
        counts = np.zeros((donors, features), dtype=np.int64)
        seen = np.zeros(donors, dtype=int)
        for _ in range(cells):
            donor, nnz = map(int, read_integers(handle, 2))
            if not 0 <= donor < donors or not 0 <= nnz <= features:
                raise ValueError("Invalid sparse column header")
            rows = read_integers(handle, nnz)
            values = read_integers(handle, nnz)
            if np.any(rows < 0) or np.any(rows >= features) or np.any(np.diff(rows) <= 0) or np.any(values < 0):
                raise ValueError("Invalid sparse count column")
            counts[donor, rows] += values
            seen[donor] += 1
        if handle.read(1) or not np.array_equal(seen, sizes):
            raise ValueError("Swiss sparse column counts disagree")
    ranks = 100 * (rankdata(counts, axis=1, method="average") - 0.5) / features
    return ranks, labels.astype(bool), batches


def check_prior_pathways(ranks, cases, batches, pathways):
    observed, allocations, _ = assignment_design(cases, batches)
    if len(allocations) != 144 or len(set(batches)) != 4:
        raise ValueError("Unexpected Swiss within-batch allocation space")
    sentinel = np.mean([ranks[:, idx].mean(axis=1) for idx in pathways.values()], axis=0)
    sentinel_effect = float(observed @ sentinel)
    sentinel_extreme = int((abs(allocations @ sentinel) >= abs(sentinel_effect) - 1e-10).sum())
    if abs(sentinel_effect - 0.7092429609566739) > 1e-9 or sentinel_extreme != 24:
        raise ValueError("Previously published Swiss five-pathway contrast failed to reproduce")
    return {"effect_pp": sentinel_effect, "extreme": sentinel_extreme}


def summarize(gene_ids, ranks, cases, batches, sentinel, reference, weights, candidates):
    observed, allocations, _ = assignment_design(cases, batches)
    excluded = set()
    for item in candidates:
        excluded.update(item["entrez_ids"] if item["family"] == "pathway" else [int(item["id"])])
    n, n_without, scores = score_profiles(gene_ids, ranks, reference, weights, excluded)
    effects = observed @ scores
    extremes = (np.abs(allocations @ scores) >= np.abs(effects) - 1e-10).sum(axis=0)
    donor_scores = [Fraction(float(x)) for x in scores[:, 0]]
    actual, alternatives = allocation_contrasts(donor_scores, cases.astype(int).tolist(), batches.tolist())
    tolerance = Fraction(1, 10**10)
    extreme_check = sum(abs(a) + tolerance >= abs(actual) for a, _ in alternatives)
    if extreme_check != extremes[0] or abs(float(actual) - effects[0]) > 1e-10:
        raise ValueError("Independent exact contrast and permutation enumeration disagree")
    low, high, _ = shift_limits(actual, alternatives)
    omitted = []
    for batch in np.unique(batches):
        keep = batches != batch
        design, _, _ = assignment_design(cases[keep], batches[keep])
        omitted.append(float(design @ scores[keep, 0]))
    shared = np.isin(reference, gene_ids)
    return {
        "cells": 3564, "donors": len(cases), "cases": int(cases.sum()), "controls": int((~cases).sum()),
        "batches": len(omitted), "assignments": len(allocations), "features_ranked": LABEL_COUNT,
        "reference_genes": len(reference), "measured_genes": n,
        "absolute_weight_fraction_measured": float(abs(weights[shared]).sum() / abs(weights).sum()),
        "primary_score_difference_pp": float(effects[0]), "primary_exact_two_sided_extreme": int(extremes[0]),
        "primary_exact_two_sided_p": float(extremes[0] / len(allocations)),
        "constant_shift_95_interval": [None if x is None else float(x) for x in (low, high)],
        "leave_one_batch_out_range": [min(omitted), max(omitted)],
        "genes_after_removing_nominations": n_without,
        "without_nominations_score_difference_pp": float(effects[1]),
        "without_nominations_extreme": int(extremes[1]),
        "prior_five_pathway_sentinel": sentinel,
        "discovery_sha256": DISCOVERY_SHA, "source_rds_sha256": DIGEST,
        "source_clinical_sha256": WORKBOOK_DIGEST, "source_hgnc_sha256": HGNC_SHA256,
        "ordered_feature_labels_sha256": LABEL_SHA256,
    }


def _checked_profiles(rds, workbook, hgnc, rscript):
    verify(rds)
    if hgnc.stat().st_size != HGNC_SIZE or sha256_file(hgnc) != HGNC_SHA256 or hgnc.name != HGNC_FILE:
        raise ValueError("Unexpected HGNC source")
    if sha256_file(LOCK_FILE) != LOCK_SHA256:
        raise ValueError("Swiss sentinel nomination lock changed")
    candidates = json.loads(LOCK_FILE.read_text())["candidates"]
    # In repository use the adjacent base-R scripts. Temporary private files are removed on exit.
    r_scripts = Path(__file__).resolve().parent / "r"
    env = dict(os.environ, R_DEFAULT_PACKAGES="NULL", OPENBLAS_NUM_THREADS="1")
    with tempfile.TemporaryDirectory(dir=rds.parent) as private:
        directory = Path(private)
        labels_path = directory / "row_labels.txt"
        subprocess.run([rscript, r_scripts / "extract_swiss_feature_labels.R", rds, labels_path],
                       env=env, check=True)
        if sha256_file(labels_path) != LABEL_SHA256:
            raise ValueError("Unexpected ordered Swiss feature axis")
        labels = labels_path.read_text().splitlines()
        if len(labels) != LABEL_COUNT:
            raise ValueError("Unexpected Swiss feature count")
        with hgnc.open(encoding="utf-8", newline="") as handle:
            mapping, _ = crosswalk(csv.DictReader(handle, delimiter="\t"), labels)
        projection = directory / "projection.tsv"
        write_projection(projection, candidate_projection(candidates, mapping))
        pathways = load_pathways(projection)
        clinical, stream = directory / "clinical.tsv", directory / "columns.bin"
        write_sample_map(workbook, clinical)
        subprocess.run([rscript, r_scripts / "export_swiss_verification_columns.R", rds,
                        clinical, projection, stream], env=env, check=True)
        ranks, cases, batches = glandular_ranks(stream)
    sentinel = check_prior_pathways(ranks, cases, batches, pathways)
    genes = np.array(sorted(int(key) for key in mapping))
    rows = [mapping[str(gene)][1] - 1 for gene in genes]
    return (genes, ranks[:, rows], cases, batches), sentinel


def load_glandular_profiles(rds, workbook, hgnc, rscript="Rscript"):
    """Return sorted Entrez IDs, donor-by-gene rank percentiles, case flags and capture batches.

    Each donor's ranks use the entire original feature axis. Only uniquely
    mapped genes are returned. All inputs and the prior pathway sentinel are checked.
    """
    profiles, _ = _checked_profiles(Path(rds), Path(workbook), Path(hgnc), rscript)
    return profiles


def analyze(rds, workbook, hgnc, rscript):
    reference, weights, candidates = weights_and_candidates()
    profiles, sentinel = _checked_profiles(rds, workbook, hgnc, rscript)
    return summarize(*profiles, sentinel, reference, weights, candidates)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rds", type=Path)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("hgnc", type=Path)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    print(json.dumps(analyze(args.rds, args.workbook, args.hgnc, args.rscript), indent=2, allow_nan=False))
