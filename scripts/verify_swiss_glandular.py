"""Independently recompute Swiss glandular aggregate endpoints from native sparse counts.

Takes the pinned RDS, clinical workbook and locked HGNC membership projection.
Temporary cell and donor data are private and removed after each execution.
"""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
from itertools import combinations, product
import json
import os
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from audit_gse266265_epithelial import verify, write_sample_map
from open_disease.geo_microarray import sha256_file

PROJECTION_SHA256 = "6609ef4bb35608c7f791391b54344e0b0e9f77a44ed3004d95925aed2779a0c2"
TOLERANCE = Fraction(1, 10**10)


def load_pathways(path):
    if sha256_file(path) != PROJECTION_SHA256:
        raise ValueError("Unexpected locked membership projection")
    groups = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            if record["family"] == "pathway" and record["row_index"]:
                groups.setdefault(record['id'], []).append(int(record['row_index']) - 1)
    return groups


def read_integers(handle, length):
    result = np.fromfile(handle, dtype="<i4", count=length)
    if len(result) != length:
        raise ValueError("Incomplete private count stream")
    return result


def twice_rank_offset(vector, selected):
    # L observations strictly below and U at or below a value imply an
    # average rank (L+1+U)/2. Twice its offset from 0.5 equals L+U.
    ordered = np.sort(vector)
    values = vector[selected]
    return np.searchsorted(ordered, values, side="left") + np.searchsorted(ordered, values, side="right")


def donor_endpoints(numerators, lengths, features, cells):
    scores = []
    offset = 0
    for members in lengths:
        numerator = int(numerators[offset:offset + members].sum()) * 50
        scores.append(Fraction(numerator, features * cells * members))
        offset += members
    return [sum(scores) / len(scores), *scores]


def reconstruct(stream, pathways):
    selected = np.concatenate(list(pathways.values()))
    lengths = [len(positions) for positions in pathways.values()]
    with stream.open('rb') as handle:
        features, donors, cells = map(int, read_integers(handle, 3))
        if features != 29_034 or donors != 13 or not donors * 20 <= cells <= 20_549:
            raise ValueError("Unexpected source dimensions")
        labels = read_integers(handle, donors).astype(bool)
        batches = read_integers(handle, donors)
        cell_counts = read_integers(handle, donors)
        sums = np.zeros((donors, features), dtype=np.int64)
        cell_ranks = np.zeros((donors, len(selected)), dtype=np.int64)
        observed_counts = np.zeros(donors, dtype=np.int64)
        single = np.zeros(features, dtype=np.int32)
        for _ in range(cells):
            who, nnz = map(int, read_integers(handle, 2))
            if not 0 <= who < donors or not 0 <= nnz <= features:
                raise ValueError("Invalid CSC column")
            rows = read_integers(handle, nnz)
            values = read_integers(handle, nnz)
            # R checks unique ordered indices and nonnegative integer values.
            # Dense zeros here retain *all* absent features in the tie group.
            single[rows] = values
            sums[who] += single
            cell_ranks[who] += twice_rank_offset(single, selected)
            observed_counts[who] += 1
            single[rows] = 0
        if handle.read(1) or not np.array_equal(observed_counts, cell_counts):
            raise ValueError("Column counts or stream length disagree")
    primary = [donor_endpoints(twice_rank_offset(total, selected), lengths, features, 1)
               for total in sums]
    sensitivity = [donor_endpoints(ranks, lengths, features, int(n))
                   for ranks, n in zip(cell_ranks, cell_counts, strict=True)]
    return labels, batches, cells, primary, sensitivity


def donor_contrast(scores, donors, cases):
    assigned = set(cases)
    controls = [row for row in donors if row not in assigned]
    return sum(scores[row] for row in assigned) / len(assigned) - sum(scores[row] for row in controls) / len(controls)


def test_endpoints(labels, batches, primary, sensitivity, pathways):
    strata = [np.flatnonzero(batches == key).tolist() for key in sorted(set(batches))]
    observed = [[row for row in stratum if labels[row]] for stratum in strata]
    assignments = [list(combinations(rows, len(actual))) for rows, actual in zip(strata, observed, strict=True)]
    count = int(np.prod([len(choices) for choices in assignments]))
    if count != 144:
        raise ValueError("The frozen allocation space has changed")
    endpoints = {}
    for j, endpoint in enumerate(['GLOBAL', *pathways]):
        scores = [row[j] for row in primary]
        cell_scores = [row[j] for row in sensitivity]
        batchwise = [donor_contrast(scores, donors, cases) for donors, cases in zip(strata, observed, strict=True)]
        difference = sum(batchwise) / len(strata)
        omit = [(sum(batchwise) - delta) / (len(strata) - 1) for delta in batchwise]
        possible = [[donor_contrast(scores, donors, choice) for choice in choices]
                    for donors, choices in zip(strata, assignments, strict=True)]
        extreme = sum(abs(sum(deltas) / len(strata)) + TOLERANCE >= abs(difference) for deltas in product(*possible))
        cell_first = sum(donor_contrast(cell_scores, donors, cases)
                         for donors, cases in zip(strata, observed, strict=True)) / len(strata)
        endpoints[endpoint] = dict(difference=float(difference), extreme=extreme, p=extreme/count, holm=None,
                                   batch_omit_low=float(min(omit)), batch_omit_high=float(max(omit)),
                                   cell_first=float(cell_first))
    adjusted = 0
    ordered = sorted(pathways, key=lambda key: endpoints[key]['extreme'])
    for order, key in enumerate(ordered):
        adjusted = max(adjusted, endpoints[key]['extreme'] * (len(ordered) - order))
        endpoints[key]['holm'] = min(count, adjusted) / count
    return count, endpoints


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rds', type=Path)
    parser.add_argument('supplement', type=Path)
    parser.add_argument('projection', type=Path)
    parser.add_argument('--rscript', default='Rscript')
    args = parser.parse_args()
    verify(args.rds)
    pathways = load_pathways(args.projection)
    with tempfile.TemporaryDirectory(dir=args.rds.parent) as directory:
        clinical = Path(directory) / 'clinical.tsv'
        stream = Path(directory) / 'columns.bin'
        write_sample_map(args.supplement, clinical)
        env = dict(os.environ, R_DEFAULT_PACKAGES='NULL', OPENBLAS_NUM_THREADS='1')
        subprocess.run([args.rscript, str(Path(__file__).parent / 'r/export_swiss_verification_columns.R'),
                        str(args.rds), str(clinical), str(args.projection), str(stream)], env=env, check=True)
        labels, batches, cells, primary, sensitivity = reconstruct(stream, pathways)
    count, endpoints = test_endpoints(labels, batches, primary, sensitivity, pathways)
    print(json.dumps(dict(donors=len(labels), cases=int(sum(labels)), controls=int(sum(~labels)),
                          batches=len(set(batches)), cells=cells, assignments=count, endpoints=endpoints), indent=2))


if __name__ == '__main__':
    main()
