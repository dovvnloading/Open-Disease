"""Project Finnish genome-wide contrasts into two separately recruited biopsy cohorts."""
from __future__ import annotations

import argparse
import csv
from fractions import Fraction
import json
from pathlib import Path

import numpy as np

from bulk_gse153739 import load_mumbai
from cross_cohort_discovery import align, assignment_design, holm
from genome_rome import load_rome
from open_disease.geo_microarray import sha256_file
from rome_scan_design import scan_design
from rome_hybridization_controls import hybridization_controls
from swiss_epithelial_mapping import HGNC_FILE, LOCK_FILE, LOCK_SHA256
from swiss_global_interval import allocation_contrasts, shift_limits

DISCOVERY = LOCK_FILE.with_name('discovery_gene.tsv')
DISCOVERY_SHA = '9df65b29615ec56767e9abaaf3844f8ebde5e21baf8cdab193ac24857a170da5'


def weights_and_candidates():
    if sha256_file(DISCOVERY) != DISCOVERY_SHA or sha256_file(LOCK_FILE) != LOCK_SHA256:
        raise ValueError('Finnish discovery weights or nomination membership changed')
    with DISCOVERY.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    genes = np.array([int(r['id']) for r in rows])
    weight = np.array([float(r['statistic']) for r in rows])
    if len(genes) != 12044 or not (np.diff(genes) > 0).all() or not np.isfinite(weight).all():
        raise ValueError('Unexpected Finnish discovery universe')
    return genes, weight, json.loads(LOCK_FILE.read_text())['candidates']


def score_profiles(genes, percentiles, reference, weight, excluded):
    shared, idx, ref_idx = np.intersect1d(genes, reference, return_indices=True)
    if not len(shared) or len(np.unique(genes)) != len(genes):
        raise ValueError('Missing overlap or duplicate gene identity')
    data = percentiles[:, idx]
    w = weight[ref_idx]
    retained = ~np.isin(shared, list(excluded))
    def project(x, v):
        norm = np.sum(np.abs(v))
        if norm == 0:
            raise ValueError('No nonzero training weights in measured set')
        return x @ (v / norm)
    return (len(shared), int(retained.sum()),
            np.column_stack([project(data, w), project(data[:, retained], w[retained]), project(data, np.sign(w))]))


def transfer(genes, percentiles, cases, blocks, reference, weight, excluded):
    n, retained, scores = score_profiles(genes, percentiles, reference, weight, excluded)
    mixed = [b for b in np.unique(blocks) if len(np.unique(cases[blocks == b])) == 2]
    keep = np.isin(blocks, mixed)
    cases, blocks, scores = cases[keep], blocks[keep], scores[keep]
    obs, allocations, _ = assignment_design(cases, blocks)
    differences = obs @ scores
    # The same participant label allocation is used for every gene in the score.
    null = allocations @ scores[:, 0]
    extreme = int((np.abs(null) >= abs(differences[0]) - 1e-10).sum())
    y = [Fraction(float(s)) for s in scores[:, 0]]
    observed, contrasts = allocation_contrasts(y, cases.astype(int).tolist(), blocks.tolist())
    low, high, _ = shift_limits(observed, contrasts)
    omission = []
    for i in range(len(cases)):
        eligible = np.arange(len(cases)) != i
        if len(np.unique(cases[eligible & (blocks == blocks[i])])) == 2:
            delete, _, _ = assignment_design(cases[eligible], blocks[eligible])
            omission.append(float(delete @ scores[eligible, 0]))
    return {'measured_genes': n, 'genes_after_omitting_nominations': retained,
            'cases': int(cases.sum()), 'controls': int((~cases).sum()), 'strata': len(mixed),
            'score_difference_pp': float(differences[0]), 'exact_two_sided_extreme': extreme,
            'assignments': len(null), 'p_two_sided': extreme / len(null),
            'constant_shift_95_interval': [None if x is None else float(x) for x in (low, high)],
            'leave_one_biopsy_out_range': [min(omission), max(omission)] if omission else None,
            'admissible_deletions': len(omission),
            'without_nominations_difference_pp': float(differences[1]), 'sign_weights_difference_pp': float(differences[2])}


def rome_nominations(genes, data, cases, blocks, candidates):
    members = {r['id']: (r['entrez_ids'] if r['family'] == 'pathway' else [int(r['id'])]) for r in candidates}
    measured = np.column_stack([align(genes, data, m).mean(axis=1) for m in members.values()])
    pathways = [i for i,r in enumerate(candidates) if r['family'] == 'pathway']
    measured = np.column_stack([measured, measured[:,pathways].mean(axis=1)])
    mixed = [b for b in np.unique(blocks) if len(np.unique(cases[blocks == b])) == 2]
    keep = np.isin(blocks, mixed)
    obs, perms, _ = assignment_design(cases[keep], blocks[keep])
    conditional = obs @ measured[keep]
    unrestricted = measured[cases].mean(axis=0) - measured[~cases].mean(axis=0)
    tail = np.sum(np.abs(perms @ measured[keep]) >= np.abs(conditional)-1e-10, axis=0)
    return {k:{'unrestricted_pp':float(unrestricted[i]), 'matched_scanner_date_pp':float(conditional[i]),
               'exploratory_two_sided_extreme':int(tail[i]), 'assignments':len(perms)}
            for i,k in enumerate([*members, 'FIVE_PATHWAY_MEAN'])}


def analyze(data):
    reference, weight, candidates = weights_and_candidates()
    excluded = set()
    for r in candidates:
        excluded.update(r['entrez_ids'] if r['family'] == 'pathway' else [int(r['id'])])
    rome = data / 'gse25628'
    genes, scores, cases, unblocked = load_rome(rome)
    audited_cases, blocks, receipt = scan_design(rome)
    if not np.array_equal(cases, audited_cases):
        raise ValueError('Eligible scan and expression matrix groups disagree')
    primary_rome = transfer(genes, scores, cases, blocks, reference, weight, excluded)
    original_rome = transfer(genes, scores, cases, unblocked, reference, weight, excluded)
    nom = rome_nominations(genes, scores, cases, blocks, candidates)
    controls = hybridization_controls(rome / "GSE25628_series_matrix.txt.gz", cases, blocks)
    with LOCK_FILE.with_name('validation.tsv').open(newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            if abs(nom[row['id']]['unrestricted_pp']-float(row['rome_effect_pp'])) > 1e-9:
                raise ValueError('Rome frozen nomination sentinel failed')
    genes, scores, cases, _ = load_mumbai(data/'gse153739', data/'gse266265'/HGNC_FILE)
    mumbai = transfer(genes, scores, cases, np.zeros(len(cases), dtype=int), reference, weight, excluded)
    for report, adj in zip([primary_rome, mumbai], holm([primary_rome['p_two_sided'], mumbai['p_two_sided']]), strict=True):
        report['p_holm_two_cohorts'] = float(adj)
    return {'discovery_weights_sha256': DISCOVERY_SHA, 'original_weight_count': len(reference),
            'prior_nominated_gene_union_size':len(excluded), 'rome_scan_metadata':receipt,
            'rome_within_mixed_scan_dates':primary_rome, 'mumbai':mumbai,
            'rome_unrestricted_sensitivity':original_rome, 'rome_published_feature_diagnostics':nom,
            'rome_posthoc_hybridization_control_diagnostic':controls}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path('data/private'))
    args = parser.parse_args()
    print(json.dumps(analyze(args.data), indent=2, allow_nan=False))
