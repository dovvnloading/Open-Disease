"""Genome-wide transport of the frozen Finnish expression contrasts into Rome.

The primary statistic is a mean product of signed, randomization-standardized
contrasts. Enumerating Rome's possible donor assignments preserves all
cross-gene dependence. No participant-level output is generated.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from cross_cohort_discovery import align, assignment_design
from genome_rome import load_rome

SCREEN_SHA = '9df65b29615ec56767e9abaaf3844f8ebde5e21baf8cdab193ac24857a170da5'
LOCK_SHA = 'ba140f2bb026b0416e1c705ab28f616f430cc4295f94d174feaf5d15947c1e31'
TOL = 1e-10


def pinned_text(path, expected):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f'Input digest changed for {path.name}')
    return data.decode('utf-8')


def score_transport(profiles, reference, cases, design):
    """Compare an external reference with every synchronized assignment.

    At fixed scales, projecting each donor once is algebraically identical to
    first forming each assignment's entire contrast vector and taking its
    dot product with the Finnish reference.
    """
    observed, weights, strata_info = design
    if len(strata_info) != 1:
        raise ValueError('Genome transport requires a single unstratified target cohort')
    profiles = np.asarray(profiles, dtype=float)
    reference = np.asarray(reference, dtype=float)
    cases = np.asarray(cases)
    if (profiles.ndim != 2 or profiles.shape != (len(cases), len(reference))
            or reference.ndim != 1 or not len(reference) or cases.dtype != np.bool_
            or not np.isfinite(profiles).all() or not np.isfinite(reference).all()):
        raise ValueError('Expected finite aligned profiles, reference effects and binary labels')
    if min(cases.sum(), (~cases).sum()) < 2:
        raise ValueError('Need at least two donors in both groups for omissions')
    centered = profiles - profiles.mean(axis=0)
    scale = centered.std(axis=0, ddof=1) * np.sqrt(1 / cases.sum() + 1 / (~cases).sum())
    if np.any(scale <= 0):
        raise ValueError('Cannot standardize a constant gene')
    standardized = centered / scale
    donor_projection = standardized @ reference / len(reference)
    contrast = float(observed @ donor_projection)
    null = weights @ donor_projection
    n = len(null)
    extreme = int(np.count_nonzero(np.abs(null) >= abs(contrast) - TOL))
    omissions = []
    for i in range(len(cases)):
        keep = np.arange(len(cases)) != i
        omissions.append(donor_projection[keep & cases].mean() - donor_projection[keep & ~cases].mean())
    pairs = donor_projection[cases, None] - donor_projection[~cases]
    sign_oriented_difference_pp = float((observed @ centered) @ np.sign(reference) / len(reference))
    result = dict(genes=len(reference), mean_product=contrast, assignments=n,
                  sign_oriented_difference_pp=sign_oriented_difference_pp,
                  case_control_pairs=int(pairs.size),
                  pairs_with_higher_case_score=int((pairs > TOL).sum()),
                  pairs_tied=int((np.abs(pairs) <= TOL).sum()),
                  absolute_extreme=extreme, exact_two_sided_p=extreme / n,
                  lower_extreme=int(np.sum(null <= contrast + TOL)),
                  upper_extreme=int(np.sum(null >= contrast - TOL)),
                  null_q025=float(np.quantile(null, .025)), null_q975=float(np.quantile(null, .975)),
                  omission_min=float(min(omissions)), omission_max=float(max(omissions)))
    return result, observed @ standardized


def overlap(ids, reference, target):
    rows = []
    rank_ref = np.lexsort((ids, -np.abs(reference)))
    rank_target = np.lexsort((ids, -np.abs(target)))
    for pct in (1, 5, 10):
        size = int(np.ceil(len(ids) * pct / 100))
        selected = rank_ref[:size]
        both = np.intersect1d(selected, rank_target[:size])
        rows.append(dict(percent=pct, n=size, concordant_finnish_subset=int(np.sum(
            reference[selected] * target[selected] > 0)), shared_in_top=int(len(both)),
            shared_same_direction=int(np.sum(reference[both] * target[both] > 0)),
            shared_opposite_direction=int(np.sum(reference[both] * target[both] < 0))))
    return rows


def analyze(evidence, rome):
    screen = list(csv.DictReader(pinned_text(evidence / 'discovery_gene.tsv', SCREEN_SHA).splitlines(), delimiter='\t'))
    ids = np.array([int(r['id']) for r in screen])
    if len(ids) != 12044 or not (np.diff(ids) > 0).all():
        raise ValueError('Shared gene universe changed')
    reference = np.array([float(r['statistic']) for r in screen])
    lock = json.loads(pinned_text(evidence / 'nominations.json', LOCK_SHA))['candidates']
    rome_ids, profiles, cases, strata = load_rome(rome)
    if len(cases) != 14 or cases.sum() != 8 or len(np.unique(strata)) != 1:
        raise ValueError('Rome design changed')
    design = assignment_design(cases, strata)
    observed, weights, _ = design
    profiles = align(rome_ids, profiles, ids)
    with (evidence / 'validation.tsv').open() as stream:
        earlier = {(r['family'], r['id']): r for r in csv.DictReader(stream, delimiter='\t')}
    prior_members = set()
    max_difference = 0.
    for row in lock:
        members = row['entrez_ids'] if row['family'] == 'pathway' else [int(row['id'])]
        prior_members.update(members)
        donor_score = align(ids, profiles, members).mean(axis=1)
        diff = float(observed @ donor_score)
        sign = 1 if row['direction'] == 'higher' else -1
        count = int(np.sum((weights @ donor_score - diff) * sign >= -TOL))
        known = earlier[(row['family'], row['id'])]
        error = abs(diff - float(known['rome_effect_pp']))
        if error > TOL or count != int(known['rome_extreme']):
            raise ValueError('Previously verified Rome effects or counts do not reproduce')
        max_difference = max(max_difference, error)
    primary, z = score_transport(profiles, reference, cases, design)
    products = reference * z
    direction = dict(same=int((products > TOL).sum()), opposite=int((products < -TOL).sum()),
                     tied=int((np.abs(products) <= TOL).sum()), tolerance=TOL)
    remove_top = np.lexsort((ids, -np.abs(reference)))[:int(np.ceil(len(ids) / 100))]
    sensitivity = {}
    for name, mask in [('without_prior_nominations_and_members', ~np.isin(ids, list(prior_members))),
                       ('without_finnish_top_one_percent', ~np.isin(np.arange(len(ids)), remove_top))]:
        sensitivity[name], _ = score_transport(profiles[:, mask], reference[mask], cases, design)
    return dict(primary=primary, signs=direction, pearson=float(np.corrcoef(reference, z)[0, 1]),
                spearman=float(spearmanr(reference, z).statistic), top_overlap=overlap(ids, reference, z),
                sensitivity=sensitivity, sentinel_count=len(lock), sentinel_max_effect_difference=max_difference)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', type=Path, required=True)
    parser.add_argument('--rome', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.write_text(json.dumps(analyze(args.evidence, args.rome), indent=2, allow_nan=False) + '\n')
