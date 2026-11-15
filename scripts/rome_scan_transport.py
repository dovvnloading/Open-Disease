"""Condition Finland-to-Rome transfer on native scanner date with fixed scales."""
from __future__ import annotations

import argparse
import csv
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np

from cross_cohort_discovery import align, assignment_design
from finland_rome_transport import SCREEN_SHA, TOL, pinned_text
from genome_rome import load_rome


def compare_mixed_dates(profiles, reference, cases, dates):
    """Compute one contrast per date where cases and controls share a scanner day."""
    dates = np.asarray(dates)
    if cases.shape != dates.shape or profiles.shape != (len(cases), len(reference)):
        raise ValueError('Profiles, case indicators and dates must align')
    mixed = [d for d in np.unique(dates) if np.unique(cases[dates == d]).size == 2]
    if not mixed:
        raise ValueError('No native scan date contains both groups')
    keep = np.isin(dates, mixed)
    centered = profiles - profiles.mean(axis=0)
    scale = np.std(centered, axis=0, ddof=1) * np.sqrt(1 / cases.sum() + 1 / (~cases).sum())
    if np.any(scale <= 0):
        raise ValueError('Constant gene cannot be standardized')
    full_standardized = centered / scale
    x = full_standardized[keep]
    case, dates = cases[keep], dates[keep]
    observed, weights, _ = assignment_design(case, dates)
    score = x @ reference / len(reference)
    joint = float(observed @ score)
    null = weights @ score
    tail = int(np.count_nonzero(abs(null) >= abs(joint) - TOL))
    signed_pp = float(observed @ centered[keep] @ np.sign(reference) / len(reference))
    groups, extreme_checks = [], []
    omitted = []
    for day in mixed:
        inside = dates == day
        affected = inside & case
        unaffected = inside & ~case
        pairs = score[affected, None] - score[unaffected]
        rec = dict(cases=int(sum(affected)), controls=int(sum(unaffected)),
                   mean_product=float(score[affected].mean() - score[unaffected].mean()),
                   case_control_pairs=int(pairs.size), pairs_case_higher=int(np.sum(pairs > TOL)),
                   pairs_tied=int(np.sum(abs(pairs) <= TOL)))
        groups.append(rec)
        # Verify with independent direct binary subset means within each date.
        one_day = score[inside]
        options = []
        for membership in combinations(range(len(one_day)), int(sum(affected))):
            mask = np.zeros(len(one_day), dtype=bool)
            mask[list(membership)] = True
            options.append(one_day[mask].mean() - one_day[~mask].mean())
        extreme_checks.append(options)
    direct = np.array([np.mean(ds) for ds in product(*extreme_checks)])
    if not np.allclose(np.sort(null), np.sort(direct), atol=TOL, rtol=0):
        raise ValueError('Independent blocked enumeration disagrees with design matrix')
    for omit in range(len(case)):
        include = np.arange(len(case)) != omit
        contrasts = []
        for d in mixed:
            here = include & (dates == d)
            if (here & case).any() and (here & ~case).any():
                contrasts.append(score[here & case].mean() - score[here & ~case].mean())
        if len(contrasts) == len(mixed):
            omitted.append(np.mean(contrasts))
    return dict(genes=len(reference), cases=int(sum(case)), controls=int(sum(~case)),
                mixed_dates=len(mixed), excluded_arrays=int(sum(~keep)),
                assignments=len(null), absolute_extreme=tail, two_sided_p=tail / len(null),
                lower_extreme=int(np.sum(null <= joint + TOL)), upper_extreme=int(np.sum(null >= joint - TOL)),
                mean_product=joint, pearson=float(np.corrcoef(observed @ x, reference)[0,1]),
                sign_oriented_difference_pp=signed_pp, case_control_pairs=sum(g['case_control_pairs'] for g in groups),
                pairs_case_higher=sum(g['pairs_case_higher'] for g in groups),
                pairs_tied=sum(g['pairs_tied'] for g in groups),
                omission_min=float(min(omitted)) if omitted else None,
                omission_max=float(max(omitted)) if omitted else None, date_specific=groups)


if __name__ == '__main__':
    from rome_scan_dates import scan_records
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ('--evidence', '--rome','--scan-audit','--out'):
        p.add_argument(flag, type=Path, required=True)
    args = p.parse_args()
    rows = list(csv.DictReader(pinned_text(args.evidence / 'discovery_gene.tsv', SCREEN_SHA).splitlines(), delimiter='\t'))
    ids = np.array([int(r['id']) for r in rows])
    reference = np.array([float(r['statistic']) for r in rows])
    rome_ids, expr, cases, _ = load_rome(args.rome)
    records = scan_records(args.scan_audit, args.rome / 'GSE25628_series_matrix.txt.gz')
    uterine = [(group, stamp) for group, stamp, _ in records if group != 'lesion']
    if not np.array_equal(cases, [group == 'case' for group, _ in uterine]):
        raise ValueError('Native chronology differs from observed cohort order')
    dates = [stamp.date() for _, stamp in uterine]
    result = compare_mixed_dates(align(rome_ids, expr, ids), reference, cases, dates)
    earlier = sum(a < b for g, a in uterine if g == 'case' for h, b in uterine
                  if h == 'control' and a.date() == b.date())
    result['same_date_case_scanned_before_control'] = earlier
    result['same_date_case_scanned_after_control'] = result['case_control_pairs'] - earlier
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
