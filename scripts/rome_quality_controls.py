"""Exact mixed-scan-date checks of three prespecified Affymetrix control measures."""
import argparse
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
import re

import numpy as np

from analyze_gse25628_fosb import checked_groups
from cross_cohort_discovery import assignment_design, holm
from rome_scan_dates import uterine_scan_dates


def technical_controls(matrix, audit):
    dates = np.array(uterine_scan_dates(audit, matrix))
    with gzip.open(matrix, 'rt', newline='') as source:
        reader = csv.reader(source, delimiter='\t')
        header = defaultdict(list)
        for row in reader:
            if row == ['!series_matrix_table_begin']:
                break
            if row and row[0].startswith(('!Sample_', '!Series_')):
                header[row[0]].append(row[1:])
        keys, groups = checked_groups(header)
        if next(reader) != ['ID_REF', *keys]:
            raise ValueError('Probe signals are misaligned with cohort labels')
        select = np.array(groups) != 'lesion'
        case = (np.array(groups) == 'case')[select]
        controls = {row[0]: np.array(row[1:], dtype=float)[select] for row in reader
                    if row and row[0].startswith('AFFX-')}
    if len(controls) != 62 or not np.isfinite(list(controls.values())).all():
        raise ValueError('Rome control probe set changed')
    gradients = [controls[k+'_3_at'] - controls[k+'_5_at'] for k in
                 ['AFFX-HSAC07/X00351','AFFX-HUMGAPDH/M33197']]
    scores = {'housekeeping_3prime_minus_5prime': np.mean(gradients, axis=0)}
    for name, targets in [('hybridization_control_probes', ['Ec-bioB','Ec-bioC','Ec-bioD','P1-cre']),
                          ('polya_control_probes', ['Bs-dap','Bs-lys','Bs-phe','Bs-thr'])]:
        target_means = []
        for target in targets:
            matches = [v for k, v in controls.items() if re.fullmatch(
                'AFFX-r2-' + re.escape(target) + '-(?:3|5|M)(?:_s)?_at', k)]
            expected = 2 if target in ('Ec-bioC','Ec-bioD','P1-cre') else 3
            if len(matches) != expected:
                raise ValueError('Expected control target probes are missing')
            target_means.append(np.mean(matches, axis=0))
        scores[name] = np.mean(target_means, axis=0)
    mixed = np.isin(dates, [d for d in np.unique(dates) if len(np.unique(case[dates == d])) == 2])
    contrast, weights, _ = assignment_design(case[mixed], dates[mixed])
    out = {}
    for name, y in scores.items():
        obs = float(contrast @ y[mixed])
        null = weights @ y[mixed]
        tail = int(np.sum(np.abs(null) >= abs(obs) - 1e-10))
        out[name] = dict(log2_difference=obs, absolute_extreme=tail, assignments=len(null), two_sided_p=tail / len(null))
    for row, adjusted in zip(out.values(), holm([r['two_sided_p'] for r in out.values()]), strict=True):
        row['holm_three_controls'] = float(adjusted)
    return out


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--matrix', type=Path, required=True)
    p.add_argument('--scan-audit', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    args.out.write_text(json.dumps(technical_controls(args.matrix, args.scan_audit), indent=2) + '\n')
