"""Post-hoc technical check using the four external hybridization controls."""
from __future__ import annotations

from collections import defaultdict
import csv
import gzip
from pathlib import Path
import re

import numpy as np

from analyze_gse25628_fosb import MATRIX_HASH, checked_groups
from cross_cohort_discovery import assignment_design
from open_disease.geo_microarray import sha256_file

# Every matching probe for the four controls on the chip, including the r2 probes.
CONTROL = re.compile(r'AFFX-(?:r2-(?:Ec-|P1-))?(bioB|bioC|bioDn?|creX?)-(3|5|M)_at', re.I)


def hybridization_controls(matrix: Path, expected_cases, blocks):
    if sha256_file(matrix) != MATRIX_HASH:
        raise ValueError('Expected pinned GSE25628 matrix')
    header = defaultdict(list)
    grouped = defaultdict(list)
    with gzip.open(matrix, 'rt', newline='') as stream:
        reader = csv.reader(stream, delimiter='\t')
        for row in reader:
            if row == ['!series_matrix_table_begin']:
                break
            if row and row[0].startswith(('!Sample_', '!Series_')):
                header[row[0]].append(row[1:])
        ids, groups = checked_groups(header)
        retained = [i for i,g in enumerate(groups) if g != 'lesion']
        cases = np.array([groups[i] == 'case' for i in retained])
        if next(reader) != ['ID_REF',*ids] or not np.array_equal(cases, expected_cases):
            raise ValueError('Expected audited order and eligibility')
        for row in reader:
            if row == ['!series_matrix_table_end']:
                break
            match = CONTROL.fullmatch(row[0])
            if match:
                name = match[1].upper().removesuffix('N').removesuffix('X')
                grouped[name].append([float(row[i+1]) for i in retained])
    if {k:len(v) for k,v in grouped.items()} != {'BIOB':6, 'BIOC':4, 'BIOD':4, 'CRE':4}:
        raise ValueError('Unexpected external hybridization control coverage')
    mixed = [b for b in np.unique(blocks) if len(np.unique(cases[blocks == b])) == 2]
    keep = np.isin(blocks, mixed)
    obs, perms, _ = assignment_design(cases[keep], blocks[keep])
    names = sorted(grouped)
    by_target = np.column_stack([np.mean(grouped[name], axis=0) for name in names])
    scores = np.column_stack([by_target, by_target.mean(axis=1)])[keep]
    effects = obs @ scores
    tails = np.sum(np.abs(perms @ scores) >= np.abs(effects) - 1e-10, axis=0)
    return {name:{'difference_processed_log2':float(effect), 'probe_sets':sum(map(len, grouped.values())) if name == 'MEAN' else len(grouped[name]),
                  'date_only_reference_extreme':int(tail), 'assignments':len(perms)}
            for name, effect, tail in zip([*names, 'MEAN'], effects, tails, strict=True)}
