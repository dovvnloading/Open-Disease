"""Check the main effects by explicit gene profiles and label enumeration.

Rome uses an independent native-header, matrix, and annotation parser. Mumbai
reuses the existing pinned transcript mapper and verifies downstream scoring.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from datetime import datetime
import gzip
from itertools import combinations, product
import json
from pathlib import Path
import re
import zlib

import numpy as np

from bulk_gse153739 import load_mumbai
from swiss_epithelial_mapping import HGNC_FILE, LOCK_FILE


def rome_reference(folder):
    with gzip.open(folder/'GPL571.annot.gz', 'rt', newline='') as fh:
        for line in fh:
            if line.strip() == '!platform_table_begin':
                break
        platform = {}
        for a in csv.DictReader(fh, delimiter='\t'):
            if a['ID'] == '!platform_table_end':
                break
            if (not a['ID'].startswith('AFFX') and re.fullmatch(r'[1-9][0-9]*', a['Gene ID'])
                    and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', a['Gene symbol'])):
                platform[a['ID']] = int(a['Gene ID'])
    with gzip.open(folder/'GSE25628_series_matrix.txt.gz','rt',newline='') as fh:
        reader = csv.reader(fh, delimiter='\t')
        clinical = []
        for row in reader:
            if row == ['!series_matrix_table_begin']:
                break
            if row and row[0] == '!Sample_characteristics_ch1':
                clinical.append(row[1:])
        accessions = next(reader)[1:]
        group, columns, stamps = [], [], []
        for i, values in enumerate(zip(*clinical, strict=True)):
            if 'disease state: Pathological (ectopic)' in values:
                continue
            if not ({'disease state: Normal (eutopic)', 'disease state: Normal (control)'} & set(values)):
                raise ValueError('Unexpected phenotype')
            columns.append(i + 1)
            group.append('disease state: Normal (eutopic)' in values)
            prefix = (folder/'ranges'/(accessions[i]+'.CEL.gz.part')).read_bytes()
            text = zlib.decompressobj(31).decompress(prefix,8192)
            stamp, = re.findall(rb'DatHeader=[^\n]+?([0-9]{2}/[0-9]{2}/[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}) ',text)
            stamps.append(datetime.strptime(stamp.decode(), '%m/%d/%y %H:%M:%S'))
        values, keys = [], []
        for row in reader:
            if row == ['!series_matrix_table_end']:
                break
            keys.append(platform.get(row[0], 0))
            values.append([float(row[k]) for k in columns])
    values = np.array(values).T
    if values.shape != (14,22277) or not np.isfinite(values).all():
        raise ValueError('Unexpected Rome matrix')
    by_gene = defaultdict(list)
    for i, key in enumerate(keys):
        if key:
            by_gene[key].append(i)
    ranks = []
    for donor in values:
        sorted_values = np.sort(donor)
        midrank = (np.searchsorted(sorted_values,donor,'left')+np.searchsorted(sorted_values,donor,'right'))*50.0/len(donor)
        ranks.append([midrank[ix].mean() for ix in by_gene.values()])
    return list(by_gene), np.array(ranks), np.array(group), stamps


def reverse_projection(genes, values, case, blocks, weights):
    keep = [(i, weights[int(g)]) for i,g in enumerate(genes) if int(g) in weights]
    idx, w = zip(*keep, strict=True)
    w = np.array(w)
    profile = values[:,idx]
    contributions = []
    observed = []
    for block in set(blocks):
        rows = np.array([i for i,b in enumerate(blocks) if b == block])
        ncase = int(case[rows].sum())
        if ncase in (0,len(rows)):
            continue
        # Average across *donors* first. Project every group gene-contrast vector
        # second. This reverses the production script's order of operations.
        def effect(selected):
            controls = np.setdiff1d(rows,selected)
            diff = profile[selected].mean(axis=0) - profile[controls].mean(axis=0)
            return float(np.inner(w,diff)/sum(abs(w)))
        observed.append(effect(rows[case[rows]]))
        contributions.append([effect(list(ix)) for ix in combinations(rows,ncase)])
    null = np.array([np.mean(x) for x in product(*contributions)])
    delta = float(np.mean(observed))
    return delta, int((abs(null) >= abs(delta)-1e-10).sum()), len(null)


def verify(data, aggregate):
    with LOCK_FILE.with_name('discovery_gene.tsv').open() as stream:
        weights = {int(r['id']):float(r['statistic']) for r in csv.DictReader(stream,delimiter='\t')}
    actual = json.loads(aggregate.read_text())
    genes, values, case, stamps = rome_reference(data/'gse25628')
    days = [s.date() for s in stamps]
    gaps = []
    for day in sorted(set(days)):
        cs = [s for i,s in enumerate(stamps) if s.date()==day and case[i]]
        ct = [s for i,s in enumerate(stamps) if s.date()==day and not case[i]]
        if cs and ct:
            gaps.append(int((min(cs)-max(ct)).total_seconds()))
    if gaps != actual['rome_scan_metadata']['earliest_case_minus_latest_control_seconds_by_mixed_date']:
        raise ValueError('Independent order audit disagrees')
    runs = {'rome_within_mixed_scan_dates':reverse_projection(genes, values, case, days, weights)}
    genes, values, case, _ = load_mumbai(data/'gse153739',data/'gse266265'/HGNC_FILE)
    runs['mumbai'] = reverse_projection(genes,values,case, [0]*len(case), weights)
    for key,(val,ext,n) in runs.items():
        record = actual[key]
        if abs(val-record['score_difference_pp']) > 1e-10 or (ext,n) != (record['exact_two_sided_extreme'],record['assignments']):
            raise ValueError('Independent score or tail disagrees')
    return {'score_effect_and_exact_tail_verified':list(runs), 'rome_native_time_gaps_verified_seconds':gaps,
            'absolute_score_tolerance_pp':1e-10}


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data',type=Path,default=Path('data/private'))
    ap.add_argument('--aggregate',type=Path,required=True)
    args = ap.parse_args()
    print(json.dumps(verify(args.data,args.aggregate),indent=2))
