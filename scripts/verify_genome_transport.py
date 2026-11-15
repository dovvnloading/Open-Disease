"""Cross-check donor projections by directly constructing permuted gene effects.

Rome raw probe ranks and annotations are reconstructed independently. Mumbai
uses the existing study loader for transcript mapping and fixed-point FPKMs.
No arrays or individual scores are emitted.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from itertools import combinations
import json
from pathlib import Path
import re

import numpy as np

from bulk_gse153739 import load_mumbai


def table(path, begin, end):
    with gzip.open(path, 'rt', newline='') as stream:
        reader = csv.reader(stream, delimiter='\t')
        metadata = []
        for row in reader:
            if row == [begin]:
                break
            metadata.append(row)
        columns = next(reader)
        rows = []
        for row in reader:
            if row == [end]:
                return metadata, columns, rows
            rows.append(row)
    raise ValueError('Missing table boundary')


def independent_rome(rome):
    _, header, rows = table(rome / 'GPL571.annot.gz', '!platform_table_begin', '!platform_table_end')
    annotation = {}
    for row in rows:
        symbol, number, probe = [row[header.index(f)] for f in ('Gene symbol', 'Gene ID', 'ID')]
        if (re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', symbol)
                and re.fullmatch(r'[1-9][0-9]*', number) and not probe.startswith('AFFX')):
            annotation[probe] = int(number)
    meta, _, rows = table(rome / 'GSE25628_series_matrix.txt.gz', '!series_matrix_table_begin', '!series_matrix_table_end')
    clinical = [r[1:] for r in meta if r and r[0] == '!Sample_characteristics_ch1']
    disease, = [r for r in clinical if all(v.startswith('disease state: ') for v in r)]
    status, = [r for r in clinical if all(v.startswith('donor: ') for v in r)]
    selected = [i for i in range(len(status)) if (status[i], disease[i]) in {
        ('donor: patient','disease state: Normal (eutopic)'), ('donor: healthy donor','disease state: Normal (control)')}]
    cases = np.array([status[i] == 'donor: patient' for i in selected])
    values = np.array([[float(r[i+1]) for r in rows] for i in selected])
    if values.shape != (14,22277) or not np.isfinite(values).all() or cases.sum() != 8:
        raise ValueError('Rome reference data changed')
    rank = np.empty_like(values)
    for i, row in enumerate(values):
        ordered = np.sort(row)
        rank[i] = (np.searchsorted(ordered, row, side='left') + np.searchsorted(ordered, row, side='right')) * 50. / len(row)
    probe_ix = {}
    for i, row in enumerate(rows):
        gene = annotation.get(row[0])
        if gene is not None:
            probe_ix.setdefault(gene, []).append(i)
    genes = np.array(sorted(probe_ix))
    return genes, np.column_stack([rank[:, probe_ix[g]].mean(axis=1) for g in genes]), cases


def direct_dot(profile, cases, reference):
    """Derive each comparison vector and its SD from the enumeration itself."""
    allocations = list(combinations(range(len(cases)), int(cases.sum())))
    labels = np.zeros((len(allocations),len(cases)), dtype=bool)
    for i, included in enumerate(allocations):
        labels[i, list(included)] = True
    target = np.flatnonzero(np.all(labels == cases, axis=1)).item()
    weights = labels / cases.sum() - (~labels) / (~cases).sum()
    products = np.zeros(len(labels))
    observed_vector = []
    for start in range(0, len(reference), 128):
        stop = min(start+128, len(reference))
        contrasts = weights @ profile[:, start:stop]
        standardized = contrasts / contrasts.std(axis=0, ddof=0)
        products += np.sum(standardized * reference[start:stop], axis=1)
        observed_vector.extend(standardized[target])
    products /= len(reference)
    value = products[target]
    return float(value), int((np.abs(products) >= abs(value) - 1e-10).sum()), np.array(observed_vector)


def check(row, value, count):
    if not np.isclose(value, row['mean_product'], rtol=0, atol=1e-11) or count != row['absolute_extreme']:
        raise ValueError('Direct enumeration disagrees with the primary implementation')


def verify(args):
    with args.screen.open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    ids = np.array([int(r['id']) for r in rows])
    finland = np.array([float(r['statistic']) for r in rows])
    genes, x, cases = independent_rome(args.rome)
    index = np.searchsorted(genes, ids)
    if not np.array_equal(genes[index], ids):
        raise ValueError('Rome feature membership changed')
    value, count, rome_z = direct_dot(x[:,index], cases, finland)
    expected = json.loads(args.rome_results.read_text())
    check(expected['primary'], value, count)
    expected_mumbai = json.loads(args.mumbai_results.read_text())
    genes, x, cases, alternate = load_mumbai(args.mumbai, args.hgnc)
    ix_f = np.flatnonzero(np.isin(ids, genes))
    ix_m = np.searchsorted(genes, ids[ix_f])
    variable = np.ptp(x[:, ix_m], axis=0) != 0
    checked_ids = np.loadtxt(args.genes, skiprows=1, dtype=np.int64)
    if not np.array_equal(ids[ix_f][variable], checked_ids):
        raise ValueError('Mumbai gene universe disagrees with the public projection')
    ix_f, ix_m = ix_f[variable], ix_m[variable]
    evaluations = 1
    for label, expr in [('primary', x), ('transcript_first', alternate)]:
        keep = np.ptp(expr[:, ix_m], axis=0) != 0
        for name, ref in [('finland', finland), ('rome',rome_z)]:
            effect, extreme, _ = direct_dot(expr[:, ix_m[keep]], cases, ref[ix_f[keep]])
            check(expected_mumbai[label][name], effect, extreme)
            evaluations += 1
    print(f'{evaluations} exact tests agreed; independent Rome ranks and Mumbai membership verified.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ('--screen', '--rome', '--rome-results','--mumbai','--mumbai-results','--hgnc','--genes'):
        p.add_argument(flag, type=Path, required=True)
    verify(p.parse_args())
