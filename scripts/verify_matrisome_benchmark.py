"""Replay the benchmark with per-array sorted probe searches and physical cell identities."""
import argparse
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
import re
import zlib

import numpy as np

from bulk_gse153739 import load_mumbai
from verify_genome_transport import independent_rome
from verify_matrisome_cell_identity import cell_vector


def annotation(path, symbols):
    mapping, targets = {}, defaultdict(set)
    with gzip.open(path, 'rt') as stream:
        for line in stream:
            if line.strip() == '!platform_table_begin':
                break
        for row in csv.DictReader(stream, delimiter='\t'):
            probe = row['ID']
            if probe == '!platform_table_end':
                break
            if (not probe.startswith('AFFX') and re.fullmatch(r'[1-9][0-9]*', row['Gene ID'])
                    and re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', row['Gene symbol'])):
                mapping[probe] = int(row['Gene ID'])
                if row['Gene symbol'] in symbols:
                    targets[row['Gene symbol']].add(mapping[probe])
    genes = {g for matches in targets.values() if len(matches) == 1 for g in matches}
    return mapping, genes


def matrix_fields(path):
    fields = defaultdict(list)
    with gzip.open(path,'rt') as stream:
        for row in csv.reader(stream, delimiter='\t'):
            if row == ['!series_matrix_table_begin']:
                break
            if row and row[0].startswith('!Sample_'):
                fields[row[0]].append(row[1:])
    return fields


def labels(path, study):
    header = matrix_fields(path)
    groups = {}
    for i, key in enumerate(header['!Sample_geo_accession'][0]):
        chars = [r[i] for r in header['!Sample_characteristics_ch1']]
        if study == 'GSE51981':
            fields = dict(s.split(': ',1) for s in chars)
            state = fields['endometriosis/no endometriosis']
            if state == 'Non-Endometriosis' and fields['presence or absence of uterine/pelvic pathology'] != 'No Uterine Pelvic Pathology':
                continue
            case = state == 'Endometriosis'
            phase = fields['tissue'].removesuffix(' Endometrial tissue').lower().replace('-','_').replace(' ','_')
        else:
            phase, state = chars[0].split(' Phase ')
            case = {'Endometriosis':True, 'Normal':False}[state]
            phase = phase.lower().replace(' ','_')
        if phase in ('proliferative','early_secretory','mid_secretory'):
            groups[key] = (case, phase)
    return groups


def physical_id(directory, calvin, study, key):
    file = Path(study, key + '.part')
    short = (directory / file).read_bytes()
    raw = zlib.decompressobj(31).decompress(short)
    if raw.startswith(b'\x3b\x01'):
        full = (calvin / file).read_bytes()
        if not full.startswith(short):
            raise ValueError('Calvin source prefixes disagree')
        raw = zlib.decompressobj(31).decompress(full)
    return cell_vector(raw)[0]


def gene_matrix(path, selected, mapping, genes):
    with gzip.open(path,'rt') as stream:
        for line in stream:
            if line.strip() == '!series_matrix_table_begin':
                break
        rows = csv.reader(stream,delimiter='\t')
        header = next(rows)
        keys = sorted(selected)
        ix = [header.index(k) for k in keys]
        values, probe_genes = [], []
        for row in rows:
            if row == ['!series_matrix_table_end']:
                break
            values.append([float(row[i]) for i in ix])
            probe_genes.append(mapping.get(row[0],0))
    x = np.array(values,dtype=np.float64).T
    if not np.isfinite(x).all() or x.shape != (len(keys),54675):
        raise ValueError('Invalid complete GPL570 reference')
    probe_genes = np.array(probe_genes)
    keep = np.isin(probe_genes, genes)
    positions = np.searchsorted(genes,probe_genes[keep])
    denom = np.bincount(positions)
    design = np.empty((len(keys),len(genes)))
    for i,row in enumerate(x):
        ordered = np.sort(row)
        ranks = (np.searchsorted(ordered, row[keep],'left') + np.searchsorted(ordered, row[keep],'right')) * (50 / len(row))
        design[i] = np.bincount(positions,weights=ranks) / denom
    return keys, design


def check_pairs(values, cases, expected):
    positive = sum(a>b for a in values[cases] for b in values[~cases])
    tied = sum(a==b for a in values[cases] for b in values[~cases])
    if (positive, tied, len(values[cases]) * len(values[~cases])) != (expected['higher'],expected['ties'],expected['pairs']):
        raise ValueError('Independent pair ordering differs from benchmark')
    return int(positive)


def verify(args):
    expected = json.loads(args.expected.read_text())
    sources = ['GSE51981','GSE6364']
    with (args.scans / 'coldata.tsv').open() as stream:
        manifest = list(csv.DictReader(stream,delimiter='\t'))
    groups = [labels(args.directory / (s+'.gz'), s) for s in sources]
    for group, batch in zip(groups,['3','2'],strict=True):
        published = {r['sample_name']: r for r in manifest if r['batch']==batch}
        for k,(case,phase) in group.items():
            if k not in published or (case,phase) != (published[k]['condition']=='endometriosis',published[k]['phase']):
                raise ValueError('Primary GEO diagnosis or phase conflicts with classifier metadata')
    physical = [{k:physical_id(args.scans,args.calvin,s,k) for k in group} for s,group in zip(sources,groups,strict=True)]
    test_images = set(physical[1].values())
    rg,rx,rc = independent_rome(args.rome)
    mg,mx,mc,_ = load_mumbai(args.mumbai, args.hgnc)
    symbols = set(json.loads((args.directory/'all_matrisome_genes.json').read_text())['symbols'])
    mapping, mapped = annotation(args.directory/'GPL570.annot.gz',symbols)
    genes = np.array(sorted(mapped & set(rg) & set(mg)))
    if (len(symbols),len(mapped),genes.tolist()) != (expected['reference_symbols'],expected['GPL570_mapped_genes'],expected['common_entrez_ids']):
        raise ValueError('Matrisome projection disagrees')
    (ref, x),(test, y) = [gene_matrix(args.directory/(s+'.gz'),g,mapping,genes) for s,g in zip(sources,groups,strict=True)]
    report = {}
    for definition, result in expected['results'].items():
        use = np.array([definition=='with_reused_scans' or physical[0][k] not in test_images for k in ref])
        # Fit directly by giving each phase total mass +1/3 to cases and -1/3 to controls.
        weights = np.zeros(len(ref))
        phases = ('early_secretory','mid_secretory','proliferative')
        for phase in phases:
            for state,sign in [(True,1.),(False,-1.)]:
                ix = use & np.array([groups[0][k]==(state,phase) for k in ref])
                if sum(ix) != result['reference_counts_by_phase'][phase]['case' if state else 'control']:
                    raise ValueError('Reference eligibility disagrees')
                weights[ix] = sign/(3*sum(ix))
        ref_score = weights @ x  # A common positive normalization has no effect on pair order.
        y_score = y@ref_score
        raw = {}
        for phase in phases:
            ix=np.array([groups[1][k][1]==phase for k in test]); case=np.array([groups[1][k][0] for k in test])
            raw[phase]=check_pairs(y_score[ix],case[ix],result['GSE6364_by_phase'][phase])
        raw['Rome']=check_pairs(rx[:,np.searchsorted(rg,genes)]@ref_score,rc,result['Rome'])
        raw['Mumbai']=check_pairs(mx[:,np.searchsorted(mg,genes)]@ref_score,mc,result['Mumbai'])
        report[definition]=raw
    print(json.dumps(report,sort_keys=True))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for flag in ('directory','scans','calvin','rome','mumbai','hgnc','expected'):
        p.add_argument(flag,type=Path)
    verify(p.parse_args())
