"""Test the frozen Finnish score against reported Melbourne histologic stages."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from urllib.request import urlopen

import numpy as np
from scipy.stats import rankdata, t

from fetch_external_bulk_sources import verify_or_fetch
from genome_signature_transfer import weights_and_candidates
from open_disease.geo_microarray import sha256_file
from swiss_epithelial_mapping import HGNC_FILE, HGNC_SHA256, HGNC_SIZE, HGNC_URL

SERIES = 'GSE234354'
MATRIX = SERIES + '_gene_count_matrix.txt.gz'
MATRIX_SHA = '1f0406b1884aa12ba379b6ac4a8ada5994c85dea203af71bd4b8006ac066c669'
MATRIX_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234354/suppl/' + MATRIX
META_URL = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234354&targ=gsm&view=brief&form=text'
META_SHA = '73499ada2625d6b00c27854d0d2a82b8c77bed72d85c8f6fa21857913b69f937'
SEED = 234354
N_PERM = 49_999


def clinical(path):
    records = []
    for block in path.read_text(encoding='utf-8').split('^SAMPLE = ')[1:]:
        names = re.findall(r'^!Sample_title = endometrial biopsy, ((?:X|URS)\d+(?:_\d+)?)$', block, re.M)
        biosamples = re.findall(r'^!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/(SAMN\d+)$', block, re.M)
        runs = re.findall(r'^!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra\?term=(SRX\d+)$', block, re.M)
        fields = re.findall(r'^!Sample_characteristics_ch1 = ([^:]+): (.+)$', block, re.M)
        metadata = dict(fields)
        if (len(names) != 1 or len(biosamples) != 1 or len(runs) != 1 or len(fields) != 5
                or set(metadata) != {'tissue', 'age', 'cycle stage', 'pod', 'batch'}
                or metadata['tissue'] != 'Human endometrium'
                or metadata['batch'] not in ('batch 1','batch 2','batch 3')
                or metadata['cycle stage'] not in {'1','2','3','4','5','6','7','NA'}):
            raise ValueError('Malformed or unexpected Melbourne specimen metadata')
        records.append((names[0], metadata['cycle stage'], metadata['batch'], biosamples[0], runs[0]))
    if len(records) != 295 or any(len({r[i] for r in records}) != len(records) for i in (0,3,4)):
        raise ValueError('Unexpected library, BioSample or sequencing-experiment count')
    return sorted(records)


def clinical_digest(records):
    return hashlib.sha256(json.dumps(records, separators=(',', ':')).encode()).hexdigest()


def mapping(hgnc):
    if sha256_file(hgnc) != HGNC_SHA256:
        raise ValueError('Expected the October 2026 HGNC snapshot')
    ids, inverse = defaultdict(set), defaultdict(set)
    with hgnc.open(newline='') as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            stable, entrez = row['ensembl_gene_id'], row['entrez_id']
            if row['status'] == 'Approved' and re.fullmatch(r'ENSG\d{11}', stable) and re.fullmatch(r'[1-9]\d*', entrez):
                ids[stable].add(int(entrez))
                inverse[int(entrez)].add(stable)
    return {e:next(iter(v)) for e,v in ids.items() if len(v) == 1 and len(inverse[next(iter(v))]) == 1}


def read_profiles(matrix, hgnc, metadata):
    """Return phase-eligible biopsies by measured Finnish gene, in memory only."""
    if sha256_file(matrix) != MATRIX_SHA:
        raise ValueError('Melbourne RNA-seq count matrix differs from pinned source')
    valid = mapping(hgnc)
    base = [rec[0].split('_')[0] for rec in metadata]
    if Counter(Counter(base).values()) != {1:249, 2:23}:
        raise ValueError('Replicate structure changed')
    selected = {rec[0]:rec[1:3] for rec in metadata if rec[1] != 'NA'}
    if len(selected) != 257 or len({k.split('_')[0] for k in selected}) != len(selected):
        raise ValueError('More than one phase-eligible measurement per original biopsy')
    aligned, values, seen = [], [], set()
    with gzip.open(matrix,'rt',newline='') as stream:
        rows = csv.reader(stream,delimiter='\t')
        header = next(rows)
        if header[0] != 'gene_id' or len(header) != 296 or set(header[1:]) != {r[0] for r in metadata}:
            raise ValueError('Count columns do not match the source biopsy keys')
        use = [i for i, s in enumerate(header) if i and s in selected]
        stages = np.array([int(selected[header[i]][0]) for i in use])
        batches = np.array([selected[header[i]][1] for i in use])
        for row in rows:
            if len(row) != 296 or row[0] in seen or not re.fullmatch(r'ENSG\d{11}', row[0]):
                raise ValueError('Invalid or duplicate gene count row')
            seen.add(row[0])
            if row[0] in valid:
                nums = np.array([row[i] for i in use], dtype=np.int64)
                if np.any(nums < 0):
                    raise ValueError('Negative RNA-seq counts')
                aligned.append(valid[row[0]])
                values.append(nums)
    if len(seen) != 58302 or len(set(aligned)) != len(aligned):
        raise ValueError('Unexpected source universe or duplicate mapped gene')
    ref, w, _ = weights_and_candidates()
    shared, ix, ir = np.intersect1d(aligned, ref,return_indices=True)
    # The ranking background contains all approval-mapped genes, including zeros.
    arr = np.asarray(values)
    rank = (rankdata(arr, axis=0,method='average')-0.5)*100/len(aligned)
    stats = {'source_genes':len(seen), 'pinned_mapped_genes':len(aligned),
             'matched_finnish_genes':len(shared), 'absolute_weight_coverage':float(abs(w[ir]).sum()/abs(w).sum())}
    return shared, rank[ix].T, stages, batches, stats


def read_scores(matrix, hgnc, metadata):
    genes, profiles, stages, batches, stats = read_profiles(matrix, hgnc, metadata)
    reference, weights, _ = weights_and_candidates()
    matched = weights[np.searchsorted(reference, genes)]
    score = (matched / np.sum(np.abs(matched))) @ profiles.T
    return score, stages, batches, stats


def infer(score, stages, batches, *, draws=N_PERM):
    batch_levels = np.unique(batches)
    groups = [np.flatnonzero(batches == b) for b in batch_levels]
    stage_levels = np.unique(stages)
    reduced = np.column_stack([batches == b for b in batch_levels]).astype(float)
    phase = np.column_stack([stages == k for k in stage_levels[1:]]).astype(float)
    # Residualize the phase indicators with respect to batch. Q spans exactly
    # the part of phase that adds information beyond the batch-only model.
    residual = phase - reduced @ np.linalg.lstsq(reduced,phase,rcond=None)[0]
    q, _ = np.linalg.qr(residual,mode='reduced')
    full_design = np.column_stack([reduced, phase])
    if np.linalg.matrix_rank(full_design) != full_design.shape[1]:
        raise ValueError('Batch and phase contrasts are aliased')
    ss = float(np.sum((score @ q)**2))
    within_ss = sum(np.sum((score[g]-score[g].mean())**2) for g in groups)
    rng = np.random.default_rng(SEED)
    extreme = 0
    for start in range(0,draws,256):
        n = min(256,draws-start)
        permutations = np.empty((n,len(score)))
        for g in groups:
            choices = np.argsort(rng.random((n,len(g))), axis=1)
            permutations[:,g] = score[g][choices]
        extreme += int(np.sum(np.sum((permutations @ q)**2,axis=1) >= ss-1e-10))
    eligible = np.isin(stages,[3,6])
    x = np.column_stack([reduced[eligible], stages[eligible] == 6])
    y = score[eligible]
    beta = np.linalg.lstsq(x,y,rcond=None)[0]
    inv = np.linalg.inv(x.T @ x)
    hat = np.sum((x @ inv)*x, axis=1)
    residual = (y-x @ beta)/(1-hat)
    se = np.sqrt(np.sum(((x @ inv)[:,-1]*residual)**2))
    df = len(y) - x.shape[1]
    delta = t.ppf(.975, df) * se
    by_batch = {}
    for name, ix in zip(batch_levels,groups,strict=True):
        cells = {}
        for s in stage_levels:
            vals = score[ix[stages[ix] == s]]
            cells[str(s)] = {'n':len(vals),'mean_score':float(vals.mean()) if len(vals) else None}
        by_batch[str(name)] = {'cells':cells,'stage6_minus_stage3':cells['6']['mean_score']-cells['3']['mean_score']}
    return {'phase_partial_r2':ss/within_ss,'random_draws':draws,'extreme_random_draws':extreme,
            'plus_one_permutation_p':(extreme+1)/(draws+1),'seed':SEED,'phase_levels':stage_levels.tolist(),
            'stage6_minus_stage3':float(beta[-1]),'contrast_hc3_95_interval':[float(beta[-1]-delta),float(beta[-1]+delta)],
            'contrast_n':len(y),'by_batch':by_batch}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=Path('data/private'))
    parser.add_argument('--fetch',action='store_true',help='Acquire missing public files after checksum validation')
    args = parser.parse_args()
    folder = args.data / 'gse234354'
    matrix = folder / MATRIX
    hgnc = args.data / 'gse266265' / HGNC_FILE
    verify_or_fetch(hgnc,{'url':HGNC_URL,'sha256':HGNC_SHA256,'bytes':HGNC_SIZE},args.fetch)
    verify_or_fetch(matrix,{'url':MATRIX_URL,'sha256':MATRIX_SHA,'bytes':14249644},args.fetch)
    meta = folder/'meta'
    if args.fetch and not meta.exists():
        with urlopen(META_URL,timeout=90) as source:
            payload = source.read(1_200_001)
        if len(payload) > 1_200_000:
            raise ValueError('Response exceeded source metadata bound')
        with TemporaryDirectory(dir=folder) as temp:
            staged = Path(temp) / 'meta'
            staged.write_bytes(payload)
            if clinical_digest(clinical(staged)) != META_SHA:
                raise ValueError('Source phenotype or assay links changed')
            staged.replace(meta)
    records = clinical(meta)
    if clinical_digest(records) != META_SHA:
        raise ValueError('GEO mapping differs from pinned design')
    score, stages, batches, stats = read_scores(matrix,hgnc,records)
    result = {'cohort':SERIES,'matrix_sha256':MATRIX_SHA,'metadata_mapping_sha256':META_SHA,
              'archived_libraries':len(records),'distinct_source_biopsies':272,'known_phase_biopsies':len(score),
              'unique_biopsies_without_reported_phase':15,'redundant_libraries':23,**stats,**infer(score,stages,batches)}
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__ == '__main__':
    main()
