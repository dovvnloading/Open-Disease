"""Rebuild the Australian score, regression and contrast interval independently."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from scipy.stats import t

from melbourne_signature_phase import MATRIX
from swiss_epithelial_mapping import HGNC_FILE, LOCK_FILE


def verify(data, aggregate):
    expected = json.loads(aggregate.read_text())
    folder = data / 'gse234354'
    counts = pd.read_csv(folder/MATRIX, sep='\t', index_col=0)
    annotation = pd.read_csv(data/'gse266265'/HGNC_FILE, sep='\t', dtype=str, keep_default_na=False)
    links = annotation.loc[annotation.status.eq('Approved') & annotation.ensembl_gene_id.str.fullmatch(r'ENSG\d{11}')
                           & annotation.entrez_id.str.fullmatch(r'[1-9]\d*'), ['ensembl_gene_id','entrez_id']].drop_duplicates()
    both = ~links.ensembl_gene_id.duplicated(keep=False) & ~links.entrez_id.duplicated(keep=False)
    links = links[both].set_index('ensembl_gene_id').entrez_id.astype(int)
    genes = counts.index.intersection(links.index)
    training = pd.read_csv(LOCK_FILE.with_name('discovery_gene.tsv'), sep='\t').set_index('id').statistic
    mapped = links.loc[genes]; common = mapped.isin(training.index)
    weights = training.loc[mapped[common]].to_numpy()
    scores, phases, batches = [], [], []
    for block in (folder/'meta').read_text().split('^SAMPLE = ')[1:]:
        stage = re.search(r'^!Sample_characteristics_ch1 = cycle stage: (\w+)$', block, re.M)[1]
        if stage == 'NA':
            continue
        key = re.search(r'^!Sample_title = endometrial biopsy, (\w+)$', block, re.M)[1]
        values = counts.loc[genes,key].to_numpy()
        ordered = np.sort(values)
        loc = np.flatnonzero(common.to_numpy())
        quantiles = (np.searchsorted(ordered,values[loc],'left')+np.searchsorted(ordered,values[loc],'right'))*50.0/len(ordered)
        scores.append(np.dot(weights,quantiles) / np.abs(weights).sum())
        phases.append(int(stage))
        batches.append(re.search(r'^!Sample_characteristics_ch1 = batch: ([^\n]+)$',block,re.M)[1])
    y, phase, batch = np.array(scores), np.array(phases), np.array(batches)
    batch_x = np.column_stack([batch==k for k in sorted(set(batch))])
    full_x = np.column_stack([batch_x, *[phase==k for k in range(2,8)]])
    def rss(x):
        err = y-x@np.linalg.lstsq(x,y,rcond=None)[0]
        return err@err
    r2 = 1-rss(full_x)/rss(batch_x)
    keep = np.isin(phase,[3,6]); contrast = (phase[keep]==6).astype(float)
    # Frisch-Waugh-Lovell coefficient plus direct sandwich covariance.
    residualized = contrast - batch_x[keep]@np.linalg.lstsq(batch_x[keep],contrast,rcond=None)[0]
    beta = (residualized @ y[keep]) / (residualized @ residualized)
    x = np.column_stack([batch_x[keep],contrast]); hat = x@np.linalg.pinv(x)
    err = (y[keep]-hat@y[keep])/(1-np.diag(hat)); inv = np.linalg.pinv(x)
    se = np.sqrt((inv @ np.diag(err*err) @ inv.T)[-1,-1]); delta = t.ppf(.975,len(err)-x.shape[1])*se
    if (len(y) != expected['known_phase_biopsies'] or common.sum() != expected['matched_finnish_genes']
            or not np.allclose([r2,beta,beta-delta,beta+delta], [expected['phase_partial_r2'],expected['stage6_minus_stage3'],
                               *expected['contrast_hc3_95_interval']],rtol=0,atol=1e-10)):
        raise ValueError('Independent rank, stage association or contrast disagrees')
    return {'recalculation':'Passed', 'phase_partial_r2':float(r2), 'stage6_minus_stage3':float(beta),
            'hc3_95_interval':[float(beta-delta),float(beta+delta)],'checked_biopsies':len(y),
            'matched_finnish_genes':int(common.sum()), 'absolute_tolerance':1e-10}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,default=Path('data/private'))
    parser.add_argument('--aggregate',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(verify(args.data,args.aggregate),indent=2))
