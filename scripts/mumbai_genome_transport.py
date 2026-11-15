"""Test both fixed whole-genome Finnish and Rome references in Mumbai biopsies."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from bulk_gse153739 import load_mumbai
from cross_cohort_discovery import align, assignment_design, holm
from finland_rome_transport import SCREEN_SHA, pinned_text, score_transport
from genome_rome import load_rome


def transport(profiles, cases, references):
    design = assignment_design(cases, np.zeros(len(cases), dtype=int))
    rows = {}
    for name, vector in references.items():
        result, z = score_transport(profiles, vector, cases, design)
        product = vector * z
        result.update(pearson=float(np.corrcoef(vector, z)[0, 1]),
                      spearman=float(spearmanr(vector, z).statistic),
                      same=int((product > 1e-10).sum()), opposite=int((product < -1e-10).sum()),
                      ties=int((np.abs(product) <= 1e-10).sum()))
        rows[name] = result
    adjusted = holm([rows[k]['exact_two_sided_p'] for k in references])
    for name, p in zip(references, adjusted, strict=True):
        rows[name]['holm_two_references'] = float(p)
    return rows


def analyze(evidence, rome, mumbai, hgnc):
    table = list(csv.DictReader(pinned_text(evidence / 'discovery_gene.tsv', SCREEN_SHA).splitlines(), delimiter='\t'))
    ids = np.array([int(r['id']) for r in table])
    finland = np.array([float(r['statistic']) for r in table])
    rome_ids, expr, case_r, _ = load_rome(rome)
    rome_expr = align(rome_ids, expr, ids)
    rome_diff = rome_expr[case_r].mean(axis=0) - rome_expr[~case_r].mean(axis=0)
    rome_z = rome_diff / (rome_expr.std(axis=0, ddof=1) * np.sqrt(1 / case_r.sum() + 1 / (~case_r).sum()))
    mumbai_ids, profiles, cases, alternative = load_mumbai(mumbai, hgnc)
    if len(cases) != 7 or cases.sum() != 4:
        raise ValueError('Mumbai design changed')
    shared, orig_ix, new_ix = np.intersect1d(ids, mumbai_ids, assume_unique=True, return_indices=True)
    is_variable = np.ptp(profiles[:, new_ix], axis=0) > 0
    genes = shared[is_variable]
    selected = new_ix[is_variable]
    finn = finland[orig_ix[is_variable]]
    rom = rome_z[orig_ix[is_variable]]
    refs = {'finland': finn, 'rome': rom}
    primary = transport(profiles[:, selected], cases, refs)
    alt_variable = np.ptp(alternative[:, selected], axis=0) > 0
    alt_refs = {k: v[alt_variable] for k,v in refs.items()}
    sensitivity = transport(alternative[:, selected[alt_variable]], cases, alt_refs)
    ids_text = 'entrez_id\n' + ''.join(f'{gene}\n' for gene in genes)
    result = dict(finnish_rome_universe=len(ids), mumbai_mapped=len(mumbai_ids), common_mapped=len(shared),
                  zero_variance_primary=int((~is_variable).sum()), primary_genes=len(genes),
                  finland_rome_pearson_on_common=float(np.corrcoef(finn, rom)[0, 1]),
                  measured_universe_sha256=hashlib.sha256(ids_text.encode()).hexdigest(),
                  cases=int(cases.sum()), controls=int((~cases).sum()),
                  excluded_from_alternative=int((~alt_variable).sum()), primary=primary, transcript_first=sensitivity)
    return result, ids_text


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ('--evidence', '--rome', '--mumbai', '--hgnc', '--out', '--genes-out'):
        parser.add_argument(flag, type=Path, required=True)
    args = parser.parse_args()
    results, ids_text = analyze(args.evidence, args.rome, args.mumbai, args.hgnc)
    args.out.write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')
    args.genes_out.write_text(ids_text)
