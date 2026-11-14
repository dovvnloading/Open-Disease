"""Synthetic checks for the frozen percentile, source joins and file guards."""
import csv
import gzip
import hashlib
import json

import pytest

from scripts import analyze_gse315857_fosb as study


LIBRARIES = ('OR9001', 'MOSE9002', 'N9003', 'OR9004', 'MOSE9005', 'OR9006', 'N9007', 'MOSE9008')
GROUPS = ('Control', 'Eosis', 'Control', 'Control', 'Eosis', 'Eosis', 'Control', 'Eosis')


def metadata(tmp_path, monkeypatch):
    blocks, assays, records = [], [], []
    for i, (lib, group) in enumerate(zip(LIBRARIES, GROUPS)):
        gsm, run, bio, experiment = [f'{prefix}99900{i}' for prefix in ('GSM', 'SRR', 'SAMN', 'SRX')]
        number = GROUPS[:i + 1].count(group)
        fields = {'!Sample_geo_accession': gsm, **study.GEO_FIXED,
                  '!Sample_title': f'{group} Prolif{number}', '!Sample_description': f'Library name: {lib}'}
        blocks.extend([f'^SAMPLE = {gsm}', *[f'{k} = {v}' for k, v in fields.items()],
                       f'!Sample_relation = SRA: https://www.ncbi.nlm.nih.gov/sra?term={experiment}',
                       f'!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/{bio}'])
        assays.append({'run_accession': run, 'sample_accession': bio, 'experiment_accession': experiment,
                       'library_name': gsm, 'study_accession': 'PRJNA1399448', 'instrument_platform': 'ILLUMINA',
                       'instrument_model': 'Illumina NovaSeq 6000', 'library_strategy': 'RNA-Seq',
                       'library_layout': 'PAIRED'})
        records.append((gsm, lib, group, run, bio, experiment))
    geo, ena = tmp_path / 'samples.soft', tmp_path / 'runs.tsv'
    geo.write_text('\n'.join(blocks) + '\n')
    with ena.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=assays[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(reversed(assays))
    monkeypatch.setattr(study, 'COHORT_HASH', study.digest_records(records))
    return geo, ena


def annotation(tmp_path, monkeypatch):
    ids = [f'ENST{9900+i:011d}.1' for i in range(12)] + [study.TARGET + '.8']
    rows = [{'name': key, 'geneName': 'FOSB', 'transcriptType': 'protein_coding',
             'tag': 'Ensembl_canonical,MANE_Select' if key == ids[-1] else 'basic'} for key in ids]
    data = {'genome': 'hg38', 'track': 'knownGene', 'chrom': 'chr19', 'start': 45450000,
            'end': 45550000, 'bigDataUrl': '/gbdb/hg38/gencode/gencodeV50.bb', 'knownGene': rows}
    records = [(r['name'], r['geneName'], r['transcriptType'], r['tag']) for r in rows]
    monkeypatch.setattr(study, 'ANNOTATION_HASH', study.digest_records(records))
    path = tmp_path / 'locus.json'
    path.write_text(json.dumps(data))
    return path


def matrix(tmp_path, monkeypatch, mutate=None):
    # Reversed header order relative to GEO. Three named rows. One tied value
    # in each control, low and high endpoints in cases. Unnamed rows never rank.
    target = [0.5, 0, 0.5, 0.5, 2.5, 2.5, 0.5, 2.5]
    rows = [['EnsembleID', *LIBRARIES[::-1]], ['ENST00000000001', *([0.5] * 8)],
            [study.TARGET, *target[::-1]], ['ENST00000000002', *([1.5] * 8)], ['', *([0] * 8)]]
    if mutate:
        mutate(rows)
    path = tmp_path / 'matrix.gz'
    with gzip.open(path, 'wt', newline='') as stream:
        csv.writer(stream, delimiter='\t').writerows(rows)
    for k, v in {'N_NAMED': 3, 'N_UNNAMED': 1, 'N_ZERO': 9, 'N_FRACTIONAL': 23,
                 'MATRIX_HASH': hashlib.sha256(path.read_bytes()).hexdigest()}.items():
        monkeypatch.setattr(study, k, v)
    return path


def test_complete_source_reordered_columns_midrank_and_leave_one_out(tmp_path, monkeypatch):
    geo, ena = metadata(tmp_path, monkeypatch)
    locus = annotation(tmp_path, monkeypatch)
    result = study.analyze(matrix(tmp_path, monkeypatch), geo, ena, locus)
    assert result['case_median_percentile'] == 100
    assert result['control_median_percentile'] == 25
    assert result['case_minus_control_percentile_points'] == 75
    assert result['leave_one_out_range_percentile_points'] == [75, 75]
    assert result['direction'] == 'descriptively_agrees'
    assert not any('SRR' in str(v) or 'GSM' in str(v) for v in result.values())


@pytest.mark.parametrize(('mod', 'error'), [
    (lambda rows: rows[0].__setitem__(1, rows[0][2]), 'columns'),
    (lambda rows: rows[1].__setitem__(0, study.TARGET), 'duplicate'),
    (lambda rows: rows[2].__setitem__(1, 'invalid-private-value'), 'Malformed'),
    (lambda rows: rows[2].__setitem__(1, 'nan'), 'finite'),
    (lambda rows: rows[3].__setitem__(1, '-1'), 'nonnegative'),
    (lambda rows: rows[1].pop(), 'Ragged'),
    (lambda rows: rows[-1].__setitem__(1, '1'), 'totals'),
    (lambda rows: rows[2].__setitem__(0, 'ENST00000000003'), 'canonical'),
])
def test_reject_corrupt_matrix(tmp_path, monkeypatch, mod, error):
    p = matrix(tmp_path, monkeypatch, mod)
    with pytest.raises(ValueError, match=error):
        study.transcript_percentiles(p, dict(zip(LIBRARIES, GROUPS)), {study.TARGET})


def test_reject_second_fosb_transcript_and_digest(tmp_path, monkeypatch):
    p = matrix(tmp_path, monkeypatch)
    groups = dict(zip(LIBRARIES, GROUPS))
    with pytest.raises(ValueError, match='canonical'):
        study.transcript_percentiles(p, groups, {study.TARGET, 'ENST00000000002'})
    monkeypatch.setattr(study, 'MATRIX_HASH', '0' * 64)
    with pytest.raises(ValueError, match='checksum'):
        study.transcript_percentiles(p, groups, {study.TARGET})


def test_source_rejects_swapped_experiments_and_ambiguous_labels(tmp_path, monkeypatch):
    geo, ena = metadata(tmp_path, monkeypatch)
    contents = ena.read_text().replace('SRX999000', 'temporary').replace('SRX999001', 'SRX999000')
    ena.write_text(contents.replace('temporary', 'SRX999001'))
    with pytest.raises(ValueError, match='join'):
        study.checked_cohort(geo, ena)
    geo.write_text(geo.read_text().replace('Control Prolif1', 'Eosis Prolif1'))
    with pytest.raises(ValueError, match='title'):
        study.checked_cohort(geo, ena)


def test_annotation_ignores_response_time_but_requires_release_and_tags(tmp_path, monkeypatch):
    path = annotation(tmp_path, monkeypatch)
    data = json.loads(path.read_text())
    data['downloadTime'] = 'a different timestamp'
    data['knownGene'].reverse()
    path.write_text(json.dumps(data))
    assert study.TARGET in study.checked_annotation(path)
    data['bigDataUrl'] = 'updated'
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='v50'):
        study.checked_annotation(path)
    data['bigDataUrl'] = '/gbdb/hg38/gencode/gencodeV50.bb'
    data['knownGene'][0]['tag'] = 'basic'
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='canonical'):
        study.checked_annotation(path)


def test_delete_one_range_can_cross_zero_and_keep_primary_direction():
    result = study.summarize(['Eosis'] * 4 + ['Control'] * 4, [0, 30, 60, 90, 20, 40, 45, 80])
    assert result['case_minus_control_percentile_points'] == 2.5
    assert result['leave_one_out_range_percentile_points'] == [-12.5, 17.5]
    assert study.summarize(GROUPS, [10]*8)['direction'] == 'no_observed_difference'
    assert study.summarize(GROUPS, [2 if g == 'Control' else 1 for g in GROUPS])['direction'] == 'descriptively_opposite'
