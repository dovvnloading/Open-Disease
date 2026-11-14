"""Apply the prospective GSE315857 transcript-rank comparison to four local inputs."""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

TARGET = 'ENST00000353609'
MATRIX_HASH = 'fd83f667ab8df44b108ab499b2abbc0b16fd3b9b73fb2c899a410d908095df6d'
BED_HASH = 'f2709a357264408a1eeeda24b0dd8f56505b06520b51a09be678cbf5571e38d5'
COHORT_HASH = 'a07ce9637047dc20843ddca0b487297f9208aa0e133b3d5666ac5ea98f2b9027'
ANNOTATION_HASH = '4239eead87c0fc3d70885e19d04f3629e9712bec32a4d94b1f7b098e8757f72b'
N_NAMED = 22961
N_UNNAMED = 4
N_ZERO = 12581
N_FRACTIONAL = 171139
GEO_FIXED = {
    '!Sample_source_name_ch1': 'Eutopic endometrium',
    '!Sample_characteristics_ch1': 'tissue: Eutopic endometrium',
    '!Sample_platform_id': 'GPL24676', '!Sample_series_id': 'GSE315857',
    '!Sample_organism_ch1': 'Homo sapiens', '!Sample_taxid_ch1': '9606',
    '!Sample_library_strategy': 'RNA-Seq', '!Sample_instrument_model': 'Illumina NovaSeq 6000',
}
PREFIX_CELLS = {('Control', 'OR'): 2, ('Control', 'N'): 2,
                ('Eosis', 'OR'): 1, ('Eosis', 'MOSE'): 3}


def digest_records(records):
    payload = json.dumps(sorted(records), separators=(',', ':'), ensure_ascii=True)
    return hashlib.sha256(payload.encode('ascii')).hexdigest()


def checked_cohort(soft, ena):
    """Validate per-assay source and cross-archive mappings without disclosing them."""
    blocks = re.split(r'^\^SAMPLE = ', soft.read_text(encoding='utf-8'), flags=re.MULTILINE)
    if blocks[0].strip() or len(blocks) != 9:
        raise ValueError('Expected eight GEO sample records')
    geo = {}
    prefix_cells = Counter()
    titles = set()
    for block in blocks[1:]:
        lines = block.splitlines()
        fields = defaultdict(list)
        for line in lines[1:]:
            key, sep, value = line.partition(' = ')
            if sep:
                fields[key].append(value)
        gsm = lines[0]
        if (not re.fullmatch(r'GSM[1-9][0-9]*', gsm) or gsm in geo or
                fields['!Sample_geo_accession'] != [gsm] or
                any(fields[key] != [value] for key, value in GEO_FIXED.items())):
            raise ValueError('Unexpected sample identity, tissue, series or platform')
        labels = fields['!Sample_title']
        descriptions = fields['!Sample_description']
        if len(labels) != 1 or len(descriptions) != 1 or labels[0] in titles:
            raise ValueError('Ambiguous sample title or library')
        titles.add(labels[0])
        title = re.fullmatch(r'(Control|Eosis) Prolif[1-4]', labels[0])
        library = re.fullmatch(r'Library name: ((OR|MOSE|N)[0-9]+)', descriptions[0])
        if not title or not library:
            raise ValueError('Unrecognized sample phenotype or library name')
        group = title[1]
        prefix_cells[group, library[2]] += 1
        links = fields['!Sample_relation']
        srx = [m[1] for link in links if (m := re.fullmatch(
            r'SRA: https://www\.ncbi\.nlm\.nih\.gov/sra\?term=(SRX\d+)', link))]
        bio = [m[1] for link in links if (m := re.fullmatch(
            r'BioSample: https://www\.ncbi\.nlm\.nih\.gov/biosample/(SAMN\d+)', link))]
        if len(srx) != 1 or len(bio) != 1:
            raise ValueError('Missing or ambiguous SRA or BioSample relation')
        geo[gsm] = (library[1], group, bio[0], srx[0])
    if prefix_cells != PREFIX_CELLS or len({r[0] for r in geo.values()}) != 8:
        raise ValueError('Unexpected source group cells or duplicated library name')
    with ena.open(encoding='utf-8', newline='') as stream:
        runs = list(csv.DictReader(stream, delimiter='\t'))
    keys = ('library_name', 'run_accession', 'experiment_accession', 'sample_accession')
    if len(runs) != 8 or any(len({r[k] for r in runs}) != 8 for k in keys):
        raise ValueError('Expected eight distinct run, experiment and sample mappings')
    records = []
    for r in runs:
        gsm = r['library_name']
        if (gsm not in geo or r['study_accession'] != 'PRJNA1399448' or
                r['instrument_platform'] != 'ILLUMINA' or r['instrument_model'] != 'Illumina NovaSeq 6000' or
                r['library_layout'] != 'PAIRED' or r['library_strategy'] != 'RNA-Seq' or
                not re.fullmatch(r'SRR\d+', r['run_accession'])):
            raise ValueError('Unexpected ENA study, platform, assay or run')
        lib, group, bio, srx = geo[gsm]
        if (bio, srx) != (r['sample_accession'], r['experiment_accession']):
            raise ValueError('ENA to GEO sample/experiment join failed')
        records.append((gsm, lib, group, r['run_accession'], bio, srx))
    if digest_records(records) != COHORT_HASH:
        raise ValueError('Cohort mapping digest differs from the audited mapping')
    return {lib: group for lib, group, _, _ in geo.values()}


def checked_annotation(path):
    if path.suffix == '.bed':
        # bigBedToBed output from the permanent, versioned UCSC file at the
        # frozen locus. Pin the complete region before interpreting BED fields.
        with path.open('rb') as raw:
            if hashlib.file_digest(raw, 'sha256').hexdigest() != BED_HASH:
                raise ValueError('Versioned UCSC locus BED differs from the frozen source')
        with path.open(encoding='utf-8', newline='') as stream:
            bed = list(csv.reader(stream, delimiter='\t'))
        if len(bed) != 95 or any(len(r) != 27 for r in bed):
            raise ValueError('Unexpected UCSC v50 locus BED schema')
        rows = [{'name': r[3], 'geneName': r[17], 'transcriptType': r[22], 'tag': r[23]}
                for r in bed if r[17] == 'FOSB']
    else:
        data = json.loads(path.read_text(encoding='utf-8'))
        expected = {'genome': 'hg38', 'track': 'knownGene', 'chrom': 'chr19',
                    'start': 45450000, 'end': 45550000, 'bigDataUrl': '/gbdb/hg38/gencode/gencodeV50.bb'}
        if any(data.get(k) != v for k, v in expected.items()):
            raise ValueError('Expected the hg38 GENCODE v50 locus annotation')
        rows = [r for r in data['knownGene'] if r['geneName'] == 'FOSB']
    records = [(r['name'], r['geneName'], r['transcriptType'], ','.join(sorted(r['tag'].split(','))))
               for r in rows]
    targets = [r for r in rows if r['name'] == TARGET + '.8']
    if (len(rows) != 13 or len(targets) != 1 or targets[0]['transcriptType'] != 'protein_coding' or
            not {'MANE_Select', 'Ensembl_canonical'}.issubset(targets[0]['tag'].split(',')) or
            digest_records(records) != ANNOTATION_HASH):
        raise ValueError('FOSB annotation or canonical transcript differs from freeze')
    return {r['name'].split('.')[0] for r in rows}


def transcript_percentiles(matrix, groups, annotated):
    with matrix.open('rb') as raw:
        if hashlib.file_digest(raw, 'sha256').hexdigest() != MATRIX_HASH:
            raise ValueError('Matrix checksum differs from the frozen release')
    named = set()
    unnamed = zeros = fractional = 0
    target = None
    with gzip.open(matrix, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.reader(stream, delimiter='\t')
        header = next(reader, [])
        if (len(header) != 9 or header[0] != 'EnsembleID' or
                set(header[1:]) != set(groups) or len(set(header[1:])) != 8):
            raise ValueError('Matrix columns do not match the audited libraries')
        columns = [[] for _ in header[1:]]
        for row in reader:
            if len(row) != len(header):
                raise ValueError('Ragged transcript row')
            key = row[0]
            try:
                values = [float(x) for x in row[1:]]
            except ValueError:
                raise ValueError('Malformed matrix abundance') from None
            if any(not math.isfinite(x) or x < 0 for x in values):
                raise ValueError('Expected complete finite nonnegative abundances')
            zeros += sum(x == 0 for x in values)
            fractional += sum(not x.is_integer() for x in values)
            if not key:
                unnamed += 1
                continue
            if not re.fullmatch(r'ENST[0-9]{11}', key) or key in named:
                raise ValueError('Unexpected or duplicate transcript ID')
            named.add(key)
            for col, value in zip(columns, values):
                col.append(value)
            if key == TARGET:
                target = values
    if (len(named), unnamed, zeros, fractional) != (N_NAMED, N_UNNAMED, N_ZERO, N_FRACTIONAL):
        raise ValueError('Whole-matrix size or audited cell totals differ')
    if named & annotated != {TARGET} or target is None:
        raise ValueError('Expected exactly the canonical FOSB transcript in the matrix')
    scores = []
    for col, value in zip(columns, target):
        below = sum(x < value for x in col)
        tied = sum(x == value for x in col)
        scores.append(100 * (below + (tied - 1) / 2) / (len(col) - 1))
    return [groups[label] for label in header[1:]], scores


def summarize(groups, scores):
    def medians(skip=None):
        return [median(v for i, (g, v) in enumerate(zip(groups, scores)) if g == group and i != skip)
                for group in ('Eosis', 'Control')]
    affected, control = medians()
    difference = affected - control
    omitted = [a - c for a, c in (medians(i) for i in range(len(scores)))]
    direction = 'no_observed_difference'
    if difference > 0:
        direction = 'descriptively_agrees'
    elif difference < 0:
        direction = 'descriptively_opposite'
    return {'case_median_percentile': affected, 'control_median_percentile': control,
            'case_minus_control_percentile_points': difference, 'direction': direction,
            'leave_one_out_range_percentile_points': [min(omitted), max(omitted)]}


def analyze(matrix, geo, ena, annotation):
    groups = checked_cohort(geo, ena)
    transcripts = checked_annotation(annotation)
    labels, scores = transcript_percentiles(matrix, groups, transcripts)
    return {'series': 'GSE315857', 'transcript': TARGET, 'case_n': 4, 'control_n': 4,
            'reference_transcripts': N_NAMED, 'unnamed_rows_excluded': N_UNNAMED,
            **summarize(labels, scores)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('matrix', 'geo', 'ena', 'annotation'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    try:
        report = analyze(args.matrix, args.geo, args.ena, args.annotation)
    except (OSError, EOFError, KeyError, TypeError, ValueError, csv.Error) as exc:
        # Never display malformed input cells, keys, private paths or source records.
        parser.exit(1, f'Analysis blocked by a source or integrity check ({type(exc).__name__}).\n')
    print(json.dumps(report, indent=2, allow_nan=False))
