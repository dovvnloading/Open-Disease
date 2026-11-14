"""Audit only project metadata for donor, library and cycle-phase compatibility.

Takes the HECA supplementary workbook, Oxford E-MTAB-14039 SDRF and Huang
GSE214411 family SOFT gzip. Prints aggregate design data; never loads assays.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
from itertools import combinations
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET

from openpyxl import load_workbook
from audit_heca_glandular import SHA256, SIZE, STATE_I, STATE_II, is_case, records
from open_disease.geo_microarray import sha256_file

OXFORD = ("https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/039/E-MTAB-14039/Files/E-MTAB-14039.sdrf.txt",
          61604, "cdcd25f5519b3ddac7b68a41140ae577b929a856c35732e598fc95fcb0bbf36b")
HUANG = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE214nnn/GSE214411/soft/GSE214411_family.soft.gz",
         3247, "9f167f0e9c87d0d86f0031989e1d074384b55fe353d8fe2db1cf0c929892bacf")

ENA_SAMPLE = (2798, "13e0ef0709fc79b0f5fa7ee0fceb7bb827080d98b01b470d355c3f2293017f13")
ENA_EXPERIMENT = (2820, "be73ace6da69b5cddaeeaa3b99a4b649396f62102050210b9519ce7035177636")


def pin(path, size, digest):
    if path.stat().st_size != size or sha256_file(path) != digest:
        raise ValueError("Unrecognized metadata source")


def atlas_design(workbook):
    pin(workbook, SIZE, SHA256)
    book = load_workbook(workbook, read_only=True, data_only=True)
    try:
        donors, _ = records(book['Table1_donors'], 'source_key', ['Dataset', 'Endometriosis_stageh'])
        cells, _ = records(book['Table3_cell_type_counts'], 'Dataset_donorID', [STATE_I, STATE_II])
        qc, _ = records(book['Table2_CellRanger'], 'Library ID', ['snRNAseq'])
    finally:
        book.close()
    counts = {row['Dataset_donorID']: row[STATE_I]+row[STATE_II] for row in cells}
    clinical = {'Oxford': {}, 'Huang': {}}
    for row in donors:
        if row['Dataset'] in ('Mareckova-cells', 'Mareckova-cells & nuclei'):
            study, prefix = 'Oxford', 'Mareckova_'
        elif row['Dataset'] == 'Huangb':
            study, prefix = 'Huang', 'Huang_'
        else:
            continue
        donor = row['source_key']
        phase = row['Binary Stagee']
        eligible = (phase == 'Proliferative' and row['Locationd'] == 'endometrium'
                    and row['Hormonal treatmentg'] == 'NA' and counts[prefix+donor] >= 20)
        clinical[study][donor] = (is_case(row), phase, eligible)
    return clinical, {row['Library ID'] for row in qc if row['snRNAseq'] == 'scRNAseq'}


def allocation_count(libraries, cases):
    donors = set.union(*libraries)
    observed = donors & cases
    return sum(all(len(set(assigned) & batch) == len(observed & batch) for batch in libraries)
               for assigned in combinations(donors, len(observed)))


def summarize_graph(libraries, cases):
    donors = set.union(*libraries)
    mixed = [batch for batch in libraries if batch & cases and batch - cases]
    if not mixed:
        raise ValueError("No mixed library")
    joined = set.union(*mixed)
    composition = Counter((len(batch & cases), len(batch-cases)) for batch in libraries)
    return {
        'eligible_donors': len(donors), 'eligible_donor_library_links': sum(map(len, libraries)),
        'eligible_libraries': len(libraries), 'libraries_by_eligible_group_sizes': [
            {'cases': k[0], 'controls': k[1], 'libraries': v} for k, v in sorted(composition.items())],
        'eligible_donors_in_multiple_libraries': sum(sum(d in lib for lib in libraries) > 1 for d in donors),
        'mixed_libraries': len(mixed), 'unique_mixed_memberships': len({frozenset(g) for g in mixed}),
        'mixed_unique_cases': len(joined & cases), 'mixed_unique_controls': len(joined - cases),
        'libraries_per_distinct_mixed_pair': sorted(Counter(map(frozenset, mixed)).values()),
        'allocations_preserving_mixed_library_case_totals': allocation_count(mixed, cases),
        'allocations_preserving_all_eligible_library_case_totals': allocation_count(libraries, cases),
        'glandular_cells_by_donor_and_library_available_in_these_tables': False,
    }


def oxford_manifest(path, atlas, sc_libraries):
    pin(path, *OXFORD[1:])
    capture_donors = defaultdict(set)
    with path.open(encoding='utf-8', newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row['Characteristics[10x technology]'] != 'scRNAseq':
                continue
            key = row['Characteristics[individual]']
            library = row['Characteristics[library_id]']
            status = {'normal': False, 'endometriosis': True}[row['Characteristics[disease]']]
            if (key in capture_donors[library] or atlas[key][:2] != (status, row['Characteristics[menstrual_cycle_stage]'])
                    or library != row['Factor Value[library_id]']):
                raise ValueError("Inconsistent Oxford metadata or duplicate donor-library incidence")
            capture_donors[library].add(key)
    if set(capture_donors) != sc_libraries or set.union(*capture_donors.values()) != set(atlas):
        raise ValueError("Unmatched Oxford library or donor key")
    eligible = {d for d, (_, _, use) in atlas.items() if use}
    cases = {d for d in eligible if atlas[d][0]}
    return summarize_graph([g & eligible for g in capture_donors.values() if g & eligible], cases)


def huang_manifest(path, atlas, ena_sample=None, ena_experiment=None):
    pin(path, *HUANG[1:])
    with gzip.open(path, 'rt', encoding='utf-8') as handle:
        text = handle.read()
    samples = {}
    for block in text.split('^SAMPLE = ')[1:]:
        gsm = block.splitlines()[0]
        if gsm in samples or not re.fullmatch(r'GSM[0-9]+', gsm):
            raise ValueError("Unexpected GEO sample key")
        properties = dict(re.findall(r'^!Sample_characteristics_ch1 = (.+?): (.*)$', block, re.MULTILINE))
        experiments = re.findall(r'^!Sample_relation = SRA: .*?\b(SRX\d+)\b.*$', block, re.MULTILINE)
        biosamples = re.findall(r'^!Sample_relation = BioSample: .*?\b(SAMN\d+)\b.*$', block, re.MULTILINE)
        if len(experiments) != 1 or len(biosamples) != 1:
            raise ValueError("Expected one SRA experiment per deposited sample")
        samples[gsm] = properties, experiments[0], biosamples[0]
    phases, matched = Counter(), Counter()
    included_experiments = set()
    divergent = []
    for gsm, (case, stage, eligible) in atlas.items():
        props, experiment, biosample = samples[gsm]
        consistent = stage.lower() == props['cycle phase'].lower()
        if case != {'endometriosis': True, 'Control': False}[props['disease']]:
            raise ValueError("Huang diagnosis mismatch")
        phases[consistent] += 1
        if not consistent:
            divergent.append((biosample, experiment, props['cycle phase']))
        if eligible:
            matched[case, consistent] += 1
            included_experiments.add(experiment)
    result = {'geo_samples': len(samples), 'atlas_donors_exact_geo_matches': len(atlas),
            'diagnoses_agree': len(atlas), 'coarse_phases_agree': phases[True],
            'candidate_cases': sum(matched[True, ok] for ok in (True, False)),
            'candidate_controls': sum(matched[False, ok] for ok in (True, False)),
            'proliferative_in_both_cases': matched[True, True], 'proliferative_in_both_controls': matched[False, True],
            'candidate_unique_sra_experiments': len(included_experiments),
            'shared_capture_identifier_in_SOFT_characteristics': False}


    if (ena_sample is None) != (ena_experiment is None):
        raise ValueError("Both ENA XML metadata records are required together")
    if ena_sample is not None:
        pin(ena_sample, *ENA_SAMPLE)
        pin(ena_experiment, *ENA_EXPERIMENT)
        if len(divergent) != 1:
            raise ValueError("Divergent phase record is not unique")
        biosample, experiment, recorded_phase = divergent[0]
        sample = ET.parse(ena_sample).find('.//SAMPLE')
        exp = ET.parse(ena_experiment).find('.//EXPERIMENT')
        if (sample is None or biosample not in set(sample.itertext()) or exp is None
                or exp.get('accession') != experiment):
            raise ValueError("ENA records are not the sample and experiment linked by GEO")
        meta = {entry.findtext('TAG'): entry.findtext('VALUE') for entry in sample.iter('SAMPLE_ATTRIBUTE')}
        if meta.get('cycle phase') != recorded_phase:
            raise ValueError("GEO and ENA phase fields disagree")
        exp_text = ' '.join(exp.itertext()).lower()
        result['discordant_ena_metadata'] = {
            'biosample_phase_agrees_with_geo': True,
            'experiment_mentions_either_phase': any(phase in exp_text for phase in ('secretory','proliferative')),
            'biosample_bytes': ENA_SAMPLE[0], 'biosample_sha256': ENA_SAMPLE[1],
            'experiment_bytes': ENA_EXPERIMENT[0], 'experiment_sha256': ENA_EXPERIMENT[1],
        }
    return result


def audit(workbook, sdrf, soft, ena_sample=None, ena_experiment=None):
    clinical, qc = atlas_design(workbook)
    return {'sources': [{'url': url, 'bytes': size, 'sha256': h} for url, size, h in (OXFORD, HUANG)],
            'oxford': oxford_manifest(sdrf, clinical['Oxford'], qc),
            'huang': huang_manifest(soft, clinical['Huang'], ena_sample, ena_experiment)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workbook', type=Path)
    parser.add_argument('sdrf', type=Path)
    parser.add_argument('soft', type=Path)
    parser.add_argument('--ena-sample', type=Path, help='Optional ENA BioSample metadata XML linked by the discordant GSM')
    parser.add_argument('--ena-experiment', type=Path, help='Corresponding ENA SRA experiment XML')
    args = parser.parse_args()
    print(json.dumps(audit(args.workbook, args.sdrf, args.soft, args.ena_sample, args.ena_experiment), indent=2))
