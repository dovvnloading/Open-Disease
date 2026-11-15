"""Recover scan-day balance from the native GPL571 CEL headers in GSE25628."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime
import gzip
import hashlib
import json
from math import comb, prod
from pathlib import Path
import re
import struct
from tempfile import TemporaryDirectory
from urllib.request import Request, urlopen
import zlib

import numpy as np

from analyze_gse25628_fosb import MATRIX_HASH, checked_groups
from open_disease.geo_microarray import sha256_file

PREFIX_BYTES = 4096
ORDERED_DIGEST = '437fe28852103bc6441560fbe6de85928bba244d06f4f4d66db4d044a701801d'


def scan_fields(prefix):
    if len(prefix) != PREFIX_BYTES:
        raise ValueError('Expected a 4 KiB compressed CEL prefix')
    raw = zlib.decompressobj(31).decompress(prefix, 8192)
    magic, version, cols, rows, cells, header_size = struct.unpack_from('<6i', raw)
    if (magic, version, cols, rows, cells) != (64, 4, 732, 732, 732**2) or not 0 < header_size <= len(raw)-24:
        raise ValueError('Unexpected GPL571 native CEL format')
    lines = raw[24:24 + header_size].decode('latin1').splitlines()
    native = [s for s in lines if s.startswith('DatHeader=') and 'HG-U133A_2.1sq' in s]
    if len(native) != 1:
        raise ValueError('Missing GPL571 DAT header')
    matched = re.search(r'\b(\d{2}/\d{2}/\d{2} \d{2}:\d{2}:\d{2})\s+(\d+)\s+M10\b', native[0])
    if not matched:
        raise ValueError('Missing acquisition date or scanner')
    return matched[2], datetime.strptime(matched[1], '%m/%d/%y %H:%M:%S')


def scan_design(directory, *, fetch=False):
    directory = Path(directory)
    matrix = directory / 'GSE25628_series_matrix.txt.gz'
    if sha256_file(matrix) != MATRIX_HASH:
        raise ValueError('Expected pinned Rome matrix')
    fields = defaultdict(list)
    with gzip.open(matrix, 'rt', newline='') as stream:
        for row in csv.reader(stream, delimiter='\t'):
            if row == ['!series_matrix_table_begin']:
                break
            if row and row[0].startswith(('!Sample_', '!Series_')):
                fields[row[0]].append(row[1:])
    ids, groups = checked_groups(fields)
    archives, = fields['!Sample_supplementary_file']
    summary = defaultdict(Counter)
    observed = set()
    eligible_times = defaultdict(lambda: defaultdict(list))
    blocks, cases = [], []
    digest = hashlib.sha256()
    folder = directory / 'ranges'
    for accession, group, archive in zip(ids, groups, archives, strict=True):
        url = archive.replace('ftp://', 'https://', 1)
        if not re.fullmatch(r'https://ftp\.ncbi\.nlm\.nih\.gov/geo/samples/GSM\d+nnn/' + accession + '/suppl/' + accession + r'\.CEL\.gz', url):
            raise ValueError('Unexpected GEO CEL URL')
        name = folder / (accession + '.CEL.gz.part')
        if fetch and not name.is_file():
            folder.mkdir(exist_ok=True)
            request = Request(url, headers={'Range': 'bytes=0-4095'})
            with urlopen(request, timeout=45) as response:
                if response.status != 206 or not re.fullmatch(r'bytes 0-4095/\d+', response.headers.get('Content-Range', '')):
                    raise ValueError('Upstream did not honor the bounded range request')
                part = response.read(PREFIX_BYTES + 1)
            scan_fields(part)
            with TemporaryDirectory(dir=folder) as temp:
                staged = Path(temp) / name.name
                staged.write_bytes(part)
                staged.replace(name)
        part = name.read_bytes()
        digest.update(part)
        scanner, timestamp = scan_fields(part)
        if (scanner, timestamp) in observed:
            raise ValueError('Repeated scanner timestamp')
        observed.add((scanner, timestamp))
        block = scanner + '/' + timestamp.date().isoformat()
        summary[block][group] += 1
        if group != 'lesion':
            eligible_times[block][group].append(timestamp)
            blocks.append(block)
            cases.append(group == 'case')
    if digest.hexdigest() != ORDERED_DIGEST:
        raise ValueError('Ordered native CEL header prefixes differ from audited files')
    counts = [c for _, c in sorted(summary.items())]
    time_gaps = []
    for block in sorted(eligible_times):
        times = eligible_times[block]
        if times['case'] and times['control']:
            time_gaps.append(int((min(times['case']) - max(times['control'])).total_seconds()))
    receipt = {'source': 'GSE25628', 'native_scanners': len({s for s, t in observed}),
               'scanner_dates': len(counts), 'source_arrays': len(observed),
               'earliest_case_minus_latest_control_seconds_by_mixed_date': time_gaps,
               'mixed_dates_all_controls_scanned_before_all_cases': sum(gap > 0 for gap in time_gaps),
               'mixed_scanner_dates': sum(c['case'] > 0 and c['control'] > 0 for c in counts),
               'by_date_chronological': [{k: c[k] for k in ('case', 'control', 'lesion')} for c in counts],
               'eligible_conditional_assignments': prod(comb(c['case']+c['control'], c['case']) for c in counts),
               'prefix_bytes_each': PREFIX_BYTES, 'ordered_prefix_sha256': digest.hexdigest()}
    return np.array(cases), np.array(blocks), receipt


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('directory', type=Path)
    ap.add_argument('--fetch', action='store_true', help='Download 4 KiB per missing GEO scan header')
    args = ap.parse_args()
    print(json.dumps(scan_design(args.directory, fetch=args.fetch)[2], indent=2))
