"""Check the genome-transfer study's pinned inputs. Pass --fetch for missing files."""
from __future__ import annotations

import argparse
from hashlib import sha256
from itertools import permutations
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen

from analyze_gse25628_fosb import ANNOTATION_HASH, MATRIX_HASH
from fetch_external_bulk_sources import EVIDENCE, verify_or_fetch
from rome_scan_design import scan_design
from swiss_epithelial_mapping import HGNC_FILE, HGNC_SIZE


def reconcile_ena(payload, expected):
    # ENA can return these seven exact rows in any order. Accept a reordering
    # only when the bytes, including all metadata fields, regain the pinned hash.
    rows = payload.splitlines(keepends=True)
    if len(rows) != 8 or len(set(rows)) != 8:
        raise ValueError('Expected header and seven distinct ENA rows')
    for body in permutations(rows[1:]):
        proposal = rows[0] + b''.join(body)
        if sha256(proposal).hexdigest() == expected:
            return proposal
    raise ValueError('ENA metadata has changed beyond a row rearrangement')


def check_sources(data, fetch):
    mumbai = json.loads((EVIDENCE / 'mumbai-source.json').read_text())
    for name, receipt in mumbai['sources'].items():
        path = data / 'gse153739' / name
        if name == 'ena-run-metadata.tsv' and fetch and not path.exists():
            with urlopen(receipt['url'], timeout=60) as response:
                payload = response.read(receipt['bytes'] + 1)
            if len(payload) != receipt['bytes']:
                raise ValueError('Unexpected ENA response length')
            restored = reconcile_ena(payload, receipt['sha256'])
            path.parent.mkdir(parents=True, exist_ok=True)
            with TemporaryDirectory(dir=path.parent) as temp:
                staged = Path(temp) / path.name
                staged.write_bytes(restored)
                staged.replace(path)
        verify_or_fetch(path, receipt, fetch)
    verify_or_fetch(data / 'gse266265' / HGNC_FILE, dict(mumbai['hgnc'], bytes=HGNC_SIZE), fetch)
    geo = 'https://ftp.ncbi.nlm.nih.gov/geo/'
    rome = data / 'gse25628'
    for suffix, length, digest in (
        ('series/GSE25nnn/GSE25628/matrix/GSE25628_series_matrix.txt.gz', 2795955, MATRIX_HASH),
        ('platforms/GPLnnn/GPL571/annot/GPL571.annot.gz', 4504226, ANNOTATION_HASH),
    ):
        verify_or_fetch(rome / suffix.rsplit('/', 1)[1], {'url':geo+suffix, 'sha256':digest, 'bytes':length}, fetch)
    return scan_design(rome, fetch=fetch)[2]


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data', type=Path, default=Path('data/private'))
    ap.add_argument('--fetch', action='store_true', help='Download and hash-check missing source files')
    args = ap.parse_args()
    print(json.dumps(check_sources(args.data, args.fetch), indent=2))
