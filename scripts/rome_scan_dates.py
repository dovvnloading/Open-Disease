"""Reconstruct native GSE25628 scan dates in the existing uterine matrix order."""
import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime
import gzip
from hashlib import sha256
import json
from pathlib import Path
import re
import struct
import zlib

from analyze_gse25628_fosb import MATRIX_HASH, checked_groups
from audit_us_bulk_sources import linked_files
from open_disease.geo_microarray import sha256_file

METADATA_HASH = "3aa501da2596734fca026219008765fc01a5b5f5571ac2b7d956e5a43f19956d"
HEADERS_HASH = "1ced4a97768d99d45161c2eda6c1ba0f01706f7524a247ca47a59e296d874418"


def scan_fields(part):
    if len(part) != 2048:
        raise ValueError("Expected a complete compressed header prefix")
    raw = zlib.decompressobj(31).decompress(part, 16384)
    magic, version, cols, rows, cells, length = struct.unpack_from("<6i", raw)
    if (magic, version, cols, rows, cells) != (64, 4, 732, 732, 732**2) or not 0 < length <= len(raw)-24:
        raise ValueError("Unexpected GPL571 XDA geometry")
    headers = [line for line in raw[24:24+length].splitlines() if line.startswith(b"DatHeader=")]
    if len(headers) != 1 or b"HG-U133A_2.1sq" not in headers[0]:
        raise ValueError("Expected one U133A 2.0 scanner header")
    match = re.search(rb"\b(\d{2}/\d{2}/\d{2} \d{2}:\d{2}:\d{2})\s+(\d+)\s+M10\b", headers[0])
    if not match:
        raise ValueError("Unrecognized DAT acquisition time and scanner")
    return datetime.strptime(match[1].decode(), "%m/%d/%y %H:%M:%S"), match[2]


def scan_records(directory, matrix_path):
    directory = Path(directory)
    soft = directory / "GSE25628.soft"
    if sha256_file(soft) != METADATA_HASH or sha256_file(matrix_path) != MATRIX_HASH:
        raise ValueError("Rome inputs differ from pinned sources")
    header = defaultdict(list)
    with gzip.open(matrix_path, "rt", newline="") as source:
        for row in csv.reader(source, delimiter="\t"):
            if row == ["!series_matrix_table_begin"]:
                break
            if row and row[0].startswith(("!Sample_", "!Series_")):
                header[row[0]].append(row[1:])
    ids, groups = checked_groups(header)
    records = {}
    digest = sha256()
    for key, _ in linked_files(soft, "GSE25628", 22):
        part = (directory / "GSE25628" / (key + ".part")).read_bytes()
        digest.update(part)
        records[key] = scan_fields(part)
    if digest.hexdigest() != HEADERS_HASH or set(ids) != records.keys() or len(set(records.values())) != 22:
        raise ValueError("Unverified scan identity, ordering or range digest")
    return [(group, *records[key]) for key, group in zip(ids, groups, strict=True)]


def uterine_scan_dates(directory, matrix_path):
    return [stamp.date() for group, stamp, _ in scan_records(directory, matrix_path) if group != "lesion"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("matrix_path", type=Path)
    args = parser.parse_args()
    rows = scan_records(args.directory, args.matrix_path)
    counts = defaultdict(Counter)
    for group, stamp, _ in rows:
        counts[stamp.date()][group] += 1
    print(json.dumps({"scanners": len({scanner for _, _, scanner in rows}),
                      "days_in_chronological_order": [counts[d] for d in sorted(counts)]}, indent=2))
