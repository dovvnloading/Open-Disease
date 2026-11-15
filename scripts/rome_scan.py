"""Read native Affymetrix scan metadata for the 14 eligible Rome biopsies."""

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime
import gzip
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import struct
import zlib

from analyze_gse25628_fosb import checked_groups, MATRIX_HASH
from bulk_gse6364 import fetch_prefix


def samples(directory):
    path = directory / "GSE25628_series_matrix.txt.gz"
    if sha256(path.read_bytes()).hexdigest() != MATRIX_HASH:
        raise ValueError("Unverified series matrix")
    header = defaultdict(list)
    with gzip.open(path, "rt", newline="") as f:
        for row in csv.reader(f, delimiter="\t"):
            if row == ["!series_matrix_table_begin"]:
                break
            if row and row[0].startswith(("!Series_", "!Sample_")):
                header[row[0]].append(row[1:])
    ids, groups = checked_groups(header)
    (urls,) = header["!Sample_supplementary_file"]
    records = []
    for sample, group, url in zip(ids, groups, urls, strict=True):
        expected = f"ftp://ftp.ncbi.nlm.nih.gov/geo/samples/{sample[:-3]}nnn/{sample}/suppl/{sample}.CEL.gz"
        if url != expected:
            raise ValueError("Sample URL inconsistent with GEO accession")
        if group != "lesion":
            records.append(
                (directory / (sample + ".header.gzpart"), url.replace("ftp:", "https:", 1), group)
            )
    return records


def load_prefix(record, fetch):
    path, url, group = record
    if not path.exists() and fetch:
        path.write_bytes(fetch_prefix(url))
    return path.read_bytes(), group


def scanner_fields(prefix):
    if len(prefix) != 2048:
        raise ValueError("Incomplete compressed CEL prefix")
    raw = zlib.decompressobj(31).decompress(prefix, 8192)
    magic, ver, cols, rows, n, size = struct.unpack_from("<6i", raw)
    if (magic, ver, cols, rows, n) != (64, 4, 732, 732, 732**2) or not 0 < size <= len(raw) - 24:
        raise ValueError("Invalid GPL571 geometry, CEL format or header size")
    headers = [
        s for s in raw[24 : 24 + size].decode("latin1").splitlines() if s.startswith("DatHeader=")
    ]
    if len(headers) != 1 or "HG-U133A_2.1sq" not in headers[0]:
        raise ValueError("Unable to corroborate chip design from native header")
    stamp = re.search(r"\b(\d{2}/\d{2}/\d{2}) (\d{2}:\d{2}:\d{2})\s+(\d+)\s+M10\b", headers[0])
    if not stamp:
        raise ValueError("Native scan date or scanner identifier is missing")
    return datetime.strptime(stamp[1], "%m/%d/%y").date(), stamp[3]


def audit(directory, fetch=False):
    records = samples(directory)
    with ThreadPoolExecutor(max_workers=4) as pool:
        prefixes = list(pool.map(lambda record: load_prefix(record, fetch), records))
    cells = defaultdict(Counter)
    stamps = []
    scanners = set()
    digest = sha256()
    for compressed, group in prefixes:
        day, scanner = scanner_fields(compressed)
        cells[day][group] += 1
        stamps.append(day)
        scanners.add(scanner)
        digest.update(compressed)
    case_controls = sorted((cell["case"], cell["control"]) for cell in cells.values())
    return stamps, dict(
        arrays=len(records),
        native_scanners=len(scanners),
        scan_dates=len(cells),
        mixed_scan_dates=sum(a > 0 and b > 0 for a, b in case_controls),
        scan_day_case_control_counts=case_controls,
        conditional_assignments=math.prod(math.comb(a + b, a) for a, b in case_controls),
        compressed_bytes=sum(len(x) for x, _ in prefixes),
        ordered_sha256=digest.hexdigest(),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--directory", type=Path, default=Path("data/private/gse25628"))
    p.add_argument("--fetch-headers", action="store_true")
    args = p.parse_args()
    print(json.dumps(audit(args.directory, args.fetch_headers)[1], indent=2))
