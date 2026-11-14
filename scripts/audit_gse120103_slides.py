"""Audit slide and extraction-protocol overlap in the New Delhi microarray study.

Use --fetch-headers to retrieve only the first 8 KiB of each official GEO gzip
into an ignored local cache. The default requires an existing cache. Never
emit individual records, paths, raw identifiers, or per-array signal values.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import re
from urllib.request import Request, urlopen
import zlib

from open_disease.geo_microarray import sha256_file

MATRIX_HASH = "0ab4900ba460232b1c76ca25028b3a034dd50491d3a7e63d4bf2fa2c4afbc16e"
PREFIX_SIZE = 8192
GROUPS = ("1A", "2A", "1B", "2B")


def sample_headers(matrix):
    if sha256_file(matrix) != MATRIX_HASH:
        raise ValueError("Unexpected processed matrix digest")
    wanted = {"!Sample_source_name_ch1", "!Sample_supplementary_file", "!Sample_geo_accession"}
    fields = {}
    with gzip.open(matrix, "rt", newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if row == ["!series_matrix_table_begin"]:
                break
            if row and row[0] in wanted:
                if row[0] in fields or len(row) != 37:
                    raise ValueError("Invalid 36-array header")
                fields[row[0]] = row[1:]
    for source, url, accession in zip(fields["!Sample_source_name_ch1"],
                                       fields["!Sample_supplementary_file"],
                                       fields["!Sample_geo_accession"], strict=True):
        group = re.search(r"\(Group ([12][AB])\)$", source)
        if group is None or re.fullmatch(r"GSM[0-9]+", accession) is None:
            raise ValueError("Invalid GEO clinical code or accession")
        expected = (r"ftp://ftp\.ncbi\.nlm\.nih\.gov/geo/samples/GSM[0-9]+nnn/"
                    + accession + r"/suppl/" + accession + r"_[\w.-]+\.txt\.gz")
        if re.fullmatch(expected, url) is None:
            raise ValueError("Unexpected supplementary location")
        yield group[1], url.replace("ftp://", "https://", 1)


def original_slide(prefix):
    # The FEPARAMS record precedes the feature table in Agilent FE text files.
    # Parse only the first three lines of the decompressed prefix.
    if len(prefix) != PREFIX_SIZE:
        raise ValueError("Incomplete header byte range")
    head = zlib.decompressobj(31).decompress(prefix, 128 * 1024).split(b"\n", 3)
    if len(head) != 4:
        raise ValueError("Missing Agilent header records")
    _, names, values = list(csv.reader(io.StringIO(b"\n".join(head[:3]).decode()), delimiter="\t"))
    if names[0] != "FEPARAMS" or values[0] != "DATA" or len(names) != len(values):
        raise ValueError("Invalid Agilent FEPARAMS record")
    names, values = names[1:], values[1:]
    if len(set(names)) != len(names):
        raise ValueError("Duplicate Agilent header field")
    fields = dict(zip(names, values, strict=True))
    barcode = re.fullmatch(r"([0-9]{12})_[0-9]_[0-9]", fields.get("FeatureExtractor_Barcode", ""))
    protocol = fields.get("Protocol_Name", "").strip()
    if barcode is None or not protocol:
        raise ValueError("Missing or unexpected array barcode or extraction protocol")
    slide = barcode[1]
    # Four-array Agilent slides append a grid position to the shared slide ID.
    # Check the shared ID independently against both original scanner fields.
    name = fields.get("FeatureExtractor_ArrayName", "")
    scan = re.split(r"[/\\]", fields.get("FeatureExtractor_ScanFileName", ""))[-1]
    if scan != name + ".tif" or not re.fullmatch(r"(?:US[0-9]+|NA)_" + slide + r"_S[0-9]+", name):
        raise ValueError("Scanner filename and slide barcode disagree")
    return slide, barcode[0], protocol


def audit(matrix, cache, fetch=False):
    if fetch:
        cache = cache.resolve()
        if not cache.is_relative_to(Path(__file__).resolve().parents[1] / "data/private"):
            raise ValueError("Save headers in an ignored data/private directory")
        cache.mkdir(parents=True, exist_ok=True)
    slides, protocols = defaultdict(Counter), defaultdict(Counter)
    manifest = hashlib.sha256()
    arrays = set()
    for group, url in sample_headers(matrix):
        path = cache / (hashlib.sha256(url.encode()).hexdigest() + ".prefix")
        if fetch and not path.exists():
            req = Request(url, headers={"Range": f"bytes=0-{PREFIX_SIZE - 1}"})
            with urlopen(req, timeout=25) as response:
                if response.status != 206 or not re.fullmatch(
                    rf"bytes 0-{PREFIX_SIZE - 1}/[0-9]+", response.headers.get("Content-Range", "")
                ):
                    raise ValueError("GEO did not return the requested bounded range")
                block = response.read(PREFIX_SIZE + 1)
                original_slide(block)
            path.write_bytes(block)
        block = path.read_bytes()
        slide, barcode, protocol = original_slide(block)
        if barcode in arrays:
            raise ValueError("Duplicated array position on a slide")
        arrays.add(barcode)
        manifest.update(hashlib.sha256(block).digest())
        slides[slide][group] += 1
        protocols[protocol][group] += 1
    total = sum(slides.values(), Counter())
    if total != Counter({group: 9 for group in GROUPS}):
        raise ValueError("Incomplete slide audit")
    return {
        "group_order": GROUPS, "n_arrays": sum(total.values()), "unique_array_positions": len(arrays), "unique_slides": len(slides),
        "slides_per_group": [sum(group in groups for groups in slides.values()) for group in GROUPS],
        "slides_spanning_clinical_groups": sum(len(counts) > 1 for counts in slides.values()),
        "protocol_group_counts": {protocol: [counts[group] for group in GROUPS]
                                  for protocol, counts in sorted(protocols.items())},
        "ordered_header_digest": manifest.hexdigest(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--fetch-headers", action="store_true")
    args = parser.parse_args()
    print(json.dumps(audit(args.matrix, args.cache, args.fetch_headers), indent=2))
