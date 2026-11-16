"""Audit Affymetrix acquisition metadata without decoding intensity cells."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import struct
from urllib.parse import quote
from urllib.request import Request, urlopen

from audit_sohler_peritoneum import SDRF_SHA256, TISSUES, summarize
from sohler_conditional_design import conditional_design, verify_dat_filename
from genome_rome import read_annotation
from open_disease.geo_microarray import sha256_file


BASE = "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/694/E-MTAB-694/Files/"
STUDY_SHA = "d2586e9fecd9461a7f5aeefec898bae5d2526dd541532ad383896717230aef31"
ANNOTATION_SHA = "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394"
PREFIXES_SHA = "89d1ba2a3e645277a652624fa07c4c147f894a5f33448ad683e4bda3d22b6244"
HEADER_BYTES = 512
ROOT = Path(__file__).resolve().parents[1]


def study_files(node):
    if isinstance(node, list):
        for child in node:
            yield from study_files(child)
    elif isinstance(node, dict):
        if node.get("type") == "file" and "path" in node:
            yield node
        for child in node.values():
            yield from study_files(child)


def metadata(directory):
    raw = (directory / "E-MTAB-694.sdrf.txt").read_bytes()
    study = (directory / "biostudies.json").read_bytes()
    if sha256(raw).hexdigest() != SDRF_SHA256 or sha256(study).hexdigest() != STUDY_SHA:
        raise ValueError("Unexpected study or SDRF version")
    rows = list(csv.DictReader(io.StringIO(raw.decode("ascii")), delimiter="\t"))
    summarize(rows)
    manifest = [f for f in study_files(json.loads(study)) if f["path"].lower().endswith(".cel")]
    sizes = {f["path"]: int(f["size"]) for f in manifest}
    files = [r["Array Data File"] for r in rows]
    if (len(manifest) != len(sizes) or len(set(files)) != len(rows) or set(files) != sizes.keys()
            or any(Path(f).name != f for f in files)
            or {r["Array Design REF"] for r in rows} != {"A-AFFY-44"}):
        raise ValueError("Ambiguous CEL file identity or design")
    return rows, sizes


def fetch_prefix(directory, filename, size):
    path = directory / (filename + ".prefix")
    if path.exists():
        return
    end = HEADER_BYTES - 1
    req = Request(BASE + quote(filename, safe=""), headers={"Range": f"bytes=0-{end}", "Accept-Encoding": "identity"})
    with urlopen(req, timeout=30) as response:
        if response.status != 206 or response.headers.get("Content-Range") != f"bytes 0-{end}/{size}":
            raise ValueError("Server did not honor the bounded metadata request")
        part = response.read(HEADER_BYTES + 1)
    scan(part)
    path.write_bytes(part)


def scan(prefix):
    if len(prefix) != HEADER_BYTES:
        raise ValueError("Wrong prefix length")
    magic, version, cols, rows, cells, header_length = struct.unpack_from("<6i", prefix)
    if (magic, version, cols, rows, cells) != (64, 4, 1164, 1164, 1164**2) or header_length < len(prefix)-24:
        raise ValueError("Unexpected format, geometry or header length")
    lines = [line for line in prefix[24:].splitlines() if line.startswith(b"DatHeader=")]
    if len(lines) != 1 or b"HG-U133_Plus_2.1sq" not in lines[0]:
        raise ValueError("Missing U133 Plus 2 scanner header")
    fields = re.findall(rb"\b(\d{2}/\d{2}/\d{2} \d{2}:\d{2}:\d{2})\s+(\d+)\s+M10\b", lines[0])
    if len(fields) != 1:
        raise ValueError("Ambiguous acquisition timestamp or scanner")
    stamp, scanner = fields[0]
    return datetime.strptime(stamp.decode(), "%m/%d/%y %H:%M:%S"), scanner


def scan_design(records):
    donors = defaultdict(dict)
    blocks = defaultdict(set)
    seen = set()
    for donor, tissue, stamp, scanner in records:
        if (stamp, scanner) in seen or tissue in donors[donor]:
            raise ValueError("Repeated scanner/timestamp or donor-tissue record")
        seen.add((stamp, scanner))
        donors[donor][tissue] = stamp, scanner
        blocks[(stamp.date(), scanner)].add(tissue)
    lesion, peritoneum, _ = TISSUES
    pairs = [(t[lesion], t[peritoneum]) for t in donors.values() if lesion in t and peritoneum in t]
    shared = [(a, b) for a, b in pairs if a[1] == b[1] and a[0].date() == b[0].date()]
    by_scanner = defaultdict(Counter)
    for _, tissue, _, scanner in records:
        by_scanner[scanner][tissue] += 1
    day_gaps = [(b[0].date() - a[0].date()).days for a, b in pairs]

    scanner_counts = sorted(by_scanner.values(), key=lambda counts: tuple(counts[tissue] for tissue in TISSUES))

    def order(group):
        return {"lesion_first": sum(a[0] < b[0] for a, b in group),
                "peritoneum_first": sum(a[0] > b[0] for a, b in group),
                "equal_clock_time": sum(a[0] == b[0] for a, b in group)}
    return {
        "arrays": len(records), "donors": len(donors),
        "physical_scanners": len({scanner for _, scanner in seen}), "scanner_days": len(blocks),
        "scanner_days_with_lesion_and_peritoneum": sum({lesion, peritoneum} <= ts for ts in blocks.values()),
        "paired_donors": len(pairs), "same_scanner_day_pairs": len(shared),
        "same_scanner_pairs": sum(a[1] == b[1] for a, b in pairs),
        "recorded_clock_order_all_paired_donors": order(pairs),
        "peritoneum_minus_lesion_recorded_date_gap_days_range": [min(day_gaps), max(day_gaps)] if pairs else [],
        "tissue_counts_by_scanner": [dict(counts) for counts in scanner_counts],
        "order_same_scanner_day": order(shared),
        "same_day_pairs_per_scanner_day_sorted": sorted(Counter((a[0].date(), a[1]) for a, _ in shared).values()),
    }


def audit(directory, fetch=False):
    rows, sizes = metadata(directory)
    if fetch:
        with ThreadPoolExecutor(max_workers=3) as pool:
            list(pool.map(lambda f: fetch_prefix(directory, f, sizes[f]), sorted(sizes)))
    records = []
    digest = sha256()
    for row in rows:  # stable public SDRF order
        prefix = (directory / (row["Array Data File"] + ".prefix")).read_bytes()
        verify_dat_filename(prefix, row["Array Data File"])
        digest.update(prefix)
        records.append((row["Characteristics[Individual]"], row["Characteristics[OrganismPart]"], *scan(prefix)))
    if digest.hexdigest() != PREFIXES_SHA:
        raise ValueError("Native header digest differs from the pinned SDRF-order census")
    ann = directory / "GPL570.annot.gz"
    if sha256_file(ann) != ANNOTATION_SHA:
        raise ValueError("Unexpected GPL570 annotation")
    mapping = read_annotation(ann)
    panel = json.loads((ROOT / "evidence/lesion-transfer/panel.json").read_text())
    probes = {str(g["entrez"]): sorted(p for p, v in mapping.items() if v == g["entrez"]) for g in panel}
    if len(mapping) != 54675 or len(probes) != 5 or not all(probes.values()):
        raise ValueError("Incomplete feature coverage")
    return {"study": "E-MTAB-694", "sdrf_sha256": SDRF_SHA256, "study_json_sha256": STUDY_SHA,
            "prefix_bytes_each": HEADER_BYTES, "sum_prefix_bytes": len(rows) * HEADER_BYTES,
            "sdrf_order_prefix_sha256": digest.hexdigest(), "archive_cel_bytes": sum(sizes.values()),
            "physical_design": "HG-U133_Plus_2", "design": scan_design(records),
            "conditional_design": conditional_design(records),
            "gpl570_sha256": ANNOTATION_SHA, "total_probe_sets": len(mapping), "clean_probes_by_entrez": probes}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--fetch-headers", action="store_true")
    args = parser.parse_args()
    print(json.dumps(audit(args.directory, args.fetch_headers), indent=2, sort_keys=True))
