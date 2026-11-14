"""Load the Burney proliferative endometrial arrays and audit physical scan metadata."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import struct
import urllib.request
import zlib

import numpy as np
from scipy.stats import rankdata

from genome_rome import read_annotation
from open_disease.geo_microarray import sha256_file
from open_disease.microarray_ranks import gene_percentiles

MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE6nnn/GSE6364/matrix/GSE6364_series_matrix.txt.gz"
BRIEF_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6364&targ=gsm&view=brief&form=text"
MATRIX_SHA = "b96db29937984402dc4fd13efdf6c609d1b58c007d11336a533fef1c86ccf7da"
BRIEF_SHA = "663f7c276d364c252b36b3a2dc8d89eee37a9cac250ce662197506746c44e413"
GPL_SHA = "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394"
PREFIX_SHA = "9906fdd3faf7e2574d252aab6e672c9a8d61875d5f1961e3c29c9e491479ee37"
EXPECTED = {("PE", "D"): 6, ("PE", "N"): 5, ("ESE", "D"): 6,
            ("ESE", "N"): 3, ("MSE", "D"): 9, ("MSE", "N"): 8}
PHASE = {"PE": "Proliferative", "ESE": "Early Secretory", "MSE": "Mid Secretory"}


def source_header(reader):
    header = {}
    for row in reader:
        if row == ["!series_matrix_table_begin"]:
            return header
        if row and row[0] in {"!Series_geo_accession", "!Series_platform_id", "!Sample_geo_accession",
                              "!Sample_platform_id", "!Sample_title", "!Sample_characteristics_ch1",
                              "!Sample_source_name_ch1", "!Sample_data_processing", "!Sample_supplementary_file"}:
            if row[0] in header:
                raise ValueError("Repeated GSE6364 sample header")
            header[row[0]] = row[1:]
    raise ValueError("Missing GSE6364 table boundary")


def selection(header):
    ids = header["!Sample_geo_accession"]
    if (header["!Series_geo_accession"] != ["GSE6364"] or header["!Series_platform_id"] != ["GPL570"]
            or len(ids) != 37 or len(set(ids)) != 37 or header["!Sample_platform_id"] != ["GPL570"] * 37
            or any(not re.fullmatch(r"GSM[0-9]+", key) for key in ids)):
        raise ValueError("GSE6364 series, platform or array identifiers changed")
    fields = ("!Sample_title", "!Sample_characteristics_ch1", "!Sample_source_name_ch1",
              "!Sample_data_processing", "!Sample_supplementary_file")
    seen, cohort, selected = set(), Counter(), []
    for i, (title, characteristic, tissue, processing, url) in enumerate(zip(*(header[k] for k in fields), strict=True)):
        match = re.fullmatch(r"(PE|ESE|MSE)_(D|N)_([A-Za-z0-9]+)", title)
        if not match:
            raise ValueError("Unrecognized case, phase or patient fields")
        phase, status, key = match.groups()
        expected = PHASE[phase] + " Phase " + ("Endometriosis" if status == "D" else "Normal")
        expected_source = ("uterus - " + PHASE[phase].capitalize() + " endometrium from "
                           + ("patient with endometriosis" if status == "D" else "normal patient"))
        if (key in seen or expected != characteristic or " ".join(tissue.split()) != expected_source
                or "Gene Chip Operating Software Version 1.2" not in processing
                or "trimmed mean target intensity of each array was arbitrarily set to 500" not in processing
                or not url.startswith("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/") or not url.endswith(".CEL.gz")):
            raise ValueError("Conflicting specimen fields, platform processing or CEL link")
        seen.add(key)
        cohort[phase, status] += 1
        if phase == "PE":
            selected.append((i, ids[i], status == "D", url.replace("ftp://", "https://", 1)))
    if cohort != EXPECTED:
        raise ValueError("Clinical group or menstrual phase counts changed")
    return selected


def verify_brief(directory, header):
    path = directory / "brief.soft"
    if sha256_file(path) != BRIEF_SHA:
        raise ValueError("Unexpected GSE6364 brief sample metadata")
    keys = ("!Sample_title", "!Sample_geo_accession", "!Sample_characteristics_ch1",
            "!Sample_source_name_ch1", "!Sample_supplementary_file")
    records = defaultdict(dict)
    current = None
    for line in path.read_text(encoding="utf-8").splitlines():
        name, sep, value = line.partition(" = ")
        if name == "^SAMPLE" and sep:
            current = value
        elif name in keys:
            if current is None or name in records[current]:
                raise ValueError("Malformed brief sample records")
            records[current][name] = value
    ids = header["!Sample_geo_accession"]
    if set(records) != set(ids) or any([records[key].get(k) for key in ids] != header[k] for k in keys):
        raise ValueError("Series matrix and GSM metadata differ")


def checked_header(directory):
    path = Path(directory) / "series.txt.gz"
    if sha256_file(path) != MATRIX_SHA:
        raise ValueError("Unrecognized GSE6364 matrix")
    with gzip.open(path, "rt", newline="") as source:
        header = source_header(csv.reader(source, delimiter="\t"))
    verify_brief(Path(directory), header)
    return header, selection(header)


def scanner_fields(prefix):
    """Read the dated DAT scanner fields of a binary XDA v4 CEL, without reading probe values."""
    if len(prefix) != 2048:
        raise ValueError("Expected a 2,048-byte compressed CEL prefix")
    raw = zlib.decompressobj(31).decompress(prefix, 8192)
    if len(raw) < 24:
        raise ValueError("Short CEL header")
    magic, version, cols, rows, n, length = struct.unpack_from("<6i", raw)
    if (magic, version, cols, rows, n) != (64, 4, 1164, 1164, 1164**2) or not 0 < length <= len(raw)-24:
        raise ValueError("Unrecognized XDA geometry or header length")
    lines = [s for s in raw[24:24+length].decode("latin-1").splitlines() if s.startswith("DatHeader=")]
    if len(lines) != 1 or "HG-U133_Plus_2.1sq" not in lines[0]:
        raise ValueError("Expected one GPL570 native DAT header")
    match = re.search(r"\b(\d{2}/\d{2}/\d{2}) (\d{2}:\d{2}:\d{2})\s+(\d+)\s+M10\b", lines[0])
    if not match:
        raise ValueError("Unsupported native scanner date or identifier")
    return datetime.strptime(match[1], "%m/%d/%y").date(), match[3]


def fetch_prefix(url):
    if not url.startswith("https://ftp.ncbi.nlm.nih.gov/geo/samples/") or not url.endswith(".CEL.gz"):
        raise ValueError("Unexpected GEO CEL source")
    request = urllib.request.Request(url, headers={"Range": "bytes=0-2047"})
    with urllib.request.urlopen(request, timeout=45) as response:
        if response.status != 206 or not re.fullmatch(r"bytes 0-2047/\d+", response.headers.get("Content-Range", "")):
            raise ValueError("Server returned an unexpected CEL range")
        part = response.read(2049)
    if len(part) != 2048:
        raise ValueError("Incomplete compressed CEL prefix")
    return part


def scan_audit(directory, selected, fetch=False):
    directory = Path(directory)
    digest = hashlib.sha256()
    scanner_ids, by_date = set(), defaultdict(Counter)
    held = []
    for _, accession, case, url in selected:
        path = directory / (accession + ".header.gzpart")
        if path.exists():
            part = path.read_bytes()
        elif fetch:
            part = fetch_prefix(url)
            held.append((path, part))
        else:
            raise ValueError("Missing native CEL scan metadata. Use --fetch-headers once.")
        date, scanner = scanner_fields(part)
        by_date[date][case] += 1
        scanner_ids.add(scanner)
        digest.update(part)
    if digest.hexdigest() != PREFIX_SHA:
        raise ValueError("Native CEL metadata digest differs")
    for path, part in held:
        path.write_bytes(part)
    counts = sorted((groups[True], groups[False]) for groups in by_date.values())
    return {"arrays": len(selected), "case_arrays": sum(x[2] for x in selected),
            "control_arrays": sum(not x[2] for x in selected), "native_scanners": len(scanner_ids),
            "scan_days": len(counts), "mixed_scan_days": sum(a > 0 and b > 0 for a,b in counts),
            "day_sizes_cases": sorted(a for a,b in counts if a), "day_sizes_controls": sorted(b for a,b in counts if b),
            "conditional_assignments": int(np.prod([math.comb(a+b, a) for a,b in counts])),
            "range_bytes": len(selected)*2048, "ordered_range_sha256": digest.hexdigest()}


def checked_platform(annotation):
    if sha256_file(annotation) != GPL_SHA:
        raise ValueError("Unexpected GPL570 annotation")
    mapping = read_annotation(annotation)
    if (len(mapping), sum(g > 0 for g in mapping.values()), len(set(mapping.values()) - {0})) != (54675, 42891, 20845):
        raise ValueError("Unexpected GPL570 feature or gene counts")
    return mapping


def read_eligible_signals(directory, annotation):
    directory = Path(directory)
    header, selected = checked_header(directory)
    scan = scan_audit(directory, selected)
    mapping = checked_platform(annotation)
    probe_ids, entrez, signals = set(), [], []
    with gzip.open(directory / "series.txt.gz", "rt", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        source_header(reader)
        if next(reader, None) != ["ID_REF", *header["!Sample_geo_accession"]]:
            raise ValueError("Assay column order differs from clinical records")
        for row in reader:
            if row == ["!series_matrix_table_end"]:
                break
            if len(row) != 38 or row[0] not in mapping or row[0] in probe_ids:
                raise ValueError("Missing or duplicate feature in GSE6364")
            probe_ids.add(row[0])
            entrez.append(mapping[row[0]])
            signals.append([float(row[idx+1]) for idx, *_ in selected])
        else:
            raise ValueError("Missing end of GSE6364 matrix")
        if any(any(row) for row in reader):
            raise ValueError("Unexpected trailing GSE6364 data")
    if probe_ids != mapping.keys() or not np.isfinite(signals).all():
        raise ValueError("Reference features missing or nonfinite in GSE6364")
    return np.array(signals).T, np.array(entrez), np.array([s[2] for s in selected]), scan


def load_us(private_study_dir, gpl570_path):
    values, entrez, cases, _ = read_eligible_signals(private_study_dir, gpl570_path)
    genes, gene_scores = gene_percentiles(values, entrez)
    return genes, gene_scores, cases


def load_us_both(private_study_dir, gpl570_path, gene_ids):
    """Return (genes, mean percentiles, case labels, median percentiles) and the scan-day audit."""
    values, entrez, cases, scan = read_eligible_signals(private_study_dir, gpl570_path)
    ranks = (rankdata(values, method="average", axis=1) - .5) * (100 / len(entrez))
    genes = np.unique(np.asarray(gene_ids, dtype=int))
    if not len(genes) or (genes <= 0).any() or not np.isin(genes, entrez).all():
        raise ValueError("A requested pathway gene is absent from GPL570")
    mean = np.empty((len(cases), len(genes)))
    median = np.empty_like(mean)
    for column, gene in enumerate(genes):
        per_probe = ranks[:, entrez == gene]
        mean[:, column] = per_probe.mean(axis=1)
        median[:, column] = np.median(per_probe, axis=1)
    return (genes, mean, cases, median), scan


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("annotation", type=Path)
    parser.add_argument("--fetch-headers", action="store_true", help="Acquire 2,048 bytes per eligible public CEL")
    args = parser.parse_args()
    _, selected = checked_header(args.directory)
    audit = scan_audit(args.directory, selected, args.fetch_headers)
    annotation = checked_platform(args.annotation)
    audit.update(reference_probes=len(annotation), mapped_probes=sum(g>0 for g in annotation.values()),
                 mapped_genes=len(set(annotation.values()) - {0}), matrix_sha256=MATRIX_SHA,
                 annotation_sha256=GPL_SHA, individual_geo_metadata_sha256=BRIEF_SHA)
    print(json.dumps(audit, indent=2))
