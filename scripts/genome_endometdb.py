"""Read genome-wide Finnish WG-6 ranks locally from the audited source release.

The load_finland API emits no identifiers and never writes matrices to disk.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path
import re

import numpy as np

from analyze_endometdb_wg6_fosb import N_FEATURES, verified_wg6_source
from audit_endometdb import pinned
from audit_endometdb_probe_keys import MANIFESTS


WG6_MANIFEST = "GPL13376-platform.gz"


def annotation_entrez(rows):
    """Represent every platform address; zero marks a reference-only probe."""
    mapping = {}
    for row in rows:
        address_text = row.get("Array_Address_Id") or ""
        if not address_text.isdecimal():
            continue
        address = int(address_text)
        if address in mapping:
            raise ValueError("Duplicate physical address in WG-6 annotation")
        entrez = row.get("Entrez_Gene_ID") or ""
        valid = (row["Species"] == "Homo sapiens" and row["Source"] == "RefSeq"
                 and re.fullmatch(r"[1-9][0-9]*", entrez))
        mapping[address] = int(entrez) if valid else 0
    return mapping


def load_wg6_annotation(directory):
    path = directory / WG6_MANIFEST
    pinned(path, MANIFESTS[WG6_MANIFEST])
    with gzip.open(path, "rt", newline="", encoding="utf-8") as fh:
        for line in fh:
            if "Probe_Id\tArray_Address_Id\t" in line:
                break
        else:
            raise ValueError("Missing WG-6 annotation header")
        fields = line.rstrip("\r\n").split("\t")
        required = {"Array_Address_Id", "Species", "Source", "Entrez_Gene_ID"}
        if not required.issubset(fields) or len(set(fields)) != len(fields):
            raise ValueError("WG-6 annotation fields changed")
        mapping = annotation_entrez(csv.DictReader(fh, delimiter="\t", fieldnames=fields))
    counts = (len(mapping), sum(g > 0 for g in mapping.values()), len(set(mapping.values()) - {0}))
    if counts != (N_FEATURES, 30_541, 24_963):
        raise ValueError("WG-6 Entrez mapping differs from the audited release")
    return mapping


def finland_annotation_ids(directory):
    """Return eligible Entrez IDs from the pinned annotation without opening clinical or expression files."""
    mapping = load_wg6_annotation(Path(directory))
    return np.array(sorted(set(mapping.values()) - {0}), dtype=np.int64)


def selected_signals(lines, design, annotation):
    """Read finite signals in aligned, local-only arrays without exposing IDs."""
    addresses = sorted(annotation)
    columns = {address: col for col, address in enumerate(addresses)}
    assays = sorted(design)
    rows = {assay: row for row, assay in enumerate(assays)}
    signals = np.full((len(assays), len(addresses)), np.nan)
    for line in lines:
        _, assay, rest = line.split(b"\t", 2)
        row = rows.get(assay)
        if row is None:
            continue
        fields = rest.split(b"\t")
        if len(fields) != 5 or not fields[0].isdigit():
            raise ValueError("Malformed selected WG-6 record")
        col = columns.get(int(fields[0]))
        if col is None or not np.isnan(signals[row, col]):
            raise ValueError("Unexpected or duplicated physical address in selected WG-6 assay")
        try:
            signal = float(fields[1])
        except ValueError:
            raise ValueError("Missing or malformed selected WG-6 signal") from None
        if not np.isfinite(signal):
            raise ValueError("Nonfinite selected WG-6 signal")
        signals[row, col] = signal
    if np.isnan(signals).any():
        raise ValueError("Incomplete selected WG-6 assay")
    cases = np.array([design[assay][1] == "case" for assay in assays], dtype=bool)
    strata = np.array(["-".join(design[assay][0]) for assay in assays])
    probe_genes = np.array([annotation[address] for address in addresses], dtype=np.int64)
    return signals, probe_genes, cases, strata


def load_finland(directory):
    """Return sorted Entrez IDs, donor × gene percentiles, case flags and run-phase strata."""
    from open_disease.microarray_ranks import gene_percentiles

    directory = Path(directory)
    annotation = load_wg6_annotation(directory)
    design, physical_addresses, lines = verified_wg6_source(directory)
    if annotation.keys() != physical_addresses:
        raise ValueError("Annotation and census feature sets disagree")
    signals, probe_genes, cases, strata = selected_signals(lines, design, annotation)
    gene_ids, scores = gene_percentiles(signals, probe_genes)
    return gene_ids, scores, cases, strata
