"""Load the pinned Rome uterine-biopsy cohort for local pathway validation."""
from __future__ import annotations

from collections import defaultdict
import csv
import gzip
from pathlib import Path
import re

import numpy as np

from analyze_gse25628_fosb import ANNOTATION_HASH, MATRIX_HASH, checked_groups
from open_disease.geo_microarray import sha256_file
from open_disease.microarray_ranks import gene_percentiles


def read_annotation(path):
    """Map every GPL571 feature to one unambiguous Entrez ID or to zero."""
    mapping = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        for row in reader:
            if row == ["!platform_table_begin"]:
                break
        header = next(reader, [])
        columns = ("ID", "Gene symbol", "Gene ID")
        if any(header.count(name) != 1 for name in columns):
            raise ValueError("Missing or ambiguous GPL571 annotation header")
        probe, symbol, entrez = (header.index(name) for name in columns)
        for row in reader:
            if row == ["!platform_table_end"]:
                break
            if len(row) != len(header) or not row[probe] or row[probe] in mapping:
                raise ValueError("Invalid or repeated platform feature")
            valid = (not row[probe].startswith("AFFX")
                     and re.fullmatch(r"[1-9][0-9]*", row[entrez])
                     and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", row[symbol]))
            mapping[row[probe]] = int(row[entrez]) if valid else 0
        else:
            raise ValueError("Missing end of GPL571 annotation")
        if any(any(row) for row in reader):
            raise ValueError("Trailing GPL571 annotation content")
    return mapping


def read_matrix(path, mapping, min_complete):
    """Select independent uterine arrays and retain a common complete feature set."""
    with gzip.open(path, "rt", encoding="utf-8", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = defaultdict(list)
        for row in reader:
            if row == ["!series_matrix_table_begin"]:
                break
            if row and row[0].startswith(("!Sample_", "!Series_")):
                header[row[0]].append(row[1:])
        ids, groups = checked_groups(header)
        if next(reader, None) != ["ID_REF", *ids]:
            raise ValueError("Rome expression columns disagree with the audited GEO order")
        included = [i for i, group in enumerate(groups) if group != "lesion"]
        is_case = np.asarray([groups[i] == "case" for i in included])

        seen = set()
        features = []
        values = []
        missing = {"", "na", "null"}
        for row in reader:
            if row == ["!series_matrix_table_end"]:
                break
            if len(row) != len(ids) + 1 or row[0] not in mapping or row[0] in seen:
                raise ValueError("Invalid or repeated Rome matrix feature")
            seen.add(row[0])
            try:
                signals = [np.nan if row[i + 1].strip().lower() in missing
                           else float(row[i + 1]) for i in included]
            except ValueError:
                raise ValueError("Malformed value in an eligible Rome array") from None
            if np.isfinite(signals).all():
                values.append(signals)
                features.append(mapping[row[0]])
            elif mapping[row[0]]:
                raise ValueError("Incomplete eligible Rome probe would change the frozen gene measurement")
        else:
            raise ValueError("Missing end of Rome matrix")
        if any(any(row) for row in reader):
            raise ValueError("Trailing Rome matrix content")
    if seen != mapping.keys() or len(features) < min_complete:
        raise ValueError("Incomplete Rome platform or insufficient common reference features")
    return np.array(values, dtype=np.float64).T, np.array(features, dtype=np.int64), is_case


def checked_annotation(directory):
    annotation = Path(directory) / "GPL571.annot.gz"
    if sha256_file(annotation) != ANNOTATION_HASH:
        raise ValueError("Rome annotation SHA-256 differs from the audited snapshot")
    mapping = read_annotation(annotation)
    if (len(mapping) != 22_277 or sum(gene > 0 for gene in mapping.values()) != 19_921
            or len(set(mapping.values()) - {0}) != 12_502):
        raise ValueError("Rome platform annotation counts differ from the audited snapshot")
    return mapping


def rome_annotation_ids(directory):
    """Return sorted Entrez IDs from platform annotation alone, before validation."""
    return np.array(sorted(set(checked_annotation(directory).values()) - {0}), dtype=np.int64)


def load_rome(directory):
    """Return gene IDs, donor-by-gene scores, case labels, and one Rome stratum."""
    matrix = Path(directory) / "GSE25628_series_matrix.txt.gz"
    if sha256_file(matrix) != MATRIX_HASH:
        raise ValueError("Rome matrix SHA-256 differs from the audited snapshot")
    mapping = checked_annotation(directory)
    values, entrez, cases = read_matrix(matrix, mapping, min_complete=20_000)
    genes, scores = gene_percentiles(values, entrez)
    return genes, scores, cases, np.repeat("proliferative", len(cases))
