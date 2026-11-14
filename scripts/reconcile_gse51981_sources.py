"""Reconcile public GEO sample metadata with the GSE51981 Series Matrix header.

Print aggregate field counts. No accession, source text, reviewed manifest or
individual-level assignment is emitted. The GEO quick view embeds table
previews; their contents are skipped. This command makes no network requests.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path

from audit_gse51981_metadata import ACCESSION, GSM, _read_metadata, sha256_file

CHAR = "!Sample_characteristics_ch1"
SOURCE_NAME = "!Sample_source_name_ch1"
TISSUE = "tissue"
DIAGNOSIS = "endometriosis/no endometriosis"
PATHOLOGY = "presence or absence of uterine/pelvic pathology"
SEVERITY = "endometriosis severity"
PHASE_CODES = {
    "Proliferative Endometrial tissue": "proliferative",
    "Early Secretory Endometrial tissue": "early_secretory",
    "Mid-Secretory Endometrial tissue": "mid_secretory",
    "Late Secretory Endometrial tissue": "late_secretory",
    "Unknown": "unknown",
}
GROUP_CODES = {
    ("Non-Endometriosis", "No Uterine Pelvic Pathology"): "control",
    ("Non-Endometriosis", "Uterine Pelvic Pathology"): "other",
}
SHARED = ("!Sample_platform_id", SOURCE_NAME, CHAR)


def read_quick_view(source: Path) -> dict:
    """Read SOFT metadata, discarding all table contents and unrelated fields."""
    records = {}
    current = None
    in_table = False
    with source.open(encoding="utf-8-sig") as stream:
        for line in stream:
            line = line.rstrip("\r\n")
            if in_table:
                if line == "!sample_table_end":
                    in_table = False
                continue
            if line == "!sample_table_begin":
                in_table = True
                continue
            key, sep, value = line.partition(" = ")
            if key == "^SAMPLE":
                if not sep or not GSM.fullmatch(value) or value in records:
                    raise ValueError("Invalid or duplicate SOFT sample identifier")
                current = records[value] = {}
            elif current is not None and key in (
                *SHARED, "!Sample_geo_accession", "!Sample_series_id", "#VALUE"
            ):
                current.setdefault(key, []).append(value)
    if in_table:
        raise ValueError("Unterminated SOFT table")
    return records


def reconcile(matrix: Path, quick_view: Path) -> dict:
    samples, meta, coverage = _read_metadata(matrix)
    records = read_quick_view(quick_view)
    if set(samples) != set(records):
        raise ValueError("SOFT and matrix sample sets differ")
    if any(field not in meta for field in SHARED) or len(meta[SOURCE_NAME]) != 1:
        raise ValueError("Missing or repeated matrix source fields")

    counts: dict[str, Counter] = {g: Counter() for g in ("case", "control", "other", "unresolved")}
    no_source_severity = 0
    log2 = 0
    for i, sample in enumerate(samples):
        record = records[sample]
        if record.get("!Sample_geo_accession") != [sample] or record.get("!Sample_series_id") != [ACCESSION]:
            raise ValueError("SOFT sample or series identifier mismatch")
        for field in SHARED:
            if sorted(record.get(field, [])) != sorted(row[i] for row in meta[field]):
                raise ValueError("SOFT and matrix sample field values differ")
        fields = {}
        for entry in record[CHAR]:
            key, sep, value = entry.partition(": ")
            if not sep or key in fields:
                raise ValueError("Malformed or duplicate sample characteristic")
            fields[key] = value

        group = GROUP_CODES.get((fields.get(DIAGNOSIS), fields.get(PATHOLOGY)), "unresolved")
        if fields.get(DIAGNOSIS) == "Endometriosis":
            group = "case"
        phase = PHASE_CODES.get(fields.get(TISSUE), "unrecognized_or_missing")
        counts[group][phase] += 1
        if (group == "case" and fields.get(SEVERITY)
                and record[SOURCE_NAME] == ["Endometriosis_Severity Not Available"]):
            no_source_severity += 1
        log2 += record.get("#VALUE") == ["base-2 logarithm GC-RMA"]
    return {
        "series": ACCESSION,
        "reconciled_arrays": len(samples),
        "group_by_source_phase": counts,
        "cases_with_severity_characteristic_but_source_name_severity_unavailable": no_source_severity,
        "log2_gcrma_value_descriptions": log2,
        "source_field_coverage": coverage,
        "matrix_input_bytes": matrix.stat().st_size,
        "matrix_input_sha256": sha256_file(matrix),
        "quick_view_bytes": quick_view.stat().st_size,
        "quick_view_sha256": sha256_file(quick_view),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True, help="Matrix or gzip prefix containing its complete header")
    parser.add_argument("--soft", type=Path, required=True, help="GEO accession viewer's sample-targeted SOFT quick view")
    args = parser.parse_args()
    print(json.dumps(reconcile(args.matrix, args.soft), sort_keys=True, indent=2))
