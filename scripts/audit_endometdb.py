"""Audit the public EndometDB v1 clinical tables locally, before reading expression.

Read two checksum-pinned byte ranges of the authors' Figshare ZIP. Only the
listed PostgreSQL COPY payloads are decoded. The default output is aggregate
JSON. The optional manifest is local to the ignored data/private directory.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import zlib
from collections import Counter, defaultdict
from pathlib import Path


# These ZIP/Deflate64 restart and pg_dump offsets are specific to this release.
# The range hashes are checked before any offset is used.
HEAD_SHA256 = "023a9ee95ecb6901cac61f1fdf82ed699e1c8d2159e61665d2ef499ebd8708fa"
TAIL_SHA256 = "75f0e4ce7f600af1716c16b5b551331725d437f0608b92dc4544024f2a351fc8"
DUMP_LENGTH = 2957024788
TABLES = (
    "code_value", "patient", "patient_sample", "patient_survey",
    "pathology_report_endometrium", "sample_analysis_result",
)
PHASES = {"Proliferative": "proliferative", "Secretory": "secretory"}
GROUPS = {("patient", "endometriosis"): "case", ("control", "healthy"): "control"}


def pinned(path, digest):
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != digest:
        raise ValueError("Archive range does not match the pinned EndometDB v1 bytes")
    return data


def read_metadata(head_file, tail_file):
    from inflate64 import Inflater

    head = Inflater().inflate(pinned(head_file, HEAD_SHA256)[48:])
    tail = Inflater().inflate(pinned(tail_file, TAIL_SHA256)[157851:])
    tail_start = DUMP_LENGTH - len(tail)
    tables = {}
    for table in TABLES:
        # Each custom-archive TOC entry stores its COPY columns and a file offset.
        marker = re.search(rb"COPY public\." + table.encode() + rb" \(([^)]+)\) FROM stdin;", head)
        offset = re.search(rb"\x01\x01\x00\x00\x00\x02(.{8})", head[marker.end():marker.end()+240], re.S)
        position = int.from_bytes(offset[1], "little")
        buf = head
        if position >= tail_start:
            position -= tail_start
            buf = tail
        columns = marker[1].decode().split(", ")
        if buf[position:position+2] != b"\x01\x00":
            raise ValueError("Expected a PostgreSQL custom archive data block")
        position += 6  # block tag and signed four-byte dump ID
        chunks = bytearray()
        while position + 5 <= len(buf) and buf[position] == 0:
            size = int.from_bytes(buf[position+1:position+5], "little")
            position += 5
            if size == 0:
                break
            chunks.extend(buf[position:position+size])
            position += size
        else:
            raise ValueError("Incomplete PostgreSQL COPY block")
        rows = zlib.decompress(chunks).decode().splitlines()
        tables[table] = [dict(zip(columns, row.split("\t"), strict=True))
                         for row in rows if row and row != r"\."]
    return tables


def index(rows, key):
    result = {r[key]: r for r in rows}
    if len(result) != len(rows):
        raise ValueError(f"Nonunique metadata join on {key}")
    return result


def eligible(tables):
    """Join within the native participant namespace and exclude repeated donors."""
    code = {k: row["title_en"] for k, row in index(tables["code_value"], "id").items()}
    patients = index(tables["patient"], "id")
    samples = index(tables["patient_sample"], "sample_id")
    surveys = index(tables["patient_survey"], "patient_id")
    histology = index(tables["pathology_report_endometrium"], "patient_sample_id")
    assays = index(tables["sample_analysis_result"], "id")
    by_donor = defaultdict(list)
    endometrium = 0
    for assay in assays.values():
        sample = samples.get(assay["sample_id"])
        if not sample or code[sample["code_tissue_type_id"]] != "Endometrium tissue":
            continue
        endometrium += 1
        donor = patients[sample["patient_id"]]
        group = GROUPS.get((code[donor["code_patient_category_id"]], code[donor["code_confirmed_disease_status_id"]]))
        if not group or surveys.get(donor["id"], {}).get("hormonal_medication") != "f":
            continue
        record = histology.get(sample["id"], {})
        phase = PHASES.get(code.get(record.get("code_histology_class_id")))
        if phase:
            by_donor[donor["id"]].append({
                "donor_id": donor["id"], "analysis_result_id": assay["id"],
                "analysis_run_id": assay["analysis_run_id"], "group": group, "phase": phase,
            })
    selected = [rows[0] for rows in by_donor.values() if len(rows) == 1]
    counts = Counter((row["group"], row["phase"]) for row in selected)
    strata = defaultdict(Counter)
    for row in selected:
        strata[row["analysis_run_id"], row["phase"]][row["group"]] += 1
    overlap = [c for c in strata.values() if c["control"] and c["case"]]
    report = {
        "native_participant_count": len(patients),
        "native_assays_total": len(assays),
        "endometrium_assays": endometrium,
        "survey_hormonal_medication": dict(Counter(r["hormonal_medication"] for r in surveys.values())),
        "participants_with_no_survey": len(patients.keys() - surveys.keys()),
        "candidate_assays": sum(map(len, by_donor.values())),
        "candidate_donors": len(by_donor),
        "repeated_donors_excluded": len(by_donor) - len(selected),
        "single_assay_donors": len(selected),
        "group_phase": {g: {p: counts[g, p] for p in PHASES.values()} for g in ("case", "control")},
        "shared_run_phase_strata": len(overlap),
        "single_assay_donors_in_shared_strata": {g: sum(c[g] for c in overlap) for g in ("case", "control")},
    }
    return report, selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--head", type=Path, required=True)
    parser.add_argument("--tail", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, help="Optional new private CSV of eligible native IDs")
    args = parser.parse_args()
    report, rows = eligible(read_metadata(args.head, args.tail))
    if args.manifest:
        private = Path(__file__).resolve().parents[1] / "data" / "private"
        if not args.manifest.resolve().is_relative_to(private) or args.manifest.suffix != ".csv":
            parser.error("Manifest must be a new CSV within the project's data/private directory")
        with args.manifest.open("x", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=("donor_id", "analysis_result_id", "analysis_run_id", "group", "phase"))
            w.writeheader()
            w.writerows(rows)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
