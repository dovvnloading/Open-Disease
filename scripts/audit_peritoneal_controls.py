"""Count peritoneal tissue comparators in the pinned EndometDB v1 metadata.

Read only the same six clinical/assay-index tables as audit_endometdb.py.
No expression table is restored. All joins use one native identifier namespace.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
import json
from pathlib import Path
import re

from audit_endometdb import GROUPS, HEAD_SHA256, PHASES, TAIL_SHA256, index, pinned, read_metadata


ARMS = ("lesion", "case_peritoneum", "control_peritoneum")
UTERUS = "control_endometrium"
# Figshare v1 codebook values: each color has replicate labels 1, 2 and 3.
LESIONS = {color + str(number) for color in ("PU", "MU", "VA") for number in (1, 2, 3)}
MEDICATION = {"f": "stored_false", "t": "stored_true"}


@dataclass(frozen=True)
class Assay:
    arm: str
    donor: str
    specimen: str
    run: str
    endometrial_histology: str | None
    assay_id: str = ""


def paired_capacity(assays, first, second, medication):
    by_donor = defaultdict(lambda: defaultdict(list))
    for a in assays:
        if a.arm in (first, second):
            by_donor[a.donor][a.arm].append(a)
    paired = {d: arms for d, arms in by_donor.items() if first in arms and second in arms}
    common_run = single_each = single_common_run = 0
    paired_medication = Counter()
    phase_counts = Counter()
    single_specimen_phase_counts = Counter()
    one_common_run_one_specimen_each = 0
    for donor, arms in paired.items():
        left, right = arms[first], arms[second]
        paired_medication[medication[donor]] += 1
        one_each = len(left) == len(right) == 1
        single_each += one_each
        common = {a.run for a in left} & {a.run for a in right}
        if not common:
            continue
        common_run += 1
        single_common_run += one_each
        # Restrict to a single common run and a single physical specimen in
        # each arm within it. Replicate assay rows of either specimen may exist.
        unique_specimens = len(common) == 1 and all(
            len({a.specimen for a in records if a.run in common}) == 1
            for records in (left, right)
        )
        one_common_run_one_specimen_each += unique_specimens
        if second == UTERUS and medication[donor] == "stored_false":
            # Require one consistent cycling phase across *all* uterine assays
            # in matching runs. Multiple assays from one donor contribute once.
            phases = {a.endometrial_histology for a in right if a.run in common}
            if len(phases) == 1 and phases <= PHASES.keys():
                phase = PHASES[phases.pop()]
                phase_counts[phase] += 1
                if unique_specimens:
                    single_specimen_phase_counts[phase] += 1
    result = {
        "donors_with_both_labels": len(paired),
        "donors_with_one_assay_each": single_each,
        "donors_with_a_shared_run": common_run,
        "donors_with_one_assay_each_and_shared_run": single_common_run,
        "donors_with_one_common_run_one_specimen_each_within_it": one_common_run_one_specimen_each,
        "hormonal_medication_among_paired_donors": dict(paired_medication),
    }
    if second == UTERUS:
        result["stored_false_shared_run_with_consistent_matched_uterine_phase"] = dict(phase_counts)
        result["stored_false_one_common_run_one_specimen_each_with_consistent_matched_uterine_phase"] = dict(
            single_specimen_phase_counts)
    return result


def native_assays(tables):
    """Decode labels and clinical joins in the native namespace, without reading signal."""
    codes = index(tables["code_value"], "id")
    patients = index(tables["patient"], "id")
    index(tables["patient_sample"], "id")  # physical specimen key must also be unique
    samples = index(tables["patient_sample"], "sample_id")
    surveys = index(tables["patient_survey"], "patient_id")
    pathology = index(tables["pathology_report_endometrium"], "patient_sample_id")
    medication = {
        key: MEDICATION.get(surveys.get(key, {}).get("hormonal_medication"), "unknown")
        for key in patients
    }

    assays = []
    excluded = Counter()
    results = index(tables["sample_analysis_result"], "id")
    for row in results.values():
        if row["result_data_type"] != "microarray":
            raise ValueError("Unexpected assay type in pinned release")
        sample = samples[row["sample_id"]]
        donor = sample["patient_id"]
        patient = patients[donor]
        category_key = patient["code_patient_category_id"]
        disease_key = patient["code_confirmed_disease_status_id"]
        category = codes[category_key]["title_en"] if category_key not in ("", r"\N", None) else "unknown"
        disease = codes[disease_key]["title_en"] if disease_key not in ("", r"\N", None) else "unknown"
        group = GROUPS.get((category, disease))
        tissue = codes[sample["code_tissue_type_id"]]["value"]
        if tissue in LESIONS and group == "case":
            arm = "lesion"
        elif tissue == "P" and group:
            arm = group + "_peritoneum"
        elif tissue == "E" and group == "control":
            arm = UTERUS
        else:
            if tissue == "P" or tissue in LESIONS:
                kind = "peritoneum" if tissue == "P" else "lesion"
                excluded[f"{kind} / {category} / {disease}"] += 1
            continue
        run = row["analysis_run_id"]
        if not run or run == r"\N":
            raise ValueError("Missing native run; do not treat unknown runs as matched")
        phase_key = pathology.get(sample["id"], {}).get("code_histology_class_id")
        phase = codes[phase_key]["title_en"] if phase_key and phase_key != r"\N" else None
        assays.append(Assay(arm, donor, sample["id"], run, phase, row["id"]))
    return assays, medication, excluded


def summarize(tables):
    assays, medication, excluded = native_assays(tables)
    codes = index(tables["code_value"], "id")
    sample_rows = index(tables["patient_sample"], "id")
    pathology = index(tables["pathology_report_endometrium"], "patient_sample_id")
    tissue_histology = Counter(codes[sample_rows[key]["code_tissue_type_id"]]["title_en"]
                               for key in pathology)
    groups = {}
    for arm in (*ARMS, UTERUS):
        members = [a for a in assays if a.arm == arm]
        donors = {a.donor for a in members}
        specimens = {a.specimen for a in members}
        groups[arm] = {
            "assays": len(members), "specimens": len(specimens), "donors": len(donors),
            "specimens_with_native_pathology_row": len(specimens & pathology.keys()),
            "donor_hormonal_medication": dict(Counter(medication[d] for d in donors)),
            "donors_seen_in_multiple_runs": sum(
                len({a.run for a in members if a.donor == d}) > 1 for d in donors),
        }
    runs = {
        run: {arm: len({a.donor for a in assays if a.arm == arm and a.run == run})
              for arm in (*ARMS, UTERUS)}
        for run in sorted({a.run for a in assays}, key=int)
    }
    return {
        "native_participants": len(tables["patient"]), "native_assays_all_tissues": len(tables["sample_analysis_result"]),
        "tissue_label_capacity": groups, "donors_by_run_within_arm": runs,
        "excluded_peritoneal_assays_by_recruitment_and_confirmed_status": dict(excluded),
        "case_lesion_and_peritoneum": paired_capacity(assays, ARMS[0], ARMS[1], medication),
        "control_peritoneum_and_uterus": paired_capacity(assays, ARMS[2], UTERUS, medication),
        "all_native_pathology_rows_by_specimen_tissue": dict(tissue_histology),
        "data_origin_codes": sorted(r["value"] for r in codes.values() if r["code_type"] == "data_origin"),
        "lesion_histopathology_codebook_entries": sum(
            r["code_type"] == "lesion_histopathology" for r in codes.values()),
    }


def catalog(head_file):
    """Inspect only COPY table definitions in the PostgreSQL header."""
    from inflate64 import Inflater

    head = Inflater().inflate(pinned(head_file, HEAD_SHA256)[48:])
    schemas = dict(re.findall(rb"COPY public\.(\w+) \(([^)]+)\) FROM stdin;", head))
    allowed = (b"patient", b"patient_sample", b"patient_survey", b"pathology_report_endometrium",
               b"sample", b"sample_analysis_result")
    links = re.findall(
        rb"ALTER TABLE ONLY public\.(\w+)\s+ADD CONSTRAINT \w+ FOREIGN KEY \((\w+)\) "
        rb"REFERENCES public\.(\w+)\((\w+)\);", head)
    metadata_links = sorted(
        f"{t.decode()}.{k.decode()} -> {u.decode()}.{v.decode()}"
        for t, k, u, v in links if t in allowed and (u in allowed or u == b"code_value"))
    return {
        "clinical_foreign_keys": metadata_links,
        "copy_table_names": sorted(k.decode() for k in schemas),
        "clinical_join_columns": {k.decode(): schemas[k].decode().split(", ") for k in allowed},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--head", type=Path, required=True)
    parser.add_argument("--tail", type=Path, required=True)
    args = parser.parse_args()
    report = summarize(read_metadata(args.head, args.tail))
    report["release_catalog"] = catalog(args.head)
    report["input"] = {
        "public_release": "https://doi.org/10.6084/m9.figshare.12195975.v1",
        "file_url": "https://ndownloader.figshare.com/files/22809185",
        "head": {"inclusive_bytes": [0, 65535], "sha256": HEAD_SHA256},
        "tail": {"inclusive_bytes": [2955439946, 2955839945], "sha256": TAIL_SHA256},
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
