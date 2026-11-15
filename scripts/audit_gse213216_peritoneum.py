"""Summarize peritoneal comparators in Fonseca Supplementary Tables S1 and S2."""

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from openpyxl import load_workbook


EXPECTED_SHA256 = "9a05e5fc20a726dca2e685ad671b18bada4545a1f353bc0fa7c4abe8fcf9087b"


def audit(path):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise ValueError("Workbook differs from the audited version")
    workbook = load_workbook(path, data_only=True, read_only=True)
    try:
        diagnoses = {
            row[0]: row[1] for row in workbook["S2"].values
            if isinstance(row[0], int)
        }
        controls = Counter()
        positive_sites = defaultdict(set)
        negative_sites = defaultdict(set)
        negative_count = 0
        control_peritoneum = 0
        donor = None
        for row in workbook["S1"].iter_rows(min_row=5, values_only=True):
            key, status, site, surgical_type, _, _, _, pathology = row
            if site is None:
                continue
            if key is not None:
                donor = key
                if status != diagnoses[donor]:
                    raise ValueError("Specimen and participant status disagree")
            is_peritoneal = surgical_type.startswith("Peritoneal")
            if diagnoses[donor] == "No":
                controls[surgical_type] += 1
                control_peritoneum += is_peritoneal
            if is_peritoneal and pathology == "Endometriosis":
                positive_sites[donor].add(site)
            # Paper identifies these four S1 NA-dagger specimens as negative on histology.
            if is_peritoneal and pathology == "NA†":
                negative_sites[donor].add(site)
                negative_count += 1
    finally:
        workbook.close()
    paired = positive_sites.keys() & negative_sites.keys()
    return {
        "workbook_sha256": digest,
        "participants_without_confirmed_endometriosis": sum(s == "No" for s in diagnoses.values()),
        "surgical_specimens_from_these_participants": dict(controls),
        "peritoneal_specimens_from_these_participants": control_peritoneum,
        "peritoneal_specimens_negative_on_pathology": negative_count,
        "distinct_donors_of_negative_peritoneal_specimens": len(negative_sites),
        "negative_peritoneal_donors_with_confirmed_disease": sum(
            diagnoses[donor] == "Yes" for donor in negative_sites
        ),
        "donors_with_positive_and_negative_peritoneal_specimens": len(paired),
        "paired_donors_with_any_identical_anatomic_location_label": sum(
            bool(positive_sites[donor] & negative_sites[donor]) for donor in paired
        ),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.workbook), indent=2, ensure_ascii=False))
