"""Compare Jiang's published E-MTAB-694 clinical table with its named SDRF donors."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import io
import json
from pathlib import Path

import openpyxl

from audit_sohler_peritoneum import SDRF_SHA256, TISSUES, summarize


JIANG_SHA256 = "7712eff6b91b7daab822b8c744b35c8837701398ba07c41e770b30d9225894ee"
HEADER = (
    "Patient I.D.", "Age [yr]", "rAFS stage", "Menstrual phase",
    "Peritoneal Lesion", "Endometrium",
)


def read_pinned(path, size, digest):
    with path.open("rb") as infile:
        content = infile.read(size + 1)
    if len(content) != size or hashlib.sha256(content).hexdigest() != digest:
        raise ValueError("Unexpected source size or SHA-256")
    return content


def compare_metadata(rows, clinical):
    """Use only the published study-local ID. Tissue flags check the join afterward."""
    summarize(rows)  # Verify the original donor/tissue and case-only source contract.
    tissue_by_donor = defaultdict(set)
    for row in rows:
        tissue_by_donor[row["Characteristics[Individual]"]].add(row["Characteristics[OrganismPart]"])
    lesion, peritoneum, uterus = TISSUES
    required = {donor for donor, tissues in tissue_by_donor.items() if tissues & {lesion, uterus}}
    seen = set()
    phases = Counter()
    paired_phases = Counter()
    stages = Counter()
    for number, stage, phase, lesion_flag, uterus_flag in clinical:
        if type(number) is not int or number <= 0:
            raise ValueError("Expected positive integer E-MTAB-694 Patient I.D.")
        donor = f"patient {number}"
        if donor in seen or donor not in required:
            raise ValueError("Duplicate or unmatched E-MTAB-694 Patient I.D.")
        seen.add(donor)
        tissues = tissue_by_donor[donor]
        flags = (lesion_flag, uterus_flag)
        if flags != tuple("+" if tissue in tissues else "--" for tissue in (lesion, uterus)):
            raise ValueError("Tissue-availability disagreement for a published donor key")
        if phase not in ("Secretory", "Proliferative") or stage not in ("I", "II", "III", "IV"):
            raise ValueError("Unexpected phase or rAFS category")
        phases[phase] += 1
        stages[stage] += 1
        if {lesion, peritoneum} <= tissues:
            paired_phases[phase] += 1
    if seen != required:
        raise ValueError("Incomplete tissue-donor coverage in clinical table")
    return {
        "clinical_donors": len(seen),
        "donors_without_clinical_rows": len(tissue_by_donor) - len(seen),
        "checked_tissue_flags": 2 * len(seen),
        "clinical_phases": dict(phases),
        "clinical_rAFS_stages": dict(stages),
        "lesion_peritoneum_pair_phases": dict(paired_phases),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdrf", type=Path)
    parser.add_argument("clinical_xlsx", type=Path)
    args = parser.parse_args()
    raw = read_pinned(args.sdrf, 25_579, SDRF_SHA256)
    rows = list(csv.DictReader(io.StringIO(raw.decode("ascii")), delimiter="\t"))
    workbook = openpyxl.load_workbook(
        io.BytesIO(read_pinned(args.clinical_xlsx, 13_416, JIANG_SHA256)),
        data_only=True, read_only=True,
    )
    try:
        table = list(workbook["Sheet1"].values)
        header = tuple(" ".join(label.split()) for label in table[1][:6])
        if "E-MTAB-694" not in table[0][0] or header != HEADER:
            raise ValueError("Unexpected source study or header")
        clinical = []
        for row in table:
            if type(row[0]) is int:
                stage, phase, lesion_flag, uterus_flag = (value.strip() for value in row[2:6])
                clinical.append((row[0], stage, phase, lesion_flag, uterus_flag))
    finally:
        workbook.close()
    print(json.dumps(compare_metadata(rows, clinical), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
