"""Count case-side assay and donor support in the E-MTAB-694 sample manifest."""
import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import io
import json
from pathlib import Path


SDRF_SHA256 = "29ef1785d73103914ee0ea6b9fd87841dd42d687d01fe410c19e8d0a2589b393"
TISSUES = ("peritoneal endometriosis lesion", "peritoneum", "endometrium")


def summarize(rows):
    donors = defaultdict(set)
    assays = Counter()
    seen = set()
    for row in rows:
        donor = row["Characteristics[Individual]"].strip()
        tissue = row["Characteristics[OrganismPart]"].strip()
        if not donor or donor != row["Factor Value[INDIVIDUAL]"].strip():
            raise ValueError("Missing or inconsistent source donor annotation")
        if row["Characteristics[DiseaseState]"].strip() != "endometriosis" or tissue not in TISSUES:
            raise ValueError("Unexpected cohort annotation")
        if (donor, tissue) in seen:
            raise ValueError("Repeated tissue within donor requires separate assay adjudication")
        seen.add((donor, tissue))
        assays[tissue] += 1
        donors[tissue].add(donor)
    lesion, peritoneum, uterus = (donors[tissue] for tissue in TISSUES)
    return {
        "arrays": sum(assays.values()),
        "arrays_by_tissue": dict(assays),
        "case_donors": len(lesion | peritoneum | uterus),
        "noncase_peritoneum_donors": 0,
        "lesion_peritoneum_pairs": len(lesion & peritoneum),
        "lesion_uterus_pairs": len(lesion & uterus),
        "uterus_peritoneum_pairs": len(uterus & peritoneum),
        "all_three_tissues": len(lesion & peritoneum & uterus),
        "peritoneum_only_donors": len(peritoneum - lesion - uterus),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sdrf", type=Path, help="Local unmodified public E-MTAB-694.sdrf.txt")
    path = parser.parse_args().sdrf
    with path.open("rb") as handle:
        raw = handle.read(150_001)
    if len(raw) > 150_000 or sha256(raw).hexdigest() != SDRF_SHA256:
        raise ValueError("Unexpected SDRF version or size")
    rows = csv.DictReader(io.StringIO(raw.decode("ascii")), delimiter="\t")
    print(json.dumps(summarize(rows), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
