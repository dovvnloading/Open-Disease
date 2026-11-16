"""Synthetic source rows distinguish matched lesions from any matched tissue."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_sohler_peritoneum import TISSUES, summarize


def source_row(donor, tissue):
    return {
        "Characteristics[Individual]": donor,
        "Factor Value[INDIVIDUAL]": donor,
        "Characteristics[OrganismPart]": tissue,
        "Characteristics[DiseaseState]": "endometriosis",
    }


def test_peritoneum_matched_to_uterus_is_not_a_lesion_pair():
    lesion, peritoneum, uterus = TISSUES
    rows = [source_row("invented-a", t) for t in (lesion, peritoneum, uterus)]
    rows += [source_row("invented-b", t) for t in (peritoneum, uterus)]
    rows += [source_row("invented-c", peritoneum)]
    counts = summarize(rows)
    assert counts["case_donors"] == 3
    assert counts["lesion_peritoneum_pairs"] == 1
    assert counts["uterus_peritoneum_pairs"] == 2
    assert counts["peritoneum_only_donors"] == 1
    assert counts["noncase_peritoneum_donors"] == 0


def test_duplicate_and_discordant_donor_annotations_fail_closed():
    row = source_row("invented", TISSUES[0])
    with pytest.raises(ValueError, match="Repeated tissue"):
        summarize([row, row])
    row["Factor Value[INDIVIDUAL]"] = ""
    with pytest.raises(ValueError, match="donor annotation"):
        summarize([row])
