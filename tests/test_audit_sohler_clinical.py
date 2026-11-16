"""Clinical labels cannot be assigned by row order or tissue similarity."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_sohler_clinical import TISSUES, compare_metadata


def mock_row(number, tissue):
    return {
        "Characteristics[Individual]": f"patient {number}",
        "Factor Value[INDIVIDUAL]": f"patient {number}",
        "Characteristics[DiseaseState]": "endometriosis",
        "Characteristics[OrganismPart]": tissue,
    }


def test_noncontiguous_explicit_ids_and_missing_peritoneum_only_donor():
    lesion, peritoneum, uterus = TISSUES
    rows = [
        mock_row(105, lesion), mock_row(107, uterus),
        mock_row(105, peritoneum), mock_row(113, peritoneum),
    ]
    clinical = [(107, "III", "Proliferative", "--", "+"), (105, "I", "Secretory", "+", "--")]
    result = compare_metadata(rows, clinical)
    assert result["checked_tissue_flags"] == 4
    assert result["donors_without_clinical_rows"] == 1
    assert result["lesion_peritoneum_pair_phases"] == {"Secretory": 1}
    invalid_tables = [
        [(104, "III", "Proliferative", "--", "+"), clinical[1]],  # Unknown ID
        [clinical[1], clinical[1]],  # Duplicate ID
        [clinical[0]],  # Missing ID
        [(105, "III", "Proliferative", "--", "+"), clinical[1]],  # Tissue flags disagree
    ]
    for wrong in invalid_tables:
        with pytest.raises(ValueError):
            compare_metadata(rows, wrong)
