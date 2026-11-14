from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_heca_glandular import audit, is_case, sc_count_key


def test_only_named_phase_footnote_is_reconciled():
    assert sc_count_key("Lai", "SyntheticKeya", "SyntheticKey") == "Lai_SyntheticKey"
    assert sc_count_key("Lai", "OtherKeya", "SyntheticKey") == "Lai_OtherKeya"
    assert sc_count_key("Wang", "SyntheticKeya", "SyntheticKey") == "Wang_SyntheticKeya"
    assert sc_count_key("Fonseca", "Fonseca_SyntheticKey", "Unused") == "Fonseca_SyntheticKey"


@pytest.mark.parametrize("stage,path,expected", [(0, "P", False), ("TBC", "E", True), (2, "E + F", True)])
def test_clinical_classification(stage, path, expected):
    assert is_case({"Endometriosis_stageh": stage, "Endometrial pathologyi": path}) == expected


@pytest.mark.parametrize("stage,path", [(0, "E + F"), (2, "C"), ("unknown", "E")])
def test_inconsistent_or_unknown_disease_fails(stage, path):
    with pytest.raises(ValueError):
        is_case({"Endometriosis_stageh": stage, "Endometrial pathologyi": path})


def test_refuse_unknown_workbook_before_reading(tmp_path):
    p = tmp_path / "source.xlsx"
    p.write_bytes(b"not the source")
    with pytest.raises(ValueError, match="Unrecognized"):
        audit(p)
