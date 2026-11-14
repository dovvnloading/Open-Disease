"""Invented native database records; these tests make no network requests."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_endometdb import eligible, pinned


def tables():
    labels = ["Endometrium tissue", "control", "healthy", "patient", "endometriosis",
              "Proliferative", "Secretory", "Insufficient sample", "Peritoneum tissue"]
    data = dict(code_value=[dict(id=label, title_en=label) for label in labels],
                patient=[], patient_sample=[], patient_survey=[],
                pathology_report_endometrium=[], sample_analysis_result=[])
    for i, (category, health, medication, phase, tissue) in enumerate([
        ("control", "healthy", "f", "Proliferative", "Endometrium tissue"),
        ("patient", "endometriosis", "f", "Proliferative", "Endometrium tissue"),
        ("control", "healthy", "t", "Proliferative", "Endometrium tissue"),
        ("patient", "healthy", "f", "Proliferative", "Endometrium tissue"),
        ("patient", "endometriosis", "f", "Insufficient sample", "Endometrium tissue"),
        ("control", "healthy", "f", "Proliferative", "Peritoneum tissue"),
        ("control", "healthy", "f", "Secretory", "Endometrium tissue"),
        ("patient", "endometriosis", None, "Proliferative", "Endometrium tissue"),
    ]):
        donor = f"donor{i}"
        sample = f"specimen{i}"
        sid = f"patient-sample{i}"
        data["patient"].append(dict(id=donor, code_patient_category_id=category, code_confirmed_disease_status_id=health))
        data["patient_sample"].append(dict(id=sid, sample_id=sample, patient_id=donor, code_tissue_type_id=tissue))
        if medication is not None:
            data["patient_survey"].append(dict(patient_id=donor, hormonal_medication=medication))
        data["pathology_report_endometrium"].append(dict(patient_sample_id=sid, code_histology_class_id=phase))
        data["sample_analysis_result"].append(dict(id=f"array{i}", sample_id=sample, analysis_run_id="1"))
    data["sample_analysis_result"].append(dict(id="replicate", sample_id="specimen6", analysis_run_id="2"))
    return data


def test_native_keys_groups_washout_unknowns_and_repeat_exclusion():
    report, rows = eligible(tables())
    assert [(x["donor_id"], x["group"]) for x in rows] == [("donor0", "control"), ("donor1", "case")]
    assert report["shared_run_phase_strata"] == 1
    assert report["participants_with_no_survey"] == report["repeated_donors_excluded"] == 1
    assert report["candidate_donors"] == 3


def test_duplicate_native_key_rejected():
    data = tables()
    data["patient_survey"].append(data["patient_survey"][0])
    with pytest.raises(ValueError, match="Nonunique"):
        eligible(data)


def test_corrupted_range_rejected(tmp_path):
    f = tmp_path / "bad"
    f.write_bytes(b"invented")
    with pytest.raises(ValueError, match="pinned"):
        pinned(f, "0" * 64)
