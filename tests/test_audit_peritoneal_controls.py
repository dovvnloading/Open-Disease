"""Fictitious native records with repeats, a clinically discordant recruit and uterus-only histology."""
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_peritoneal_controls import Assay, UTERUS, paired_capacity, summarize


def fixture():
    labels = {
        "P": "Peritoneum tissue", "H": "Peritoneal flush", "E": "Endometrium tissue",
        "PU1": "Red peritoneal endometriosis tissue, Replicate 1",
        "VA2": "White peritoneal endometriosis tissue, Replicate 2",
        "patient": "patient", "control": "control", "disease": "endometriosis",
        "healthy": "healthy", "unexpected": "unexpected endometriosis", "phase": "Secretory",
    }
    tables = {
        "code_value": [{"id": k, "value": k, "title_en": title, "code_type": "fixture"}
                       for k, title in labels.items()],
        "patient": [], "patient_sample": [], "sample_analysis_result": [], "patient_survey": [],
        "pathology_report_endometrium": [],
    }
    for donor, group, disease in [("a", "patient", "disease"), ("b", "control", "healthy"),
                                  ("c", "control", "unexpected")]:
        tables["patient"].append({"id": donor, "code_patient_category_id": group,
                                  "code_confirmed_disease_status_id": disease})
    tables["patient_survey"].append({"patient_id": "b", "hormonal_medication": "f"})
    for i, (donor, tissue, runs) in enumerate([("a", "PU1", ["10", "10"]), ("a", "P", ["11"]),
                                              ("a", "VA2", ["11"]), ("b", "P", ["10"]),
                                              ("b", "E", ["10"]), ("c", "P", ["10"]),
                                              ("c", "PU1", ["10"]), ("a", "H", ["10"])]):
        sample = f"specimen-{i}"
        tables["patient_sample"].append({"id": f"row-{i}", "sample_id": sample, "patient_id": donor,
                                         "code_tissue_type_id": tissue})
        for j, run in enumerate(runs):
            tables["sample_analysis_result"].append({"id": f"assay-{i}-{j}", "sample_id": sample,
                                                       "analysis_run_id": run,
                                                       "result_data_type": "microarray"})
        if tissue == "E":
            tables["pathology_report_endometrium"].append({"patient_sample_id": f"row-{i}",
                                                           "code_histology_class_id": "phase"})
    return tables


def test_repeats_clinical_vs_recruitment_arms_and_site_specific_histology():
    report = summarize(fixture())
    groups = report["tissue_label_capacity"]
    lesion = groups["lesion"]
    assert (lesion["assays"], lesion["specimens"], lesion["donors"]) == (3, 2, 1)
    assert lesion["donors_seen_in_multiple_runs"] == 1
    assert lesion["donor_hormonal_medication"] == {"unknown": 1}
    assert groups["case_peritoneum"]["assays"] == groups["control_peritoneum"]["assays"] == 1
    assert report["case_lesion_and_peritoneum"]["donors_with_a_shared_run"] == 1
    assert report["case_lesion_and_peritoneum"]["donors_with_one_assay_each_and_shared_run"] == 0
    assert report["case_lesion_and_peritoneum"]["donors_with_one_common_run_one_specimen_each_within_it"] == 1
    assert sum(report["excluded_peritoneal_assays_by_recruitment_and_confirmed_status"].values()) == 2
    assert report["all_native_pathology_rows_by_specimen_tissue"] == {"Endometrium tissue": 1}
    assert groups["control_peritoneum"]["specimens_with_native_pathology_row"] == 0
    phases = report["control_peritoneum_and_uterus"]["stored_false_shared_run_with_consistent_matched_uterine_phase"]
    assert phases == {"secretory": 1}


def test_cross_run_only_conflicting_matched_phase_and_unknown_medication():
    rows = [
        Assay("P", "cross", "p1", "9", None), Assay(UTERUS, "cross", "e1", "10", "Secretory"),
        Assay("P", "mixed", "p2", "10", None), Assay(UTERUS, "mixed", "e2", "10", "Proliferative"),
        Assay(UTERUS, "mixed", "e3", "10", None),
        Assay("P", "missing", "p3", "11", None), Assay(UTERUS, "missing", "e4", "11", "Secretory"),
    ]
    result = paired_capacity(rows, "P", UTERUS, {"cross": "stored_false", "mixed": "stored_false",
                                                  "missing": "unknown"})
    assert result["donors_with_both_labels"] == 3
    assert result["donors_with_a_shared_run"] == 2
    assert result["donors_with_one_common_run_one_specimen_each_within_it"] == 1
    assert result["stored_false_shared_run_with_consistent_matched_uterine_phase"] == {}


def test_unknown_run_cannot_manufacture_pair():
    tables = fixture()
    tables["sample_analysis_result"][0]["analysis_run_id"] = r"\N"
    with pytest.raises(ValueError, match="Missing native run"):
        summarize(tables)
