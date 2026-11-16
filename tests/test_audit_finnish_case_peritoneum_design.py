from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_finnish_case_peritoneum_design import CasePair, case_pairs, chip_gate, donor_uterine_phase, noncase_reference
from audit_peritoneal_controls import Assay


def test_shared_run_specimens_distinguished_from_technical_replicates():
    rows = [
        Assay("lesion", "ok", "l", "1", None, "l1"),
        Assay("lesion", "ok", "l", "1", None, "l2"),
        Assay("lesion", "ok", "extra", "2", None, "l3"),
        Assay("case_peritoneum", "ok", "p", "1", None, "p1"),
        Assay("lesion", "ambiguous", "l", "1", None, "x1"),
        Assay("lesion", "ambiguous", "ll", "1", None, "x2"),
        Assay("case_peritoneum", "ambiguous", "p", "1", None, "x3"),
        Assay("lesion", "unpaired", "l", "2", None, "y1"),
        Assay("case_peritoneum", "unpaired", "p", "1", None, "y2"),
    ]
    pairs, report = case_pairs(rows)
    assert pairs == [CasePair("ok", "1", ("l1", "l2"), ("p1",), "l", "p")]
    assert report["excluded_donors"] == {"multiple_specimens_in_shared_run": 1, "no_shared_run": 1}


def test_incomplete_or_discordant_chip_pairs_do_not_count():
    pairs = [CasePair("a", "1", ("l1", "l2"), ("p1",), "l", "p"),
             CasePair("b", "2", ("l3",), ("p2",), "other l", "other p")]
    statuses = {"l1": "repeated_addresses", "l2": "WG6", "p1": "WG6", "l3": "WG6", "p2": "HT12_R1"}
    with pytest.raises(ValueError, match="classification"):
        chip_gate(pairs, {"l1": "WG6"}, {}, {})
    retained, report = chip_gate(pairs, statuses, {"a": "stored_true", "b": "stored_false"}, {})
    assert retained == [CasePair("a", "1", ("l2",), ("p1",), "l", "p")]
    assert report["donors_passing_wg6_both_arms"] == 1
    assert report["retained_wg6_assays_by_tissue"] == {"lesion": 1, "case_peritoneum": 1}
    assert report["retained_donors_by_survey_hormone_flag"] == {"stored_true": 1}
    assert report["retained_uterine_histology_by_donor_unverified_as_same_operation"] == {"no_uterine_report": 1}


def test_discordant_donor_uterine_histology_is_never_assigned_to_a_case_biopsy():
    tables = {
        "code_value": [{"id": "e", "value": "E"}, {"id": "p", "value": "P"},
                       {"id": "s", "title_en": "Secretory"}, {"id": "f", "title_en": "Proliferative"}],
        "patient_sample": [
            {"id": "one", "patient_id": "a", "code_tissue_type_id": "e"},
            {"id": "two", "patient_id": "a", "code_tissue_type_id": "e"},
            {"id": "three", "patient_id": "b", "code_tissue_type_id": "p"},
        ],
        "pathology_report_endometrium": [
            {"patient_sample_id": "one", "code_histology_class_id": "s"},
            {"patient_sample_id": "two", "code_histology_class_id": "f"},
            {"patient_sample_id": "three", "code_histology_class_id": "s"},
        ],
    }
    assert donor_uterine_phase(tables) == {"a": "mixed_or_inconsistent"}


def test_noncase_reference_deduplicates_specimens_and_donors_across_supported_runs():
    pairs = [CasePair("case" + run, run, (), (), "l" + run, "p" + run) for run in ("1", "2", "3")]
    rows = [
        Assay("control_peritoneum", "prior", "spA", "1", None, "a1"),
        Assay("control_peritoneum", "prior", "spA", "1", None, "a2"),
        Assay("control_peritoneum", "prior", "spA", "2", None, "a3"),
        Assay("control_peritoneum", "new", "spB", "1", None, "b1"),
        Assay("control_peritoneum", "other_chip", "spC", "3", None, "c1"),
        Assay("control_peritoneum", "outside_runs", "spD", "4", None, "d1"),
    ]
    statuses = {a.assay_id: "WG6" for a in rows}
    statuses["c1"] = "HT12_R1"
    flags = {p.donor: "stored_true" for p in pairs} | {"prior": "stored_false", "new": "unknown"}
    summary = noncase_reference(pairs, rows, statuses, flags, {"prior", "unused"})
    assert [r["noncase_peritoneum_donors"] for r in summary["run_summary"].values()] == [2, 1, 0]
    assert summary["case_pairs_in_runs_with_wg6_noncase_peritoneum"] == 2
    assert summary["case_pairs_with_same_run_noncase_sharing_reported_nonmissing_flag"] == {}
    assert summary["noncase_donors_in_multiple_case_runs"] == 1
    assert summary["distinct_noncase_donors"] == summary["distinct_noncase_specimens"] == 2
    assert summary["retained_noncase_wg6_assays"] == 4
    assert summary["overlap_with_prior_paired_organ_control_donors"] == 1
    assert summary["noncase_donors_outside_prior_paired_organ_controls"] == 1
    assert summary["distinct_noncase_donors_by_survey_flag"] == {"stored_false": 1, "unknown": 1}
