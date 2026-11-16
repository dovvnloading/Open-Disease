"""All clinical and array keys below are fictitious. No expression values are needed."""
from dataclasses import replace
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from normal_peritoneum_design import Pair, candidate_pairs, frozen_pairs_json, verified_pairs
from test_audit_peritoneal_controls import fixture


def test_clinical_selection_keeps_only_assays_within_the_unique_run_and_specimen():
    t = fixture()
    # A repeated assay of the uterine specimen inside the paired run and an
    # unrelated uterine assay outside it. Both retain their original native key.
    uterus = next(a for a in t["sample_analysis_result"] if a["sample_id"] == "specimen-4")
    t["sample_analysis_result"].extend([dict(uterus, id="tech-rep"), dict(uterus, id="other", analysis_run_id="20")])
    # Unexpected endometriosis among control recruits must fail even if both
    # organs are available in the same run.
    t["patient_sample"].append(dict(id="alien", sample_id="extra", patient_id="c", code_tissue_type_id="E"))
    t["sample_analysis_result"].append(dict(uterus, id="not-eligible", sample_id="extra"))
    selected, result = candidate_pairs(t)
    assert result["candidate_pairs"] == result["healthy_donors_with_both_tissues"] == 1
    p, = selected
    assert p.donor == "b" and p.sensitivity_phase == "secretory"
    assert p.uterus_assays == ("assay-4-0", "tech-rep")
    assert result["candidate_assays_by_tissue"] == {"peritoneum": 1, "uterus": 2}

    # Order in the source is immaterial to the private freeze.
    for rows in t.values():
        rows.reverse()
    assert frozen_pairs_json(candidate_pairs(t)[0]) == frozen_pairs_json(selected)


@pytest.mark.parametrize("ambiguity", ["two shared runs", "two uterine specimens", "no shared run"])
def test_design_ambiguity_is_resolved_before_chip_classification(ambiguity):
    t = fixture()
    p, u = [a for a in t["sample_analysis_result"] if a["sample_id"] in {"specimen-3", "specimen-4"}]
    if ambiguity == "two shared runs":
        for a in (p, u):
            t["sample_analysis_result"].append(dict(a, id="copy-" + a["id"], analysis_run_id="20"))
    elif ambiguity == "two uterine specimens":
        t["patient_sample"].append(dict(id="new", sample_id="extra", patient_id="b", code_tissue_type_id="E"))
        t["sample_analysis_result"].append(dict(u, id="new-assay", sample_id="extra"))
    else:
        u["analysis_run_id"] = "20"
    selected, report = candidate_pairs(t)
    assert not selected
    assert sum(report["excluded_donors_by_design"].values()) == 1


def test_verified_pairs_preserve_good_repeats_and_lock_blinded_exclusions():
    p = Pair("a", "r", "specp", "specu", ("p1", "p2", "p3"), ("u1",), "secretory", "stored_false")
    q = Pair("b", "r", "specq", "specv", ("q1",), ("v1",), "proliferative", "unknown")
    s = {"p1": "WG6", "p2": "WG6", "p3": "repeated_addresses", "u1": "WG6", "q1": "missing", "v1": "HT12_R1"}
    kept, report = verified_pairs([p, q], s)
    assert kept == [replace(p, peritoneum_assays=("p1", "p2"))]
    assert report["excluded_donors_by_chip"] == {"no_verified_peritoneum_or_uterus": 1}
    assert report["retained_assays_by_tissue"] == {"peritoneum": 2, "uterus": 1}
    assert report["retained_sensitivity_phase_counts"] == {"secretory": 1}
    assert not report["passes_pair_count_gate"]
    assert q.sensitivity_phase is None
    for wrong in [dict(s, other="WG6"), dict(s, p1="unchecked"), {k: v for k, v in s.items() if k != "p1"}]:
        with pytest.raises(ValueError, match="classify exactly"):
            verified_pairs([p, q], wrong)
    with pytest.raises(ValueError, match="Repeated donor or assay"):
        verified_pairs([p, p], s)


def test_six_donor_resolution_gate_does_not_use_replicate_count():
    pairs = [Pair(str(i), "r", "p"+str(i), "u"+str(i), ("a"+str(i),), ("b"+str(i),), None, "stored_false")
             for i in range(6)]
    statuses = {a: "WG6" for pair in pairs for a in pair.peritoneum_assays + pair.uterus_assays}
    assert verified_pairs(pairs, statuses)[1]["passes_pair_count_gate"]
    statuses["b2"] = "ambiguous_manifest"
    assert not verified_pairs(pairs, statuses)[1]["passes_pair_count_gate"]


def test_missing_diagnosis_or_recruitment_excludes_otherwise_paired_control():
    for key in ("code_patient_category_id", "code_confirmed_disease_status_id"):
        tables = fixture()
        control = next(p for p in tables["patient"] if p["id"] == "b")
        control[key] = r"\N"
        selected, report = candidate_pairs(tables)
        assert not selected and report["healthy_donors_with_either_tissue"] == 0
