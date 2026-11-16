"""Only invented donor, specimen, assay and signal rows are used here."""
from dataclasses import replace
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import finnish_case_peritoneum_score as experiment
import normal_peritoneum_score as ranks
from audit_finnish_case_peritoneum_design import CasePair


def design():
    return [CasePair(f"toy-{i}", run, (f"l{i}",), (f"p{i}",), f"ls{i}", f"ps{i}")
            for i, run in enumerate(["10"] * 5 + ["9", "12"])]


def fixture(monkeypatch, *, reverse=False, flat=False):
    pairs = design()
    annotation = {addr: gene for gene, (_, addr) in ranks.TARGETS.items()}
    annotation[1] = 0
    monkeypatch.setattr(ranks, "N_FEATURES", len(annotation))
    lines = [b"0\tunselected\topaque"]
    for p in pairs:
        for label, assay in (("lesion", p.lesion[0]), ("peritoneum", p.peritoneum[0])):
            high = (label == "lesion") != reverse and not flat
            for addr, gene in annotation.items():
                signal = 2
                if high:
                    signal = 2 + ranks.TARGETS.get(gene, (0,))[0]
                lines.append(f"0\t{assay}\t{addr}\t{signal}\t.\t.\t.".encode())
    return pairs, annotation, lines


@pytest.mark.parametrize(("reverse", "flat", "decision", "p_value"),
                         [(False, False, "higher_in_labelled_lesion", 2 / 128),
                          (True, False, "lower_in_labelled_lesion", 2 / 128),
                          (False, True, "unresolved", 1.0)])
def test_synthetic_rank_controls_exact_inference_and_run_deletion(monkeypatch, reverse, flat, decision, p_value):
    pairs, annotation, lines = fixture(monkeypatch, reverse=reverse, flat=flat)
    report = experiment.compare(pairs, lambda: reversed(lines), annotation)
    assert report["decision"] == decision and report["two_sided_p"] == p_value
    assert report["leave_run_10_out"] == {"donors": 2, "mean_difference_pp": report["mean_difference_pp"]}
    assert report["independent_verification_passed"]
    assert not any(p.donor in json.dumps(report) for p in pairs)


def test_run_deletion_uses_the_two_other_donors(monkeypatch):
    pairs, annotation, lines = fixture(monkeypatch)
    changed = []
    for row in lines:
        fields = row.split(b"\t")
        if fields[1] in {b"l5", b"l6"}:
            fields[3] = b"2"
        changed.append(b"\t".join(fields))
    report = experiment.compare(pairs, lambda: iter(changed), annotation)
    assert report["mean_difference_pp"] > 0
    assert report["leave_run_10_out"] == {"donors": 2, "mean_difference_pp": 0}
    assert report["positive_zero_negative"] == [5, 2, 0]


def test_outcome_execution_requires_existing_exact_private_lock(monkeypatch, tmp_path):
    pairs = design()
    aggregate = {"synthetic": True}
    aggregate_path = tmp_path / "aggregate.json"
    aggregate_path.write_text(json.dumps(aggregate))
    monkeypatch.setattr(experiment, "DESIGN", aggregate_path)
    monkeypatch.setattr(experiment, "audit_with_pairs", lambda path: (aggregate, pairs))
    with pytest.raises(ValueError, match="Missing or changed"):
        experiment.preflight(tmp_path, create=False)
    experiment.preflight(tmp_path, create=True)
    experiment.preflight(tmp_path, create=False)
    pairs[0] = replace(pairs[0], lesion=pairs[0].peritoneum, peritoneum=pairs[0].lesion)
    with pytest.raises(ValueError, match="Missing or changed"):
        experiment.preflight(tmp_path, create=False)
    aggregate_path.write_text("{}")
    with pytest.raises(ValueError, match="differs"):
        experiment.preflight(tmp_path, create=True)


def test_incomplete_case_pair_gate_precedes_signal_read():
    with pytest.raises(ValueError, match="seven distinct"):
        experiment.compare(design()[:6], lambda: pytest.fail("Signal stream opened"), {})
    pairs = design()
    pairs[0] = replace(pairs[0], lesion=pairs[1].lesion)
    with pytest.raises(ValueError, match="14 distinct"):
        experiment.compare(pairs, lambda: pytest.fail("Signal stream opened"), {})


def test_independent_score_disagreement_and_malformed_selected_signal_block_release(monkeypatch):
    pairs, annotation, lines = fixture(monkeypatch)
    monkeypatch.setattr(experiment, "verified_assay_scores", lambda *args: {})
    with pytest.raises(ValueError, match="Independent per-array"):
        experiment.compare(pairs, lambda: iter(lines), annotation)
    with pytest.raises(ValueError, match="Nonfinite"):
        experiment.compare(pairs, lambda: iter([b"0\tl0\t1\tnan\t.\t.\t."]), annotation)
