"""Exercise the published aggregate workflow and its input checks."""

import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import matrisome_ranking_uncertainty as study


def test_reproduce_published_sensitivity_results():
    report = study.analyse(study.read_pinned(study.BENCHMARK), study.read_pinned(study.SCAN_DESIGN))
    expected = study.EVIDENCE / "matrisome-ranking-uncertainty" / "results.json"
    assert report == json.loads(expected.read_text())
    assert report["us_phase_sensitivity"]["equal_phase_fraction"] == pytest.approx(143 / 270)
    assert report["array_label_tests"]["GSE6364"]["assignments"] == 462 * 84 * 24310
    options = report["rome_scan_day_sensitivity"]["compatible_within_day_tests"]
    assert [(row["possible_within_day_wins"], row["extreme"]) for row in options] == [
        (11, 62), (12, 34), (13, 16), (14, 6), (15, 2)
    ]


def test_refuse_cross_group_ties_or_inconsistent_pairs():
    row = {"cases": 2, "comparators": 3, "higher": 3, "pairs": 6, "ties": 0}
    with pytest.raises(ValueError, match="ties"):
        study.stratum(dict(row, ties=1))
    with pytest.raises(ValueError, match="pair counts"):
        study.stratum(dict(row, pairs=5))


def test_refuse_changed_scan_margins_and_changed_source_bytes(tmp_path, monkeypatch):
    design = study.read_pinned(study.SCAN_DESIGN)
    with pytest.raises(ValueError, match="margins"):
        study.rome_scan_sensitivity((7, 6, 40), design)
    monkeypatch.setattr(study, "EVIDENCE", tmp_path)
    source = tmp_path / study.BENCHMARK
    source.parent.mkdir()
    source.write_bytes(b"{}")
    with pytest.raises(ValueError, match="audited release"):
        study.read_pinned(study.BENCHMARK)
