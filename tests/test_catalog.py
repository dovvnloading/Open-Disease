import json
from pathlib import Path

import pytest

from scripts.repo_safety import check_paths
from scripts.verify_catalog import audit

ROOT = Path(__file__).parents[1]


def test_canonical_catalog_links_executed_analysis():
    catalog = ROOT / "data/sources.json"
    assert audit(catalog) == {"sources": 12, "status": "exploratory_computation", "analyzed": 6}
    sources = json.loads(catalog.read_text())["sources"]
    assert [s["id"] for s in sources if s["human_data_used_by_open_disease"]] == [
        "GSE266265", "GSE25628", "GSE141549", "GSE120103", "GSE315857", "Fung2018-S20"
    ]


def sample_catalog(tmp_path, *, used=False, claim=False, report=None):
    folder = tmp_path / "data"
    folder.mkdir()
    sample = {
        "id": "toy", "url": "https://example.org/toy", "kind": "microarray",
        "access": "public", "role": "example", "cautions": ["synthetic"],
        "human_data_used_by_open_disease": used,
    }
    if report is not None:
        sample["analysis_report"] = report
    data = {
        "schema_version": "2026-11-12", "sources": [sample],
        "automatic_downloads_enabled": False, "catalog_is_not_redistribution_permission": True,
        "biological_results_generated": claim,
    }
    path = folder / "sources.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_legacy_metadata_only_catalog_is_accepted(tmp_path):
    assert audit(sample_catalog(tmp_path)) == {"sources": 1, "status": "metadata_only", "analyzed": 0}


@pytest.mark.parametrize(("used", "claim", "report", "message"), [
    (False, True, None, "biological_results_generated"),
    (True, False, "docs/research-notes/present.md", "biological_results_generated"),
    (True, True, None, "analysis_report"),
    (True, True, "docs/research-notes/missing.md", "analysis_report"),
    (True, True, "docs/research-notes/../../outside.md", "analysis_report"),
    (False, False, "docs/research-notes/present.md", "analysis_report"),
    ("true", True, None, "boolean"),
])
def test_catalog_rejects_inconsistent_or_unsupported_claims(tmp_path, used, claim, report, message):
    notes = tmp_path / "docs/research-notes"
    notes.mkdir(parents=True)
    (notes / "present.md").write_text("Documented synthetic result")
    (tmp_path / "outside.md").write_text("Outside the report directory")
    path = sample_catalog(tmp_path, used=used, claim=claim, report=report)
    with pytest.raises(ValueError, match=message):
        audit(path)


def test_guard_rejects_unsafe_tracked_paths():
    names = ["README.md", "data/raw/participant.csv", "results/diagnosis.csv",
             "misc/sample.h5ad", ".env", "docs/PROTOCOL.md"]
    assert check_paths(names) == names[1:5]
