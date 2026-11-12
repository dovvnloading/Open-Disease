import json
from pathlib import Path

import pytest

from scripts.repo_safety import check_paths
from scripts.verify_catalog import audit

ROOT = Path(__file__).parents[1]


def test_source_manifest_is_metadata_only():
    status = audit(ROOT / "catalog/sources.json")
    assert status == {"sources": 3, "status": "metadata_only", "analyzed": 0}


def test_manifest_cannot_pretend_human_data_were_analyzed(tmp_path):
    manifest = json.loads((ROOT / "catalog/sources.json").read_text())
    manifest["clinical_data_ingested"] = True
    path = tmp_path / "sources.json"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="ingestion"):
        audit(path)


def test_safety_guard_catches_restricted_file_paths():
    paths = ["README.md", "data/raw/participant.csv", "results/diagnosis.csv",
             "misc/sample.h5ad", ".env", "docs/PROTOCOL.md"]
    assert check_paths(paths) == paths[1:5]
