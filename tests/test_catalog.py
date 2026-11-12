import json
from pathlib import Path

import pytest

from scripts.repo_safety import check_paths
from scripts.verify_catalog import audit

ROOT = Path(__file__).parents[1]


def test_canonical_catalog_is_metadata_only():
    assert audit(ROOT / "data/sources.json") == {
        "sources": 6, "status": "metadata_only", "analyzed": 0
    }


def test_catalog_cannot_pretend_human_data_were_analyzed(tmp_path):
    manifest = json.loads((ROOT / "data/sources.json").read_text())
    manifest["biological_results_generated"] = True
    path = tmp_path / "sources.json"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="biological_results_generated"):
        audit(path)


def test_guard_rejects_unsafe_tracked_paths():
    names = ["README.md", "data/raw/participant.csv", "results/diagnosis.csv",
             "misc/sample.h5ad", ".env", "docs/PROTOCOL.md"]
    assert check_paths(names) == names[1:5]
