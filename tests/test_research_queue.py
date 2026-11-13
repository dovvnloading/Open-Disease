"""Task orchestration tests use only toy records; no data or network."""
from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from research_queue import load_tasks, new_run, ready_tasks  # noqa: E402


def task(ident="toy-one", deps=None):
    return {
        "schema_version": 1, "id": ident, "title": "Toy task", "program": "test",
        "priority": "P0", "status": "ready",
        "issue_url": "https://github.com/dovvnloading/Open-Disease/issues/11",
        "question": "Can we verify this toy?", "scope": "Synthetic only",
        "dependencies": [] if deps is None else deps,
        "sources": ["https://example.org/toy"], "steps": ["Run toy"],
        "deliverables": ["Toy report"], "acceptance": ["No real data"],
        "budget": {"max_runtime_minutes": 1, "max_ram_gb": 1,
                   "network": "none", "max_download_mb": 0},
    }


def write(folder, data):
    (folder / (data["id"] + ".json")).write_text(json.dumps(data), encoding="utf-8")


def test_ready_depends_on_completed_tasks(tmp_path):
    write(tmp_path, task())
    write(tmp_path, task("toy-two", ["toy-one"]))
    tasks = load_tasks(tmp_path)
    assert [x["id"] for x in ready_tasks(tasks)] == ["toy-one"]
    first = task()
    first["status"] = "done"
    write(tmp_path, first)
    assert [x["id"] for x in ready_tasks(load_tasks(tmp_path))] == ["toy-two"]


def test_reject_dependency_cycle_and_unknown_task(tmp_path):
    write(tmp_path, task("toy-one", ["toy-two"]))
    write(tmp_path, task("toy-two", ["toy-one"]))
    with pytest.raises(ValueError, match="cycle"):
        load_tasks(tmp_path)
    write(tmp_path, task("toy-two", ["nonexistent"]))
    with pytest.raises(ValueError, match="unknown dependency"):
        load_tasks(tmp_path)


def test_reject_invalid_urls_and_offline_download(tmp_path):
    bad = task()
    bad["sources"] = ["http://invalid.example"]
    write(tmp_path, bad)
    with pytest.raises(ValueError, match="HTTPS"):
        load_tasks(tmp_path)
    bad = deepcopy(task())
    bad["budget"]["max_download_mb"] = 1
    write(tmp_path, bad)
    with pytest.raises(ValueError, match="zero download"):
        load_tasks(tmp_path)


def test_new_run_requires_ready_task_and_ignored_output(tmp_path):
    record = task()
    tasks = {record["id"]: record}
    with pytest.raises(ValueError, match="under ignored results"):
        new_run(tasks, "toy-one", "agent", tmp_path / "unsafe.json")
    with pytest.raises(ValueError, match="not ready"):
        new_run(tasks, "unknown", "agent", tmp_path / "unsafe.json")


def test_new_run_creates_unpublished_provenance_stub(tmp_path, monkeypatch):
    import research_queue
    monkeypatch.setattr(research_queue, "ROOT", tmp_path)
    record = task()
    tasks = {record["id"]: record}
    destination = tmp_path / "results" / "test-research-queue" / "run.json"
    new_run(tasks, "toy-one", "test-agent", destination)
    saved = json.loads(destination.read_text(encoding="utf-8"))
    assert saved["status"] == "in_progress"
    assert saved["clinical_use"] is False
    with pytest.raises(FileExistsError):
        new_run(tasks, "toy-one", "test-agent", destination)
