"""Task orchestration tests use only toy records; no data or network."""
from __future__ import annotations

import json
import sys
from copy import deepcopy
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import research_queue  # noqa: E402
from research_queue import load_tasks, new_run, ready_tasks  # noqa: E402


@pytest.fixture(autouse=True)
def history(tmp_path, monkeypatch):
    repo = tmp_path / "repo"
    report = repo / "docs" / "previous.md"
    report.parent.mkdir(parents=True)
    report.write_text("Prior analysis")
    path = repo / "evidence-history.json"
    path.write_text(json.dumps({
        "legacy_tasks": ["old-task", "finished-task"],
        "legacy_blocked_tasks": ["old-task"],
        "studies": [{"id": "Original", "aliases": ["alias", "GSE24680", "E-MTAB-1234", "https://doi.org/10.1234/Old"],
                     "report": "docs/previous.md"}],
    }))
    monkeypatch.setattr(research_queue, "ROOT", repo)
    monkeypatch.setattr(research_queue, "HISTORY", path)
    return path


def task(ident="toy-one", deps=None):
    source = f"https://example.org/{ident}"
    return {
        "schema_version": 2, "id": ident, "title": "Toy task", "program": "test",
        "priority": "P0", "status": "ready",
        "issue_url": "https://github.com/dovvnloading/Open-Disease/issues/15",
        "question": "Can we verify this toy?", "scope": "Synthetic only",
        "dependencies": [] if deps is None else deps,
        "novelty": {
            "prior_work": "Previous toy used different samples.",
            "contribution": "A second toy source checks the direction.",
            "falsifier": "A reversed result would contradict the hypothesis.",
            "evidence": [{"study_id": ident, "source": source, "use": "new"}],
        },
        "sources": [source], "steps": ["Run toy"],
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
    first["completion_pr"] = "https://github.com/dovvnloading/Open-Disease/pull/14"
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
    assert saved["novelty"] == record["novelty"]
    with pytest.raises(FileExistsError):
        new_run(tasks, "toy-one", "test-agent", destination)


def test_done_requires_a_traceable_pr(tmp_path):
    completed = task()
    completed["status"] = "done"
    write(tmp_path, completed)
    with pytest.raises(ValueError, match="completion_pr"):
        load_tasks(tmp_path)
    completed["completion_pr"] = "https://github.com/dovvnloading/Open-Disease/pull/14"
    write(tmp_path, completed)
    assert load_tasks(tmp_path)["toy-one"]["status"] == "done"


def complete(record):
    record["status"] = "done"
    record["completion_pr"] = "https://github.com/dovvnloading/Open-Disease/pull/14"


def test_archived_tasks_remain_readable_but_cannot_run_without_novelty(tmp_path):
    old = task("old-task")
    del old["novelty"]
    old.update(schema_version=1, status="blocked")
    write(tmp_path, old)
    assert not ready_tasks(load_tasks(tmp_path))
    with pytest.raises(ValueError, match="not ready"):
        new_run(load_tasks(tmp_path), "old-task", "agent", tmp_path / "run.json")
    complete(old)
    write(tmp_path, old)
    with pytest.raises(ValueError, match="changed-status"):
        load_tasks(tmp_path)
    old["status"] = "ready"
    write(tmp_path, old)
    with pytest.raises(ValueError, match="reactivated"):
        load_tasks(tmp_path)
    migrated = task("old-task")
    complete(migrated)
    write(tmp_path, migrated)
    assert load_tasks(tmp_path)["old-task"]["status"] == "done"
    finished = task("finished-task")
    complete(finished)
    finished["schema_version"] = 1
    write(tmp_path, finished)
    assert load_tasks(tmp_path)["finished-task"]["status"] == "done"
    finished["status"] = "blocked"
    write(tmp_path, finished)
    with pytest.raises(ValueError, match="changed-status"):
        load_tasks(tmp_path)


@pytest.mark.parametrize("status", ["ready", "blocked", "done"])
def test_new_tasks_cannot_opt_out_with_an_old_schema(tmp_path, status):
    record = task()
    complete(record)
    record.update(schema_version=1, status=status)
    write(tmp_path, record)
    with pytest.raises(ValueError, match="schema_version 2"):
        load_tasks(tmp_path)


@pytest.mark.parametrize("field,value", [
    ("prior_work", ""), ("contribution", " "), ("falsifier", None), ("evidence", []),
])
def test_novelty_cannot_be_blank_even_for_blocked_work(tmp_path, field, value):
    record = task()
    record["status"] = "blocked"
    record["novelty"][field] = value
    write(tmp_path, record)
    with pytest.raises(ValueError, match=field):
        load_tasks(tmp_path)


def test_reused_alias_and_primary_source_need_consistent_id_and_reason(tmp_path):
    record = task()
    item = record["novelty"]["evidence"][0]
    item["study_id"] = "aLiAs"
    write(tmp_path, record)
    with pytest.raises(ValueError, match="already appears"):
        load_tasks(tmp_path)
    item["use"] = "reused"
    write(tmp_path, record)
    with pytest.raises(ValueError, match="reuse_reason"):
        load_tasks(tmp_path)
    record["novelty"]["reuse_reason"] = "Tests the disputed denominator in the prior result."
    write(tmp_path, record)
    assert ready_tasks(load_tasks(tmp_path)) == [record]
    record["novelty"]["evidence"].append(dict(item, study_id="ORIGINAL"))
    write(tmp_path, record)
    with pytest.raises(ValueError, match="duplicate study"):
        load_tasks(tmp_path)
    record["novelty"]["evidence"] = [dict(item, study_id="a new name", use="new",
                                          source="https://doi.org/10.1234/old")]
    record["sources"].append("https://doi.org/10.1234/old")
    write(tmp_path, record)
    with pytest.raises(ValueError, match="historical study Original"):
        load_tasks(tmp_path)


def test_primary_source_must_be_listed_and_usage_valid(tmp_path):
    record = task()
    item = record["novelty"]["evidence"][0]
    item["source"] += "#new"
    write(tmp_path, record)
    with pytest.raises(ValueError, match="cite"):
        load_tasks(tmp_path)
    item.update(source=record["sources"][0], use="reference?")
    write(tmp_path, record)
    with pytest.raises(ValueError, match="new or reused"):
        load_tasks(tmp_path)


@pytest.mark.parametrize("source", [
    "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?targ=self&acc=GSE24680&form=text",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24680/matrix/GSE24680_series_matrix.txt.gz",
    "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1234/?format=json",
    "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/234/E-MTAB-1234/Files/E-MTAB-1234.sdrf.txt",
    "https://DX.DOI.ORG/10.1234/%4fLD?from=paper#supplement",
])
def test_alternate_archive_and_doi_links_identify_relabelled_known_evidence(tmp_path, source):
    record = task()
    record["sources"] = [source]
    record["novelty"]["evidence"] = [{"study_id": "invented", "source": source, "use": "new"}]
    write(tmp_path, record)
    with pytest.raises(ValueError, match="historical study Original"):
        load_tasks(tmp_path)


def test_first_use_shared_across_proposals_and_completed_reuse(tmp_path):
    first, later = task("first"), task("later")
    item = later["novelty"]["evidence"][0]
    item["study_id"] = "FIRST"
    for record in (first, later):
        write(tmp_path, record)
    with pytest.raises(ValueError, match="first use"):
        load_tasks(tmp_path)
    complete(first)
    write(tmp_path, first)
    with pytest.raises(ValueError, match="first use"):
        load_tasks(tmp_path)
    item["use"] = "reused"
    later["novelty"]["reuse_reason"] = "Corrects the original donor pairing."
    complete(later)
    write(tmp_path, later)
    assert len(load_tasks(tmp_path)) == 2
    first["status"] = "ready"
    write(tmp_path, first)
    with pytest.raises(ValueError, match="already used in completed task"):
        load_tasks(tmp_path)


def test_primary_url_relabel_is_caught_but_shared_background_citation_is_allowed(tmp_path):
    first, later = task("first"), task("later")
    url = first["sources"][0]
    later["sources"].append(url)
    write(tmp_path, first)
    write(tmp_path, later)
    assert len(load_tasks(tmp_path)) == 2
    later["novelty"]["evidence"][0]["source"] = url
    write(tmp_path, later)
    with pytest.raises(ValueError, match="cohort-specific"):
        load_tasks(tmp_path)


def test_distinct_cohorts_in_an_unindexed_paper_can_cite_published_anchors(tmp_path):
    for ident in ["cohort-a", "cohort-b"]:
        record = task(ident)
        url = "https://doi.org/10.2345/large-study#" + ident
        record["sources"] = [url]
        record["novelty"]["evidence"][0]["source"] = url
        write(tmp_path, record)
    assert len(ready_tasks(load_tasks(tmp_path))) == 2


@pytest.mark.parametrize("completed_id", ["a-earlier", "z-later"])
def test_completion_doi_and_supplement_cannot_label_the_same_input_new_twice(tmp_path, history, completed_id):
    first, other = task(completed_id), task("m-proposed")
    complete(first)
    first["sources"] = ["https://doi.org/10.2345/example"]
    supplement = "https://dx.doi.org/10.2345/EXAMPLE?version=public#supplement"
    other["sources"] = [supplement]
    for record in (first, other):
        record["novelty"]["evidence"][0]["source"] = record["sources"][0]
        write(tmp_path, record)
    with pytest.raises(ValueError, match="bare DOI"):
        load_tasks(tmp_path)
    other["novelty"]["evidence"][0].update(study_id=completed_id, use="reused")
    other["novelty"]["reuse_reason"] = "Resolves the original coding assumption."
    write(tmp_path, other)
    tasks = load_tasks(tmp_path)
    _, names = research_queue.load_history(history)
    rows = research_queue.history_lines(supplement, tasks, names)
    assert len(rows) == 2 and any("done (new)" in row for row in rows)


def test_ready_prefers_new_evidence_inside_each_priority(tmp_path):
    old, fresh, urgent = task("a-reuse"), task("z-new"), task("urgent")
    old["novelty"]["evidence"][0]["use"] = "reused"
    old["novelty"]["reuse_reason"] = "Resolves an incorrect measurement in the report."
    urgent["novelty"] = old["novelty"]
    urgent["sources"] = old["sources"]
    old["priority"] = fresh["priority"] = "P1"
    for record in (old, fresh, urgent):
        write(tmp_path, record)
    assert [t["id"] for t in ready_tasks(load_tasks(tmp_path))] == ["urgent", "z-new", "a-reuse"]


def test_history_search_reports_aliases_canonical_ids_and_task_status(tmp_path, history):
    record = task()
    item = record["novelty"]["evidence"][0]
    item.update(study_id="ORIGINAL", use="reused")
    record["novelty"]["reuse_reason"] = "Corrects a prior conclusion."
    write(tmp_path, record)
    tasks = load_tasks(tmp_path, history_path=history)
    _, names = research_queue.load_history(history)
    assert research_queue.history_lines("Ali", tasks, names) == [
        "Original  historical  docs/previous.md",
        "ORIGINAL  ready (reused)  research/tasks/toy-one.json",
    ]
    for alternate in [
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24680/matrix/",
            "https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1234/?type=study",
            "https://dx.doi.org/10.1234/OLD?read=paper#supplement",
    ]:
        assert research_queue.history_lines(alternate, tasks, names) == research_queue.history_lines("Ali", tasks, names)
    assert research_queue.history_lines("unlisted", tasks, names) == []
