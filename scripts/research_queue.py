"""Offline research-task queue and local run scaffolding; no work is launched automatically."""
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
TASK_DIR = ROOT / "research" / "tasks"
HISTORY = ROOT / "research" / "evidence-history.json"
ID = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
STATES = {"ready", "blocked", "done"}
PRIORITIES = {"P0", "P1", "P2"}
NETWORK = {"none", "public-opt-in"}
GEO_ID = re.compile(r"(?<![a-z0-9])GSE[0-9]+(?![a-z0-9])", re.I)
ARRAYEXPRESS_ID = re.compile(r"(?<![a-z0-9])E-[A-Z]+-[0-9]+(?![a-z0-9])", re.I)


def nonblank(value: object) -> bool:
    return isinstance(value, str) and bool(value.strip())


def source_key(url: str) -> str:
    """Identify a study behind alternative official archive or DOI resolver URLs."""
    link = urlsplit(url)
    host = link.hostname or ""
    path = unquote(link.path)
    if host in {"doi.org", "dx.doi.org", "www.doi.org"}:
        # DOI names ignore case. Keep section anchors for distinct cohorts in one paper.
        doi = "https://doi.org" + path.rstrip("/").casefold()
        return doi + ("#" + link.fragment if link.fragment else "")
    if (host == "ncbi.nlm.nih.gov" or host.endswith(".ncbi.nlm.nih.gov")) and path.startswith("/geo/"):
        ids = GEO_ID.findall(path)
        for value in parse_qs(link.query).get("acc", []):
            if GEO_ID.fullmatch(value):
                ids.append(value)
    elif (host == "ebi.ac.uk" or host.endswith(".ebi.ac.uk")) and path.startswith(("/arrayexpress/", "/biostudies/")):
        ids = ARRAYEXPRESS_ID.findall(path)
    else:
        return url
    accessions = {value.casefold() for value in ids}
    if len(accessions) > 1:
        raise ValueError(f"Cite one study accession per evidence source: {url}")
    return next(iter(accessions), url)


def load_history(path: Path) -> tuple[dict[str, str], dict[str, dict]]:
    """Index the historical studies by canonical ID and known aliases."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or not isinstance(data.get("studies"), list):
        raise ValueError(f"{path}: expected an evidence history with a studies list")
    legacy = data.get("legacy_tasks")
    if (not isinstance(legacy, list)
            or any(not isinstance(i, str) or not ID.fullmatch(i) for i in legacy)
            or len(set(legacy)) != len(legacy)):
        raise ValueError(f"{path}: invalid legacy_tasks list")
    blocked = data.get("legacy_blocked_tasks")
    if (not isinstance(blocked, list) or any(not isinstance(i, str) or i not in legacy for i in blocked)
            or len(set(blocked)) != len(blocked)):
        raise ValueError(f"{path}: legacy_blocked_tasks must be distinct registered legacy tasks")
    by_name = {}
    for study in data["studies"]:
        if not isinstance(study, dict) or not nonblank(study.get("id")):
            raise ValueError(f"{path}: missing historical study ID")
        report = study.get("report")
        if (not nonblank(report) or Path(report).is_absolute() or ".." in Path(report).parts
                or not (ROOT / report).is_file()):
            raise ValueError(f"{study['id']}: history report must be an existing repository file")
        aliases = study.get("aliases")
        if not isinstance(aliases, list) or not all(nonblank(alias) for alias in aliases):
            raise ValueError(f"{study['id']}: invalid historical study aliases")
        for name in [study["id"], *aliases]:
            keys = {name.strip().casefold()}
            if name.startswith("https://"):
                keys.add(source_key(name).casefold())
            for key in keys:
                if key in by_name and by_name[key] is not study:
                    raise ValueError(f"{path}: duplicate historical study ID or alias {name}")
                by_name[key] = study
    return {ident: "blocked" if ident in blocked else "done" for ident in legacy}, by_name


def validate_novelty(tasks: dict[str, dict], by_name: dict[str, dict]) -> None:
    first_use = {}
    completed_reuse = {}
    source_ids = {}
    for ident, task in tasks.items():
        if task["schema_version"] == 1:
            continue
        novelty = task.get("novelty")
        if not isinstance(novelty, dict):
            raise ValueError(f"{ident}: missing novelty assessment")
        for field in ("prior_work", "contribution", "falsifier"):
            if not nonblank(novelty.get(field)):
                raise ValueError(f"{ident}: missing novelty.{field}")
        evidence = novelty.get("evidence")
        if not isinstance(evidence, list) or not evidence:
            raise ValueError(f"{ident}: novelty.evidence must list the study inputs")
        seen = set()
        for item in evidence:
            if not isinstance(item, dict) or not nonblank(item.get("study_id")):
                raise ValueError(f"{ident}: each evidence input needs a study_id")
            source = item.get("source")
            if source not in task["sources"]:
                raise ValueError(f"{ident}: cite the evidence source in task sources")
            study = item["study_id"].strip().casefold()
            historic = by_name.get(study)
            key = source_key(source)
            source_record = by_name.get(key.casefold()) or by_name.get(key.partition("#")[0].casefold())
            if source_record and historic != source_record:
                raise ValueError(f"{ident}: source belongs to historical study {source_record['id']}; "
                                 "use its ID and explain the reuse")
            if historic:
                study = historic["id"].strip().casefold()
            if study in seen:
                raise ValueError(f"{ident}: duplicate study ID {study} in novelty.evidence")
            seen.add(study)
            # A paper may supply more than one cohort. Use distinct cohort-specific
            # source links in that case and cite the shared paper in task.sources.
            if key in source_ids and source_ids[key] != study:
                raise ValueError(f"{ident}: source already has study ID {source_ids[key]}; "
                                 "use the same ID or cite a cohort-specific source URL")
            source_ids[key] = study
            if item.get("use") == "reused":
                if not nonblank(novelty.get("reuse_reason")):
                    raise ValueError(f"{ident}: reused evidence requires novelty.reuse_reason")
                if task["status"] == "done":
                    completed_reuse[study] = ident
            elif item.get("use") == "new":
                if historic:
                    raise ValueError(f"{ident}: {study} already appears in {historic['report']}; "
                                     "declare reused evidence and explain the reuse")
                if study in first_use:
                    raise ValueError(f"{ident}: {first_use[study]} already claims first use of {study}")
                first_use[study] = ident
            else:
                raise ValueError(f"{ident}: evidence use must be new or reused")
    # A bare DOI cannot distinguish cohorts. Check after collecting all claims
    # so earlier task filenames do not get precedence over later ones.
    for key, study in source_ids.items():
        bare, anchor, _ = key.partition("#")
        if anchor and bare.startswith("https://doi.org/") and bare in source_ids and source_ids[bare] != study:
            raise ValueError(f"{key} uses study ID {study}, while the bare DOI identifies {source_ids[bare]}; "
                             "use the same ID or cite cohort-specific evidence links for both")
    for study, ident in first_use.items():
        if study in completed_reuse and tasks[ident]["status"] != "done":
            raise ValueError(f"{ident}: {study} was already used in completed task "
                             f"{completed_reuse[study]}; declare reused evidence")


def history_lines(query: str, tasks: dict[str, dict], by_name: dict[str, dict]) -> list[str]:
    """Show reports and tasks that refer to a study ID or a registered alias."""
    if not nonblank(query):
        raise ValueError("Supply a study ID or alias to look up")
    query = query.strip()
    key = source_key(query).casefold() if query.casefold().startswith(("https://", "http://")) else query.casefold()
    terms = {query.casefold(), key}
    bare = key.partition("#")[0] if key.startswith("https://doi.org/") else None

    def matches(name: str) -> bool:
        name = name.casefold()
        return name == bare or any(term in name for term in terms)

    matched = {study["id"] for name, study in by_name.items() if matches(name)}
    lines = []
    for ident in sorted(matched):
        study = by_name[ident.strip().casefold()]
        lines.append(f'{ident}  historical  {study["report"]}')
    for ident, task in sorted(tasks.items()):
        if task["schema_version"] == 1:
            continue
        for item in task["novelty"]["evidence"]:
            name = item["study_id"].strip()
            historic = by_name.get(name.casefold())
            source = item["source"]
            if matches(name) or matches(source) or matches(source_key(source)) or (historic and historic["id"] in matched):
                lines.append(f'{name}  {task["status"]} ({item["use"]})  research/tasks/{ident}.json')
    return lines


def load_tasks(directory: Path = TASK_DIR, history_path: Path | None = None) -> dict[str, dict]:
    """Validate all task records, references, and dependency cycles."""
    legacy, by_name = load_history(HISTORY if history_path is None else history_path)
    tasks: dict[str, dict] = {}
    files = sorted(directory.glob("*.json"))
    if not files:
        raise ValueError(f"No research tasks in {directory}")
    for path in files:
        task = json.loads(path.read_text(encoding="utf-8"))
        if (not isinstance(task, dict) or type(task.get("schema_version")) is not int
                or task["schema_version"] not in (1, 2)):
            raise ValueError(f"{path}: expected task schema_version 2 or an archived version 1")
        ident = task.get("id")
        if not isinstance(ident, str) or not ID.fullmatch(ident) or path.name != ident + ".json":
            raise ValueError(f"{path}: invalid task id or filename")
        if ident in tasks:
            raise ValueError(f"Duplicate task {ident}")
        for key in ("title", "program", "question", "scope", "issue_url"):
            if not isinstance(task.get(key), str) or not task[key].strip():
                raise ValueError(f"{ident}: missing {key}")
        if not task["issue_url"].startswith("https://github.com/dovvnloading/Open-Disease/issues/"):
            raise ValueError(f"{ident}: issue_url must link to a project issue")
        if task.get("status") not in STATES or task.get("priority") not in PRIORITIES:
            raise ValueError(f"{ident}: invalid status or priority")
        if task["schema_version"] == 1 and legacy.get(ident) != task["status"]:
            raise ValueError(f"{ident}: new, reactivated or changed-status tasks require schema_version 2 and novelty")
        if task["status"] == "done":
            completed_in = task.get("completion_pr")
            if not isinstance(completed_in, str) or not re.fullmatch(
                r"https://github\.com/dovvnloading/Open-Disease/pull/[1-9][0-9]*", completed_in
            ):
                raise ValueError(f"{ident}: done task requires a project completion_pr")
        for key in ("dependencies", "steps", "deliverables", "acceptance", "sources"):
            val = task.get(key)
            if not isinstance(val, list) or (key != "dependencies" and not val):
                raise ValueError(f"{ident}: missing {key}")
            if any(not isinstance(x, str) or not x.strip() for x in val) or len(val) != len(set(val)):
                raise ValueError(f"{ident}: invalid or duplicate {key}")
        if any(not u.startswith("https://") for u in task["sources"]):
            raise ValueError(f"{ident}: source URLs must be HTTPS")
        budget = task.get("budget")
        if not isinstance(budget, dict) or budget.get("network") not in NETWORK:
            raise ValueError(f"{ident}: invalid compute/network budget")
        for key in ("max_runtime_minutes", "max_ram_gb", "max_download_mb"):
            num = budget.get(key)
            if type(num) is not int or num < (0 if key == "max_download_mb" else 1):
                raise ValueError(f"{ident}: invalid budget {key}")
        if budget["network"] == "none" and budget["max_download_mb"]:
            raise ValueError(f"{ident}: offline task must have zero download budget")
        tasks[ident] = task

    visited: set[str] = set()
    active: set[str] = set()

    def visit(ident: str) -> None:
        if ident in active:
            raise ValueError(f"Task dependency cycle at {ident}")
        if ident in visited:
            return
        active.add(ident)
        for dep in tasks[ident]["dependencies"]:
            if dep not in tasks:
                raise ValueError(f"{ident}: unknown dependency {dep}")
            visit(dep)
        active.remove(ident)
        visited.add(ident)

    for ident in tasks:
        visit(ident)
    validate_novelty(tasks, by_name)
    return tasks


def ready_tasks(tasks: dict[str, dict]) -> list[dict]:
    def order(task: dict) -> tuple:
        has_new = any(item["use"] == "new" for item in task["novelty"]["evidence"])
        return task["priority"], not has_new, task["id"]

    return sorted(
        (task for task in tasks.values()
         if task["status"] == "ready" and task["schema_version"] == 2
         and all(tasks[dep]["status"] == "done" for dep in task["dependencies"])),
        key=order,
    )


def new_run(tasks: dict[str, dict], task_id: str, agent: str, out: Path) -> Path:
    if task_id not in tasks or task_id not in {t["id"] for t in ready_tasks(tasks)}:
        raise ValueError(f"Task {task_id} is not ready; inspect its status/dependencies")
    if not agent.strip():
        raise ValueError("Agent identifier must not be blank")
    destination = out.resolve()
    if not destination.is_relative_to((ROOT / "results").resolve()) or destination.suffix != ".json":
        raise ValueError("Local run record must be a JSON file under ignored results/")
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite {destination}")
    now = datetime.now(timezone.utc).isoformat(timespec="seconds")
    record = {
        "schema_version": 1, "task_id": task_id, "status": "in_progress",
        "agent": agent, "started_utc": now, "finished_utc": None,
        "novelty": tasks[task_id]["novelty"],
        "code_commit": None, "reproduce_commands": [], "seed": None,
        "inputs": [], "outputs": [], "checks": [],
        "compute": {"runtime_seconds": None, "peak_ram_mb": None, "download_bytes": 0},
        "outcome": "not_evaluated", "summary": "", "limitations": [],
        "clinical_use": False,
    }
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    return destination


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--history-file", type=Path, default=HISTORY,
                        help="Historical study IDs, aliases, reports and legacy tasks")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("validate", help="Check task integrity and dependencies")
    sub.add_parser("ready", help="Show eligible work without starting it")
    history = sub.add_parser("history", help="Find prior reports and tasks by study ID or alias")
    history.add_argument("query")
    run = sub.add_parser("new-run", help="Create ignored local provenance stub")
    run.add_argument("--task", required=True)
    run.add_argument("--agent", required=True)
    run.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    tasks = load_tasks(history_path=args.history_file)
    if args.command == "validate":
        print(f"Validated {len(tasks)} task records, {len(ready_tasks(tasks))} ready")
    elif args.command == "ready":
        for task in ready_tasks(tasks):
            evidence = task["novelty"]["evidence"]
            new = sum(item["use"] == "new" for item in evidence)
            print(f'{task["priority"]}  {task["id"]}: {task["title"]} '
                  f'(new: {new}, reused: {len(evidence) - new})')
    elif args.command == "history":
        _, by_name = load_history(args.history_file)
        matches = history_lines(args.query, tasks, by_name)
        print("\n".join(matches) if matches else "No matching study ID or alias in the queue history.")
    else:
        print(new_run(tasks, args.task, args.agent, args.out))


if __name__ == "__main__":
    main()
