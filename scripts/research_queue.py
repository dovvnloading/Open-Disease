"""Offline research-task queue and local run scaffolding; no work is launched automatically."""
from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TASK_DIR = ROOT / "research" / "tasks"
ID = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
STATES = {"ready", "blocked", "done"}
PRIORITIES = {"P0", "P1", "P2"}
NETWORK = {"none", "public-opt-in"}


def load_tasks(directory: Path = TASK_DIR) -> dict[str, dict]:
    """Validate all task records, references, and dependency cycles."""
    tasks: dict[str, dict] = {}
    files = sorted(directory.glob("*.json"))
    if not files:
        raise ValueError(f"No research tasks in {directory}")
    for path in files:
        task = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(task, dict) or task.get("schema_version") != 1:
            raise ValueError(f"{path}: expected task schema_version 1")
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
    return tasks


def ready_tasks(tasks: dict[str, dict]) -> list[dict]:
    return sorted(
        (task for task in tasks.values()
         if task["status"] == "ready"
         and all(tasks[dep]["status"] == "done" for dep in task["dependencies"])),
        key=lambda task: (task["priority"], task["id"]),
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
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("validate", help="Check task integrity and dependencies")
    sub.add_parser("ready", help="Show eligible work without starting it")
    run = sub.add_parser("new-run", help="Create ignored local provenance stub")
    run.add_argument("--task", required=True)
    run.add_argument("--agent", required=True)
    run.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    tasks = load_tasks()
    if args.command == "validate":
        print(f"Validated {len(tasks)} task records, {len(ready_tasks(tasks))} ready")
    elif args.command == "ready":
        for task in ready_tasks(tasks):
            print(f'{task["priority"]}  {task["id"]}: {task["title"]}')
    else:
        print(new_run(tasks, args.task, args.agent, args.out))


if __name__ == "__main__":
    main()
