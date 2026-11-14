"""Validate the source catalog and its links to executed analyses offline."""
from __future__ import annotations

import json
from pathlib import Path


def audit(path: str | Path = "data/sources.json") -> dict:
    catalog = Path(path)
    data = json.loads(catalog.read_text(encoding="utf-8"))
    if not isinstance(data.get("schema_version"), str):
        raise ValueError("Invalid catalog schema")
    for field, expected in [("catalog_is_not_redistribution_permission", True),
                            ("automatic_downloads_enabled", False)]:
        if data.get(field) is not expected:
            raise ValueError(f"Incorrect catalog safety declaration: {field}")
    records = data.get("sources")
    if not isinstance(records, list) or not records:
        raise ValueError("Catalog must contain sources")
    ids = [r.get("id") for r in records]
    if any(not isinstance(i, str) or not i for i in ids) or len(ids) != len(set(ids)):
        raise ValueError("Source IDs must be unique and nonempty")
    analyzed = 0
    repo_root = catalog.resolve().parent.parent
    report_root = (repo_root / "docs/research-notes").resolve()
    for record in records:
        for field in ("url", "kind", "role", "access"):
            if not isinstance(record.get(field), str) or not record[field].strip():
                raise ValueError(f"Invalid {field} for source {record['id']}")
        if not record["url"].startswith("https://"):
            raise ValueError("Source URLs must be HTTPS")
        used = record.get("human_data_used_by_open_disease")
        if type(used) is not bool:
            raise ValueError("human_data_used_by_open_disease must be a boolean")
        if used:
            report = record.get("analysis_report")
            if not isinstance(report, str) or not report.startswith("docs/research-notes/"):
                raise ValueError("Analyzed sources must link to an analysis_report under docs/research-notes")
            target = (repo_root / report).resolve()
            if (not target.is_relative_to(report_root) or target.suffix != ".md"
                    or not target.is_file() or not target.stat().st_size):
                raise ValueError("analysis_report must point to an existing nonempty research note")
            analyzed += 1
        elif "analysis_report" in record:
            raise ValueError("A metadata-only source cannot declare an executed analysis_report")
        if not isinstance(record.get("cautions"), list) or not record["cautions"]:
            raise ValueError("Each source must list its limitations")
    if data.get("biological_results_generated") is not bool(analyzed):
        raise ValueError("biological_results_generated must agree with the analyzed sources")
    status = "exploratory_computation" if analyzed else "metadata_only"
    return {"sources": len(records), "status": status, "analyzed": analyzed}


if __name__ == "__main__":
    print(json.dumps(audit(), sort_keys=True))
