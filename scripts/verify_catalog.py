"""Validate the canonical metadata-only data/sources.json catalog offline."""
from __future__ import annotations

import json
from pathlib import Path


def audit(path: str | Path = "data/sources.json") -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data.get("schema_version"), str):
        raise ValueError("Invalid catalog schema")
    for field, expected in [("catalog_is_not_redistribution_permission", True),
                            ("automatic_downloads_enabled", False),
                            ("biological_results_generated", False)]:
        if data.get(field) is not expected:
            raise ValueError(f"Incorrect catalog safety declaration: {field}")
    records = data.get("sources")
    if not isinstance(records, list) or not records:
        raise ValueError("Catalog must contain sources")
    ids = [r.get("id") for r in records]
    if any(not isinstance(i, str) or not i for i in ids) or len(ids) != len(set(ids)):
        raise ValueError("Source IDs must be unique and nonempty")
    for record in records:
        for field in ("url", "kind", "role", "access"):
            if not isinstance(record.get(field), str) or not record[field].strip():
                raise ValueError(f"Invalid {field} for source {record['id']}")
        if not record["url"].startswith("https://"):
            raise ValueError("Source URLs must be HTTPS")
        if record.get("human_data_used_by_open_disease") is not False:
            raise ValueError("Real human-data analysis must not be claimed by this metadata catalog")
        if not isinstance(record.get("cautions"), list) or not record["cautions"]:
            raise ValueError("Each source must list its limitations")
    return {"sources": len(records), "status": "metadata_only", "analyzed": 0}


if __name__ == "__main__":
    print(json.dumps(audit(), sort_keys=True))
