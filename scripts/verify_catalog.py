"""Audit metadata-only manifest. Deliberately performs no network or data access."""
from __future__ import annotations
import json
from pathlib import Path


def audit(path: str | Path = "catalog/sources.json") -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("schema_version") != 1 or data.get("clinical_data_ingested") is not False:
        raise ValueError("Unsupported manifest or human-data ingestion assertion")
    sources = data.get("datasets")
    if not isinstance(sources, list) or not sources:
        raise ValueError("Registry must list sources")
    ids = [entry.get("id") for entry in sources]
    if any(not isinstance(i, str) or not i.strip() for i in ids) or len(ids) != len(set(ids)):
        raise ValueError("Dataset IDs must be nonempty and unique")
    for entry in sources:
        for field in ("source_url", "publication_url", "role", "license_status", "human_data_access_status", "source_note"):
            if not isinstance(entry.get(field), str) or not entry[field].strip():
                raise ValueError(f"Missing {field} for {entry['id']}")
        if not entry["source_url"].startswith("https://"):
            raise ValueError("Only public HTTPS source references are allowed")
        if entry.get("analyzed_in_this_repo") is not False:
            raise ValueError("Do not claim real patient data were analyzed")
        sha = entry.get("sha256")
        if sha is not None and (not isinstance(sha, str) or len(sha) != 64
                                or any(c not in "0123456789abcdef" for c in sha)):
            raise ValueError("Invalid SHA256")
    return {"sources": len(sources), "status": "metadata_only", "analyzed": 0}


if __name__ == "__main__":
    print(json.dumps(audit(), sort_keys=True))
