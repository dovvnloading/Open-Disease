"""Audit GSE51981 GEO Series Matrix headers, not human expression values.

No network, automatic clinical label assignment, or individual-record output.
A green structural report is NOT a clinical metadata, rights, or donor audit.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
from urllib.parse import parse_qs, urlsplit

ACCESSION = "GSE51981"
PLATFORM = "GPL570"
EXPECTED_SAMPLES = 148
GSM = re.compile(r"GSM[1-9][0-9]*\Z")
PHASES = {"proliferative", "early_secretory", "mid_secretory"}
GROUPS = {"case", "control", "other"}
ROOT = Path(__file__).resolve().parents[1]
MAX_METADATA_CHARS = 20_000_000
REQUIRED_MANIFEST = {
    "geo_accession", "donor_id", "include", "reviewed", "disease_group",
    "phase", "tissue", "pathology", "hormone_treatment",
    "exclusion_reason", "review_source_url",
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _row(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\r\n")], delimiter="\t", quotechar='"'))


def _read_metadata(source: Path) -> tuple[list[str], dict[str, list[list[str]]], dict]:
    """Stop at matrix-table BEGIN, never deserialize expression intensities."""
    meta: dict[str, list[list[str]]] = {}
    series: dict[str, list[list[str]]] = {}
    consumed = 0
    found_begin = False
    opener = gzip.open if source.suffix.lower() == ".gz" else open
    with opener(source, "rt", encoding="utf-8-sig", newline="") as stream:
        for line in stream:
            consumed += len(line)
            if consumed > MAX_METADATA_CHARS:
                raise ValueError("Metadata header exceeds bounded inspection size")
            if line.strip() == "!series_matrix_table_begin":
                found_begin = True
                break
            if line.startswith(("!Sample_", "!Series_")):
                fields = _row(line)
                if len(fields) < 2:
                    raise ValueError("Malformed GEO metadata line")
                target = meta if fields[0].startswith("!Sample_") else series
                target.setdefault(fields[0], []).append(fields[1:])
    if not found_begin:
        raise ValueError("Missing matrix table begin; source is not a Series Matrix")
    if series.get("!Series_geo_accession") != [[ACCESSION]]:
        raise ValueError("Series accession must be precisely GSE51981")
    accessions = meta.get("!Sample_geo_accession", [])
    if len(accessions) != 1:
        raise ValueError("Expected exactly one GSM accession row")
    samples = accessions[0]
    if len(samples) != EXPECTED_SAMPLES or len(set(samples)) != EXPECTED_SAMPLES or any(
        not GSM.fullmatch(s) for s in samples
    ):
        raise ValueError("Expected 148 unique valid GEO sample accessions")
    for field, rows in meta.items():
        if any(len(row) != EXPECTED_SAMPLES for row in rows):
            raise ValueError("Malformed or misaligned per-sample metadata row")
    platforms = meta.get("!Sample_platform_id", [])
    if len(platforms) != 1 or set(platforms[0]) != {PLATFORM}:
        raise ValueError("Expected one GPL570 platform ID for every sample")
    if "!Series_platform_id" in series and set(sum(series["!Series_platform_id"], [])) != {PLATFORM}:
        raise ValueError("Series platform is not GPL570")
    coverage = Counter()
    unparseable = 0
    repeated_source_keys = 0
    for index in range(EXPECTED_SAMPLES):
        keys: set[str] = set()
        for fields in meta.get("!Sample_characteristics_ch1", []):
            field = fields[index]
            if ":" not in field:
                unparseable += 1
                continue
            k, _ = field.split(":", 1)
            key = " ".join(k.lower().strip().replace("_", " ").split())
            if key in keys:
                repeated_source_keys += 1
            keys.add(key)
        # Only presence of *source field names*, never infer its clinical value.
        tests = {
            "tissue": any("tissue" in k for k in keys),
            "endometriosis": any("endometriosis" in k for k in keys),
            "other_pathology": any("pathology" in k for k in keys),
            "explicit_cycle_phase": any("phase" in k or "cycle" in k for k in keys),
            "explicit_person_id": any(word in k for k in keys for word in
                                      ("donor", "participant", "subject", "patient")),
            "hormone_or_treatment": any("hormon" in k or "treatment" in k for k in keys),
            "data_processing": bool(meta.get("!Sample_data_processing", []) and
                                    meta["!Sample_data_processing"][0][index].strip()),
        }
        coverage.update(k for k, present in tests.items() if present)
    return samples, meta, {
        "source_field_presence_n": {key: coverage[key] for key in (
            "tissue", "endometriosis", "other_pathology", "explicit_cycle_phase",
            "explicit_person_id", "hormone_or_treatment", "data_processing")},
        "unparseable_characteristics_entries": unparseable,
        "repeated_characteristic_keys": repeated_source_keys,
    }


def _flag(s: str) -> bool:
    if s.strip().lower() not in {"true", "false"}:
        raise ValueError("Manifest include/reviewed flags must be explicit true/false")
    return s.strip().lower() == "true"


def _is_sample_source(url: str, gsm: str) -> bool:
    try:
        p = urlsplit(url)
    except ValueError:
        return False
    return (p.scheme == "https" and p.hostname == "www.ncbi.nlm.nih.gov"
            and p.path == "/geo/query/acc.cgi" and gsm in parse_qs(p.query).get("acc", []))


def check_manifest(path: Path, samples: list[str]) -> dict:
    """Structural gate only: cannot authenticate source labels or person identities."""
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames is None or not REQUIRED_MANIFEST.issubset(reader.fieldnames):
            raise ValueError("Reviewed manifest lacks required audit fields")
        rows = list(reader)
    ids = [r["geo_accession"] for r in rows]
    if len(rows) != EXPECTED_SAMPLES or len(set(ids)) != EXPECTED_SAMPLES or set(ids) != set(samples):
        raise ValueError("Manifest must have exactly the same 148 GSM identifiers")
    groups: Counter[str] = Counter()
    reviewed = 0
    excluded = 0
    pending = 0
    included_ids: set[str] = set()
    for row in rows:
        inc, rev = _flag(row["include"]), _flag(row["reviewed"])
        if rev:
            reviewed += 1
        if not inc:
            excluded += 1
            if not rev:
                pending += 1
            if not isinstance(row["exclusion_reason"], str) or not row["exclusion_reason"].strip():
                raise ValueError("All excluded manifest rows need reasons")
            continue
        if not rev:
            raise ValueError("A sample included in a model has not been reviewed")
        for col in ("donor_id", "disease_group", "phase", "tissue", "pathology",
                    "hormone_treatment", "review_source_url"):
            if not isinstance(row[col], str) or not row[col].strip():
                raise ValueError("Included entry has missing mandatory review field")
        if row["disease_group"] not in GROUPS or row["phase"] not in PHASES:
            raise ValueError("Included entry has unapproved group or phase")
        if row["tissue"] != "eutopic_endometrium" or row["hormone_treatment"] != "none":
            raise ValueError("Included entry has incompatible tissue or treatment")
        if not _is_sample_source(row["review_source_url"], row["geo_accession"]):
            raise ValueError("Included entry must cite its own primary GSM page")
        if row["donor_id"] in included_ids:
            raise ValueError("Duplicate donor key among included samples")
        included_ids.add(row["donor_id"])
        groups[row["disease_group"]] += 1
    return {
        "reviewed_n": reviewed, "excluded_n": excluded, "pending_review_n": pending,
        "included_n": sum(groups.values()),
        "included_groups_n": {k: groups[k] for k in ("case", "control", "other")},
        "review_manifest_sha256": sha256_file(path),
        "status": "structural_preflight_only_not_independent_clinical_review",
        "source_label_authenticity_verified_by_script": False,
    }


def audit(matrix: Path, *, manifest: Path | None = None) -> dict:
    samples, _, coverage = _read_metadata(matrix)
    report = {
        "accession": ACCESSION,
        "platform": PLATFORM,
        "n_geo_accessions": len(samples),
        "source_file_sha256": sha256_file(matrix),
        "source_bytes": matrix.stat().st_size,
        "source_kind": "GEO_series_matrix_metadata_only",
        "expression_intensities_parsed_or_analyzed": False,
        "source_field_coverage": coverage,
        "paper_phase_assigned_counts_not_local_observations": {
            "case": 75, "control": 34, "other": 35,
            "indeterminate_case": 2, "indeterminate_other": 2,
            "source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/",
        },
        "manifest": check_manifest(manifest, samples) if manifest else None,
        "source_rights_verified": False,
        "donor_independence_verified_by_script": False,
        "preanalysis_approval": "NO_GO_manual_source_and_independence_audit_outstanding",
        "biological_findings": False,
    }
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True, help="Locally obtained official GEO matrix (.gz)")
    parser.add_argument("--manifest", type=Path, help="Optional local 148-row reviewed manifest")
    parser.add_argument("--out", type=Path, required=True, help="Aggregate-only report under ignored results/")
    args = parser.parse_args()
    dest = args.out.resolve()
    if dest.suffix != ".json" or not dest.is_relative_to((ROOT / "results").resolve()):
        parser.error("Reports must be JSON inside ignored results/")
    if dest.exists():
        parser.error("Will not overwrite existing report")
    if dest == args.matrix.resolve() or (args.manifest is not None and dest == args.manifest.resolve()):
        parser.error("Output must not overwrite any input")
    report = audit(args.matrix, manifest=args.manifest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"accession": ACCESSION, "n_samples": report["n_geo_accessions"],
                      "status": report["preanalysis_approval"]}))


if __name__ == "__main__":
    main()
