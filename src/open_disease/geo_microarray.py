"""Opt-in acquisition and strict, offline preparation of a GEO Series Matrix.

This module deliberately does NOT infer clinical labels from submitter text, analyze
biological outcomes, or contact GEO unless download_series_matrix is called.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import tempfile
from urllib.request import Request, urlopen

ACCESSION = "GSE51981"
SERIES_MATRIX_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51981/"
    "matrix/GSE51981_series_matrix.txt.gz"
)
MAX_DOWNLOAD_BYTES = 150_000_000
MANIFEST_COLUMNS = [
    "geo_accession", "donor_id", "source_title", "source_name",
    "source_characteristics_json", "source_description", "platform_id",
    "include", "reviewed", "disease_group", "phase", "tissue",
    "pathology", "hormone_treatment", "exclusion_reason", "review_source_url",
    "review_notes",
]


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _publish_without_overwrite(temporary: Path, destination: Path) -> None:
    """Publish a complete temp file atomically, only if destination is absent.

    Both paths are in the same directory/filesystem, allowing a hard link with
    exclusive-create semantics. Unlike os.replace, it cannot clobber a file
    created by a concurrent process after our initial existence check.
    """
    try:
        os.link(temporary, destination)
    except FileExistsError as exc:
        raise ValueError(f"Output already exists; refusing to overwrite: {destination}") from exc
    temporary.unlink()


def download_series_matrix(destination: str | Path, *, timeout: int = 90,
                           max_bytes: int = MAX_DOWNLOAD_BYTES) -> dict:
    """Explicitly download only the pinned public GEO matrix; fail closed on errors.

    The local sidecar records raw compressed bytes, UTC retrieval time, and URL.
    Neither a download nor a partial file overwrites an existing local input.
    """
    dest = Path(destination)
    sidecar = dest.with_name(dest.name + ".provenance.json")
    if dest.exists() or sidecar.exists():
        raise ValueError("Destination or provenance already exists; choose a new path")
    if timeout <= 0 or max_bytes <= 0:
        raise ValueError("timeout and max_bytes must be positive")
    dest.parent.mkdir(parents=True, exist_ok=True)
    temp_path: Path | None = None
    sidecar_created = False
    try:
        with tempfile.NamedTemporaryFile(dir=dest.parent, prefix=".geo-", delete=False) as temp:
            temp_path = Path(temp.name)
            request = Request(SERIES_MATRIX_URL, headers={"User-Agent": "Open-Disease/0.1 (public GEO research)"})
            with urlopen(request, timeout=timeout) as response:
                digest = hashlib.sha256()
                size = 0
                first_bytes = b""
                while True:
                    chunk = response.read(1024 * 1024)
                    if not chunk:
                        break
                    size += len(chunk)
                    if size > max_bytes:
                        raise ValueError("Download exceeded max_bytes; partial file discarded")
                    if len(first_bytes) < 2:
                        first_bytes = (first_bytes + chunk)[:2]
                    digest.update(chunk)
                    temp.write(chunk)
            if first_bytes != b"\x1f\x8b":
                raise ValueError("GEO response is not a gzipped Series Matrix")
        record = {
            "accession": ACCESSION,
            "source_url": SERIES_MATRIX_URL,
            "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
            "compressed_bytes": size,
            "sha256": digest.hexdigest(),
            "content_kind": "processed_geo_series_matrix_gzip_not_rnaseq_counts",
        }
        # Write provenance first, then publish the complete downloaded input.
        with sidecar.open("x", encoding="utf-8") as out:
            sidecar_created = True
            json.dump(record, out, indent=2, sort_keys=True)
            out.write("\n")
        _publish_without_overwrite(temp_path, dest)
        return record
    except Exception:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)
        # Never delete another process's provenance sidecar during a race.
        if sidecar_created:
            sidecar.unlink(missing_ok=True)
        raise


def _geo_row(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\r\n")], delimiter="\t", quotechar='"'))


def _metadata_column(meta: dict[str, list[list[str]]], field: str,
                     index: int) -> list[str]:
    return [values[index] for values in meta.get(field, [])]


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open("r", encoding="utf-8-sig", newline="")


def prepare_series_matrix(source: str | Path, *, manifest: str | Path,
                          expression: str | Path, audit: str | Path,
                          expected_samples: int = 148) -> dict:
    """Extract a probe x GSM table and an *unreviewed* manifest, without network.

    All outputs must be new paths; errors leave no published partial products.
    The manifest intentionally includes no automatically populated diagnoses.
    """
    src = Path(source)
    outputs = [Path(manifest), Path(expression), Path(audit)]
    if len(set(map(str, outputs))) != 3 or src in outputs:
        raise ValueError("Input and output paths must be distinct")
    if any(path.exists() for path in outputs):
        raise ValueError("An output already exists; refusing to overwrite")
    if expected_samples < 1:
        raise ValueError("expected_samples must be positive")
    for path in outputs:
        path.parent.mkdir(parents=True, exist_ok=True)

    pending: list[Path] = []
    published: list[Path] = []
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=outputs[1].parent, prefix=".geo-expression-",
                                         delete=False) as handle:
            expression_tmp = Path(handle.name)
            pending.append(expression_tmp)
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            meta: dict[str, list[list[str]]] = {}
            series: dict[str, list[str]] = {}
            samples: list[str] | None = None
            probes: set[str] = set()
            n_probes = 0
            missing_values = 0
            ended = False
            with _open_text(src) as stream:
                for line in stream:
                    token = line.strip()
                    if samples is None:
                        if token == "!series_matrix_table_begin":
                            header = _geo_row(next(stream, ""))
                            if len(header) < 2 or header[0] != "ID_REF":
                                raise ValueError("Expected ID_REF / GSM matrix header")
                            samples = header[1:]
                            if (len(samples) != expected_samples or
                                    len(set(samples)) != len(samples) or
                                    any(not re.fullmatch(r"GSM[0-9]+", s) for s in samples)):
                                raise ValueError("Unexpected count or invalid/duplicate GSM columns")
                            accessions = meta.get("!Sample_geo_accession", [])
                            if len(accessions) != 1 or accessions[0] != samples:
                                raise ValueError("GSM matrix header differs from !Sample_geo_accession")
                            for key, values in meta.items():
                                if any(len(row) != len(samples) for row in values):
                                    raise ValueError(f"Malformed sample metadata column: {key}")
                            if ACCESSION not in series.get("!Series_geo_accession", []):
                                raise ValueError("Series accession is not GSE51981")
                            writer.writerow(header)
                        elif line.startswith("!Sample_"):
                            fields = _geo_row(line)
                            meta.setdefault(fields[0], []).append(fields[1:])
                        elif line.startswith("!Series_"):
                            fields = _geo_row(line)
                            series.setdefault(fields[0], []).extend(fields[1:])
                    else:
                        if token == "!series_matrix_table_end":
                            ended = True
                            break
                        values = _geo_row(line)
                        if len(values) != len(samples) + 1:
                            raise ValueError("Expression row length does not match GSM header")
                        probe = values[0]
                        if not probe or probe in probes:
                            raise ValueError("Empty or duplicate ID_REF probe")
                        probes.add(probe)
                        normalized = [probe]
                        for value in values[1:]:
                            if value.strip().lower() in {"", "na", "null"}:
                                normalized.append("NA")
                                missing_values += 1
                            else:
                                try:
                                    valid = math.isfinite(float(value))
                                except ValueError:
                                    valid = False
                                if not valid:
                                    raise ValueError(f"Invalid numeric intensity for probe {probe}")
                                normalized.append(value)
                        writer.writerow(normalized)
                        n_probes += 1
            if samples is None or not ended or n_probes == 0:
                raise ValueError("Incomplete or empty GEO Series Matrix table")
        # Produce human-review materials only *after* successfully parsing the table.
        manifest_rows = []
        for i, sample in enumerate(samples):
            attributes = {key: _metadata_column(meta, key, i) for key in sorted(meta)}
            def first(field: str) -> str:
                values = attributes.get(field, [])
                return values[0] if values else ""
            manifest_rows.append({
                "geo_accession": sample,
                "donor_id": "",
                "source_title": first("!Sample_title"),
                "source_name": first("!Sample_source_name_ch1"),
                "source_characteristics_json": json.dumps(
                    {k: v for k, v in attributes.items() if k.startswith("!Sample_characteristics")},
                    ensure_ascii=False),
                "source_description": " | ".join(attributes.get("!Sample_description", [])),
                "platform_id": first("!Sample_platform_id"),
                "include": "false", "reviewed": "false", "disease_group": "",
                "phase": "", "tissue": "", "pathology": "",
                "hormone_treatment": "", "exclusion_reason": "pending_manual_review",
                "review_source_url": "", "review_notes": "",
            })
        digest = sha256_file(src)
        report = {
            "accession": ACCESSION,
            "source_url": SERIES_MATRIX_URL,
            "source_file_sha256": digest,
            "source_bytes": src.stat().st_size,
            "n_samples": len(samples),
            "n_probes": n_probes,
            "n_missing_intensities": missing_values,
            "included_before_review": 0,
            "status": "metadata_unreviewed_no_analysis",
            "expression_units": "GEO_processed_array_values_verify_log2_in_source",
            "manual_review_required": True,
        }
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="",
                                         dir=outputs[0].parent, prefix=".geo-manifest-",
                                         delete=False) as handle:
            manifest_tmp = Path(handle.name)
            pending.append(manifest_tmp)
            writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS, lineterminator="\n")
            writer.writeheader()
            writer.writerows(manifest_rows)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8",
                                         dir=outputs[2].parent, prefix=".geo-audit-",
                                         delete=False) as handle:
            audit_tmp = Path(handle.name)
            pending.append(audit_tmp)
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
        for tmp, dest in zip([manifest_tmp, expression_tmp, audit_tmp], outputs):
            _publish_without_overwrite(tmp, dest)
            published.append(dest)
        return report
    except Exception:
        for path in pending + published:
            path.unlink(missing_ok=True)
        raise
