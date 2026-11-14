"""Audit a locally obtained GEO GPL570 annotation table; never infer clinical outcomes.

Accepts tab-separated platform table (plain or gzip), including a GEO SOFT
platform_table_begin/end wrapper and the NCBI-generated annotated platform.
Emits aggregate category counts and the input checksum.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import NamedTuple, TextIO

SUBMITTER_COLUMNS = ("ID", "Gene Symbol", "ENTREZ_GENE_ID")
GEO_ANNOTATION_COLUMNS = ("ID", "Gene symbol", "Gene ID")
EMPTY = {"", "---", "NA", "N/A", "NULL", "NONE"}
SYMBOL = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
ENTREZ = re.compile(r"^[1-9][0-9]*$")


def tokens(value: str, *, id_field: bool) -> list[str] | None:
    """Return distinct tokens; None means formatting is not safely resolvable."""
    value = value.strip()
    if value.upper() in EMPTY:
        return []
    chunks = [part.strip() for part in value.split("///")]
    # The archived Affymetrix release may list numeric Entrez IDs separated
    # by whitespace. Treat all as distinct; never silently keep the first.
    if any(not chunk for chunk in chunks):
        return None
    if id_field:
        chunks = [part for chunk in chunks for part in chunk.split()]
    valid = ENTREZ if id_field else SYMBOL
    if any(not valid.fullmatch(chunk) for chunk in chunks):
        return None
    return list(dict.fromkeys(chunks))


def classify(symbol_field: str, entrez_field: str) -> str:
    return mapping_category(tokens(symbol_field, id_field=False), tokens(entrez_field, id_field=True))


def mapping_category(symbols: list[str] | None, entrez: list[str] | None) -> str:
    if symbols is None or entrez is None:
        return "invalid_annotation"
    if len(symbols) > 1 or len(entrez) > 1:
        return "multi_gene"
    if not symbols and not entrez:
        return "unmapped"
    if not symbols or not entrez:
        return "incomplete"
    return "single_symbol_single_entrez"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for part in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(part)
    return digest.hexdigest()


class TableHeader(NamedTuple):
    width: int
    probe_column: int
    symbol_column: int
    entrez_column: int
    wrapped: bool
    platform: str | None
    annotation_date: str | None


def read_table_header(stream: TextIO) -> TableHeader:
    """Select the column vocabulary from the source format and check its identity."""
    is_annotation = False
    platform = None
    date = None
    inside_table = False
    for line in stream:
        stripped = line.strip()
        if not stripped:
            continue
        key, sep, value = stripped.partition("=")
        key, value = key.strip(), value.strip()
        if stripped == "^Annotation":
            if is_annotation or platform is not None or inside_table:
                raise ValueError("Repeated or mixed SOFT declarations")
            is_annotation = True
        elif stripped.startswith("^PLATFORM") and key != "^PLATFORM":
            raise ValueError("Malformed SOFT platform accession")
        elif key in ("^PLATFORM", "!Annotation_platform"):
            if (key == "!Annotation_platform") != is_annotation or inside_table:
                raise ValueError("Platform declaration conflicts with SOFT format")
            if not sep or not re.fullmatch(r"GPL[0-9]+", value):
                raise ValueError("Malformed SOFT platform accession")
            if value != "GPL570":
                raise ValueError(f"Unexpected SOFT platform: {value}")
            if platform is not None:
                raise ValueError("Repeated SOFT platform declaration")
            platform = value
        elif key == "!Annotation_date":
            if not is_annotation or date is not None or not value or inside_table:
                raise ValueError("Invalid or repeated SOFT annotation date")
            date = datetime.strptime(value, "%b %d %Y").date().isoformat()
        elif stripped == "!platform_table_begin":
            if inside_table:
                raise ValueError("Repeated SOFT platform table begin")
            inside_table = True
        elif stripped.startswith(("^", "!", "#")):
            continue
        else:
            if (platform or is_annotation) and not inside_table:
                raise ValueError("Missing SOFT platform table begin")
            if is_annotation and (platform is None or date is None):
                raise ValueError("Missing annotated platform accession or annotation date")
            columns = GEO_ANNOTATION_COLUMNS if is_annotation else SUBMITTER_COLUMNS
            header = next(csv.reader([line], delimiter="\t"))
            if len(header) != len(set(header)) or any(col not in header for col in columns):
                raise ValueError(f"Expected unique {', '.join(columns)} columns")
            probe, symbol, entrez = (header.index(col) for col in columns)
            return TableHeader(len(header), probe, symbol, entrez, inside_table, platform, date)
    raise ValueError("Missing GPL570 annotation table")


def inspect(path: str | Path) -> dict:
    """Scan public platform annotation into a reproducible aggregate-only audit."""
    source = Path(path)
    counts: Counter[str] = Counter()
    seen: set[str] = set()
    single_mapped_ids: set[str] = set()
    probes_with_multiple_entrez = 0
    opener = gzip.open if source.suffix.lower() == ".gz" else open
    with opener(source, "rt", encoding="utf-8-sig", newline="") as stream:
        header = read_table_header(stream)
        for line in stream:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped == "!platform_table_end":
                if header.annotation_date is not None and any(rest.strip() for rest in stream):
                    raise ValueError("Unexpected records after annotated platform table")
                break
            if stripped.startswith(("!", "^", "#")):
                raise ValueError("Unexpected metadata within platform table")
            row = next(csv.reader([line], delimiter="\t"))
            if len(row) != header.width:
                raise ValueError("Malformed annotation table row")
            probe = row[header.probe_column].strip()
            if not probe or probe in seen:
                raise ValueError("Empty or duplicate platform probe ID")
            seen.add(probe)
            symbols = tokens(row[header.symbol_column], id_field=False)
            entrez = tokens(row[header.entrez_column], id_field=True)
            category = mapping_category(symbols, entrez)
            counts[category] += 1
            if entrez and len(entrez) > 1:
                probes_with_multiple_entrez += 1
            elif category == "single_symbol_single_entrez":
                single_mapped_ids.add(entrez[0])
        else:
            if header.wrapped:
                raise ValueError("Missing, empty or incomplete GPL570 annotation table")
    if not seen:
        raise ValueError("Missing, empty or incomplete GPL570 annotation table")
    return {
        "platform": "GPL570",
        "declared_platform_id": header.platform,
        "soft_platform_declaration_checked": header.platform == "GPL570" and header.wrapped,
        "table_schema": "geo_annotated_platform" if header.annotation_date else "submitter_platform",
        "declared_annotation_date": header.annotation_date,
        "input_sha256": sha256(source),
        "input_bytes": source.stat().st_size,
        "total_probes": len(seen),
        "distinct_entrez_ids_among_single_mapped_probes": len(single_mapped_ids),
        "probes_with_multiple_valid_entrez_ids": probes_with_multiple_entrez,
        "classification": {key: counts[key] for key in (
            "single_symbol_single_entrez", "multi_gene", "unmapped",
            "incomplete", "invalid_annotation")},
        "source_verified": False,
        "annotation_release": "not_auto_verified_record_release_separately",
        "biological_findings": False,
        "policy": "Only a unique symbol AND Entrez ID forms a straightforward mapping; do not infer causality. All ambiguous, unannotated and malformed mappings are retained in category totals.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True,
                        help="Locally acquired, versioned GPL570 platform table TSV or SOFT[.gz]")
    parser.add_argument("--out", type=Path, required=True,
                        help="New JSON report under ignored results/; no overwrite")
    parser.add_argument("--expected-sha256", metavar="HEX",
                        help="Optional independently pinned SHA-256 of the exact input bytes")
    args = parser.parse_args()
    if (args.expected_sha256 is not None and
            not re.fullmatch(r"[0-9A-Fa-f]{64}", args.expected_sha256)):
        parser.error("--expected-sha256 must be 64 hexadecimal characters")
    root = Path(__file__).resolve().parents[1]
    out = args.out.resolve()
    if out.suffix != ".json" or not out.is_relative_to((root / "results").resolve()):
        parser.error("Write reports only to results/*.json (ignored by Git)")
    report = inspect(args.input)
    if (args.expected_sha256 is not None and
            report["input_sha256"] != args.expected_sha256.lower()):
        parser.error("Input SHA-256 does not match --expected-sha256; no report written")
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"total_probes": report["total_probes"],
                      "classification": report["classification"]}, sort_keys=True))


if __name__ == "__main__":
    main()
