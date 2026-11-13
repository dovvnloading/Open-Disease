"""Audit a locally obtained GEO GPL570 annotation table; never infer clinical outcomes.

Accepts tab-separated platform table (plain or gzip), including a GEO SOFT
platform_table_begin/end wrapper. Outputs only aggregate category counts and
input checksum, never individual probe-to-gene assignments.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

REQUIRED = ("ID", "Gene Symbol", "ENTREZ_GENE_ID")
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
    if id_field:
        chunks = [part for chunk in chunks for part in chunk.split()]
    valid = ENTREZ if id_field else SYMBOL
    if any(not valid.fullmatch(chunk) for chunk in chunks):
        return None
    return list(dict.fromkeys(chunks))


def classify(symbol_field: str, entrez_field: str) -> str:
    symbols = tokens(symbol_field, id_field=False)
    entrez = tokens(entrez_field, id_field=True)
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


def inspect(path: str | Path) -> dict:
    """Scan public platform annotation into a reproducible aggregate-only audit."""
    source = Path(path)
    counts: Counter[str] = Counter()
    seen: set[str] = set()
    opener = gzip.open if source.suffix.lower() == ".gz" else open
    with opener(source, "rt", encoding="utf-8-sig", newline="") as stream:
        header = None
        idx = None
        inside_soft = False
        saw_end = False
        for line in stream:
            stripped = line.strip()
            if not stripped:
                continue
            if header is None:
                if stripped == "!platform_table_begin":
                    inside_soft = True
                    continue
                if stripped.startswith(("^", "!", "#")):
                    continue
                header = next(csv.reader([line], delimiter="\t"))
                if len(header) != len(set(header)) or any(x not in header for x in REQUIRED):
                    raise ValueError("Expected unique ID, Gene Symbol and ENTREZ_GENE_ID columns")
                idx = {column: header.index(column) for column in REQUIRED}
                continue
            if stripped == "!platform_table_end":
                saw_end = True
                break
            if stripped.startswith(("!", "^", "#")):
                raise ValueError("Unexpected metadata within platform table")
            row = next(csv.reader([line], delimiter="\t"))
            if len(row) != len(header):
                raise ValueError("Malformed annotation table row")
            probe = row[idx["ID"]].strip()
            if not probe or probe in seen:
                raise ValueError("Empty or duplicate platform probe ID")
            seen.add(probe)
            counts[classify(row[idx["Gene Symbol"]], row[idx["ENTREZ_GENE_ID"]])] += 1
    if header is None or not seen or (inside_soft and not saw_end):
        raise ValueError("Missing, empty or incomplete GPL570 annotation table")
    return {
        "platform": "GPL570",
        "input_sha256": sha256(source),
        "input_bytes": source.stat().st_size,
        "total_probes": len(seen),
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
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = args.out.resolve()
    if out.suffix != ".json" or not out.is_relative_to((root / "results").resolve()):
        parser.error("Write reports only to results/*.json (ignored by Git)")
    report = inspect(args.input)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps({"total_probes": report["total_probes"],
                      "classification": report["classification"]}, sort_keys=True))


if __name__ == "__main__":
    main()
