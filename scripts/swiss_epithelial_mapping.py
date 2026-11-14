"""Project the locked nominations onto the Swiss epithelial gene labels without outcomes.

The reference is the 2026-10-06 HGNC complete set (permanent quarterly archive).
Only approved symbols and unversioned Ensembl gene accessions participate.
No expression values or sample characteristics enter this mapping.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import re

from open_disease.geo_microarray import sha256_file

ROOT = Path(__file__).resolve().parents[1]
LOCK_FILE = ROOT / "evidence/cross-cohort-2026-11-14/nominations.json"
LOCK_SHA256 = "ba140f2bb026b0416e1c705ab28f616f430cc4295f94d174feaf5d15947c1e31"
HGNC_FILE = "hgnc_complete_set_2026-10-06.tsv"
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/"
    + HGNC_FILE
)
HGNC_SIZE = 16_973_125
HGNC_SHA256 = "8bf6f686e640dbc864a346484626953621b0500d900d26cd085a35c2782d725d"
LABEL_SHA256 = "4f0d87790afbc35d8af755a3148ed09883ca57d9a309f5e54dfdc87feb03ad14"
LABEL_COUNT = 29_034
ENTREZ = re.compile(r"[1-9][0-9]*")
ENSEMBL = re.compile(r"(ENSG[0-9]{11})(?:\.[0-9]+)?")
COLUMNS = ("family", "id", "entrez_id", "feature_label", "row_index")


def crosswalk(hgnc, labels):
    """Keep one approved HGNC entry and one observed Swiss row per Entrez ID."""
    symbols = defaultdict(set)
    accessions = defaultdict(set)
    entries = defaultdict(set)
    for record in hgnc:
        entrez = record["entrez_id"]
        if record["status"] != "Approved" or not ENTREZ.fullmatch(entrez):
            continue
        entries[entrez].add(record["hgnc_id"])
        symbols[record["symbol"]].add(entrez)
        ensembl = record["ensembl_gene_id"]
        if re.fullmatch(r"ENSG[0-9]{11}", ensembl):
            accessions[ensembl].add(entrez)

    if any(not label for label in labels) or len(set(labels)) != len(labels):
        raise ValueError("Feature labels must be nonempty and unique")
    matches = defaultdict(list)
    counts = {
        "source_rows": len(labels), "no_match": 0, "ambiguous_label": 0,
        "multiple_hgnc_entries": 0, "multiple_feature_rows": 0,
        "matched_symbol": 0, "matched_ensembl": 0,
    }
    for index, label in enumerate(labels, start=1):
        accession = ENSEMBL.fullmatch(label)
        if accession:
            candidates = accessions.get(accession.group(1), set())
            kind = "matched_ensembl"
        else:
            candidates = symbols.get(label, set())
            kind = "matched_symbol"
        if not candidates:
            counts["no_match"] += 1
        elif len(candidates) > 1:
            counts["ambiguous_label"] += 1
        else:
            entrez, = candidates
            if len(entries[entrez]) != 1:
                counts["multiple_hgnc_entries"] += 1
            else:
                matches[entrez].append((label, index, kind))

    mapping = {}
    for entrez, rows in matches.items():
        if len(rows) != 1:
            counts["multiple_feature_rows"] += len(rows)
            continue
        label, index, kind = rows[0]
        counts[kind] += 1
        mapping[entrez] = (label, index)
    return mapping, counts


def candidate_projection(candidates, mapping):
    rows = []
    seen = set()
    for item in candidates:
        family, candidate = item["family"], item["id"]
        if candidate in seen:
            raise ValueError(f"Duplicate candidate {candidate}")
        seen.add(candidate)
        if family == "gene":
            members = [candidate]
        elif family == "pathway":
            members = list(map(str, item["entrez_ids"]))
        else:
            raise ValueError(f"Unexpected candidate family {family}")
        if len(set(members)) != item["members"] or len(members) != item["members"]:
            raise ValueError(f"Inconsistent locked membership for {candidate}")
        for entrez in members:
            if not ENTREZ.fullmatch(entrez):
                raise ValueError(f"Invalid locked Entrez ID {entrez}")
            label, index = mapping.get(entrez, ("", ""))
            rows.append((family, candidate, entrez, label, index))
    return rows


def load_projection(directory, lock_path=LOCK_FILE, labels_path=None):
    directory = Path(directory)
    hgnc_path = directory / HGNC_FILE
    labels_path = directory / "row_labels.txt" if labels_path is None else Path(labels_path)
    if hgnc_path.stat().st_size != HGNC_SIZE or sha256_file(hgnc_path) != HGNC_SHA256:
        raise ValueError("HGNC quarterly archive size or SHA-256 differs")
    if sha256_file(labels_path) != LABEL_SHA256:
        raise ValueError("Swiss ordered feature label SHA-256 differs")
    if sha256_file(lock_path) != LOCK_SHA256:
        raise ValueError("Original nominations SHA-256 differs")
    labels = labels_path.read_text(encoding="utf-8").splitlines()
    if len(labels) != LABEL_COUNT:
        raise ValueError("Unexpected Swiss feature count")
    with hgnc_path.open(encoding="utf-8", newline="") as handle:
        mapping, counts = crosswalk(csv.DictReader(handle, delimiter="\t"), labels)
    candidates = json.loads(Path(lock_path).read_text(encoding="utf-8"))["candidates"]
    return candidate_projection(candidates, mapping), counts


def write_projection(path, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(rows)


def source_receipt(counts, rows, table):
    """Describe the exact input artifacts and outcome-blind coverage."""
    coverage = {}
    measured = set()
    pathway_genes = set()
    for family, key, gene, _, index in rows:
        item = coverage.setdefault(key, {"family": family, "mapped": 0, "original": 0})
        item["original"] += 1
        if index:
            item["mapped"] += 1
            measured.add(gene)
            if family == "pathway":
                pathway_genes.add(gene)
    return {
        "source": "HUGO Gene Nomenclature Committee, complete approved set",
        "reference_date": "2026-10-06", "reference_url": HGNC_URL,
        "reference_license": "CC0", "reference_bytes": HGNC_SIZE,
        "reference_sha256": HGNC_SHA256,
        "ordered_feature_labels_sha256": LABEL_SHA256,
        "original_nominations_sha256": LOCK_SHA256,
        "projection_sha256": sha256_file(table),
        "source_row_counts": counts, "coverage": coverage,
        "distinct_mapped_genes": len(measured),
        "distinct_mapped_pathway_genes": len(pathway_genes),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("directory", type=Path, help="Directory containing HGNC archive and row_labels.txt")
    ap.add_argument("--out", type=Path, required=True, help="Gene-label-only projected memberships TSV")
    ap.add_argument("--receipt", type=Path, help="Public aggregate provenance and coverage JSON")
    args = ap.parse_args()
    rows, counts = load_projection(args.directory)
    write_projection(args.out, rows)
    report = json.dumps(source_receipt(counts, rows, args.out), indent=2)
    if args.receipt:
        args.receipt.write_text(report + "\n", encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
