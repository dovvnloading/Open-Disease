"""Verify local GTEx v8 slices, or retrieve the pinned public intervals on request."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import Request, urlopen

from gwas_3p21_tabix import regional_rows

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "evidence/gwas-3p21-colocalization-audit/molecular_sources.json"


def checked(content, digest):
    if hashlib.sha256(content).hexdigest() != digest:
        raise ValueError("Payload differs from the recorded SHA-256")
    return content


def small_request(url, limit, accept):
    with urlopen(Request(url, headers={"Accept": accept}), timeout=45) as response:
        content = response.read(limit + 1)
    if len(content) > limit:
        raise ValueError("Response exceeded the declared byte limit")
    return content


def acquire(directory, manifest, fetch):
    directory.mkdir(parents=True, exist_ok=True)
    for tissue, source in manifest["regional_matrices"].items():
        table = directory / (tissue + ".tsv")
        if fetch and not table.exists():
            index = checked(small_request(source["url"] + ".tbi", source["tbi_bytes"],
                                          "application/octet-stream"), source["tbi_sha256"])
            header, rows, receipts = regional_rows(
                source["url"], index, "3", *manifest["grch38_envelope"], max_bytes=7_000_000
            )
            if receipts != source["ranges"] or len(rows) != source["rows"]:
                raise ValueError("Regional range or record count changed")
            temporary = table.with_suffix(".tsv.partial")
            try:
                with temporary.open("w", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(header)
                    writer.writerows(rows)
                checked(temporary.read_bytes(), source["tsv_sha256"])
                temporary.replace(table)
            finally:
                temporary.unlink(missing_ok=True)
        checked(table.read_bytes(), source["tsv_sha256"])
    for name, source in manifest["small_sources"].items():
        target = directory / name
        if fetch and not target.exists():
            content = small_request(source["url"], source["max_bytes"], source["accept"])
            if name == "GRCh38.fa.txt":
                content = content.rstrip()
            target.write_bytes(checked(content, source["sha256"]))
        checked(target.read_bytes(), source["sha256"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--fetch", action="store_true", help="Fetch missing public intervals")
    args = parser.parse_args()
    acquire(args.directory, json.loads(MANIFEST.read_text()), args.fetch)
