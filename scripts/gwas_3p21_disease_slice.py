"""Extract and verify the GCST90205183 sensitivity interval without a whole-file download."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request

from gwas_3p21_tabix import regional_rows


BASE = (
    "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
    "GCST90205001-GCST90206000/GCST90205183/harmonised/GCST90205183.h.tsv.gz"
)
INTERVAL = ("3", 48_614_715, 50_614_717)
INDEX_SHA256 = "bc8904cc2a25332ea9ed952f5031be3c9e7c4510cf9174ab8989e6d6ede772ff"
SLICE_SHA256 = "8f996e0bf80ae9d79b10450e4449b96cc92a4f907c43ac471ae65ad47f6f2739"


def verify(payload, digest):
    if hashlib.sha256(payload).hexdigest() != digest:
        raise ValueError("Source digest differs from the verified release")
    return payload


def extract(directory, fetch):
    directory.mkdir(parents=True, exist_ok=True)
    tranche = directory / "gwas_region.tsv"
    if not fetch:
        verify(tranche.read_bytes(), SLICE_SHA256)
        return tranche
    with urllib.request.urlopen(BASE + ".tbi", timeout=45) as response:
        index = verify(response.read(2_000_001), INDEX_SHA256)
    header, rows, ranges = regional_rows(BASE, index, *INTERVAL, max_bytes=10_000_000)
    out = io.StringIO()
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    payload = verify(out.getvalue().encode(), SLICE_SHA256)
    (directory / "gwas_range_receipt.json").write_text(json.dumps(ranges, indent=2) + "\n")
    tranche.write_bytes(payload)
    return tranche


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Local untracked output directory")
    parser.add_argument("--fetch", action="store_true", help="Read the public TBI and bounded BGZF ranges")
    args = parser.parse_args()
    print(extract(args.directory, args.fetch))
