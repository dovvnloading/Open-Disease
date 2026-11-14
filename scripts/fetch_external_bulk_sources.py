"""Verify external-bulk inputs. Pass --fetch to acquire missing public files."""

import argparse
from hashlib import sha256
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen


ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / "evidence/external-bulk-2026-11-14"
PRIVATE = ROOT / "data/private"


def sources():
    mumbai = json.loads((EVIDENCE / "mumbai-source.json").read_text())
    us = json.loads((EVIDENCE / "us-source.json").read_text())
    for name, receipt in mumbai["sources"].items():
        yield PRIVATE / "gse153739" / name, receipt
    hgnc = mumbai["hgnc"]
    yield PRIVATE / "gse266265" / hgnc["url"].rsplit("/", 1)[1], dict(hgnc, bytes=16973125)

    sites = {
        "GPL570.annot.gz": ROOT / "data/raw/GPL570.annot.gz",
        "GSE4888_brief.soft": PRIVATE / "gse6364/talbi-brief.soft",
        "GSE51981_quick.soft": PRIVATE / "gse51981-audit/gsm-quick.soft",
    }
    for name, receipt in us["sources"].items():
        yield sites.get(name, PRIVATE / "gse6364" / name), receipt


def verify_or_fetch(path, record, fetch):
    length = record["bytes"]
    if path.is_file():
        if path.stat().st_size != length or sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"Existing file differs from source receipt: {path}")
        return 0
    if not fetch:
        raise FileNotFoundError(f"Missing {path}. Use --fetch to acquire missing inputs.")
    with urlopen(record["url"], timeout=90) as response:
        payload = response.read(length + 1)
    if len(payload) != length or sha256(payload).hexdigest() != record["sha256"]:
        raise ValueError(f"Upstream response differs from source receipt: {record['url']}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=path.parent) as staging:
        temporary = Path(staging) / path.name
        temporary.write_bytes(payload)
        temporary.replace(path)
    return length


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fetch", action="store_true", help="Download missing checksum-pinned source files")
    args = parser.parse_args()
    manifest = list(sources())
    total = sum(record["bytes"] for _, record in manifest)
    if total != 40280689 or len(manifest) != 10:
        raise ValueError("Unexpected external-bulk source count or transfer budget")
    transferred = sum(verify_or_fetch(path, record, args.fetch) for path, record in manifest)
    print(f"Verified {len(manifest)} files ({total:,} bytes), downloaded {transferred:,} bytes.")
