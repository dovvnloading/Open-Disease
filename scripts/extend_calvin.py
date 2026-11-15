"""Opt in to retrieving 6,144 more compressed bytes from each of 50 public Calvin CEL files."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from pathlib import Path
import urllib.request
import zlib

from audit_us_bulk_sources import linked_files

PINS = {
    "GSE51981": ("767499b710e33598cffcc7cdaf414925630ec241d90fcc04bcfec93b903bd5d8", 148),
    "GSE29981": ("228c904d5a4d984709ea1d658aba24a8e6bef1852606187917c87d3e76a4a75c", 20),
}
PREFIXES = "5b73bf809a819992c14b6d7b7421875d4ef664cff72b9aa07583843b6b336e3e"


def collect(directory, out):
    jobs = []
    for study, (digest, count) in PINS.items():
        soft = directory / (study + ".soft")
        if hashlib.sha256(soft.read_bytes()).hexdigest() != digest:
            raise ValueError("Unrecognized GEO metadata snapshot")
        (out / study).mkdir(parents=True, exist_ok=True)
        for key, url in linked_files(soft, study, count):
            prefix = (directory / study / (key + ".part")).read_bytes()
            if len(prefix) != 2048:
                raise ValueError("Unexpected compressed header range")
            if zlib.decompressobj(31).decompress(prefix, 16384).startswith(b"\x3b\x01"):
                jobs.append((prefix, url, out / study / (key + ".part")))
    if len(jobs) != 50:
        raise ValueError("Expected 50 Calvin supplementary files")

    def fetch(job):
        first, url, path = job
        if path.is_file():
            if path.stat().st_size != 8192 or not path.read_bytes().startswith(first):
                raise ValueError("Invalid cached CEL prefix")
            return 0
        request = urllib.request.Request(url, headers={"Range": "bytes=2048-8191"})
        with urllib.request.urlopen(request, timeout=45) as response:
            if response.status != 206 or not response.headers.get("Content-Range", "").startswith("bytes 2048-8191/"):
                raise ValueError("Range response not honored")
            extra = response.read(6145)
        if len(extra) != 6144:
            raise ValueError("Unexpected range length")
        tmp = path.with_suffix(".tmp")
        tmp.write_bytes(first + extra)
        tmp.replace(path)
        return len(extra)

    with ThreadPoolExecutor(max_workers=4) as pool:
        total = sum(pool.map(fetch, jobs))
    summary = hashlib.sha256()
    for _, _, path in sorted(jobs, key=lambda job: job[2]):
        summary.update(path.read_bytes())
    if summary.hexdigest() != PREFIXES:
        raise ValueError("Calvin prefix aggregate checksum mismatch")
    return total


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("out", type=Path)
    print("Additional bytes retrieved", collect(**vars(p.parse_args())))
