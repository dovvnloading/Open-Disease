"""Audit actual cell-state support per donor and library from HECA observation fields.

Use a local H5AD, a private cache of pinned byte ranges, or explicitly request
the 1.8 MB of ranges with --fetch. Only four categorical annotation columns
are decoded.
"""
from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import io
import json
from pathlib import Path
import urllib.request

import numpy as np

from audit_heca_design_manifests import OXFORD, atlas_design, pin
from audit_heca_glandular import STATE_I, STATE_II
from open_disease.geo_microarray import sha256_file
from verify_heca_glandular import verify_state_counts

RECEIPT = Path(__file__).resolve().parents[1] / "evidence/heca-glandular-2026-11-14/annotation-ranges.json"
RECEIPT_SHA = "edd0209f4d090720b0fd62b90d748450fb5485e4d587fc1302bdf197cad037dd"


class SourceRanges(io.RawIOBase):
    """Seek through small, authenticated parts of a much larger source file."""
    def __init__(self, size, parts):
        self.size = size
        self.parts = parts
        self.starts = sorted(parts)
        self.position = 0

    def tell(self):
        return self.position

    def seek(self, offset, whence=0):
        self.position = (0, self.position, self.size)[whence] + offset
        if self.position < 0:
            raise ValueError("Negative seek")
        return self.position

    def readinto(self, buffer):
        length = max(0, min(len(buffer), self.size - self.position))
        written = 0
        while written < length:
            i = bisect_right(self.starts, self.position) - 1
            if i < 0:
                raise ValueError("Read outside verified annotation ranges")
            start = self.starts[i]
            part = self.parts[start]
            offset = self.position - start
            n = min(len(part) - offset, length - written)
            if n <= 0:
                raise ValueError("Read outside verified annotation ranges")
            buffer[written:written+n] = part[offset:offset+n]
            written += n
            self.position += n
        return length


def load_parts(manifest, directory, fetch):
    parts = {}
    directory.mkdir(parents=True, exist_ok=True, mode=0o700)
    for item in manifest["ranges"]:
        start, stop = item["start"], item["stop"]
        path = directory / f"{start}-{stop}.bin"
        if path.exists():
            chunk = path.read_bytes()
        elif fetch:
            req = urllib.request.Request(manifest["url"], headers={
                "Range": f"bytes={start}-{stop-1}", "If-Match": manifest["etag"]})
            with urllib.request.urlopen(req, timeout=45) as result:
                expected = f"bytes {start}-{stop-1}/{manifest['file_size']}"
                if (result.status != 206 or result.headers.get("Content-Range") != expected
                        or result.headers.get("ETag") != manifest["etag"]):
                    raise ValueError("Server did not return the pinned partial object")
                chunk = result.read(stop - start + 1)
        else:
            raise ValueError("Provide verified cached ranges, --h5ad, or --fetch")
        if len(chunk) != stop-start or sha256(chunk).hexdigest() != item["sha256"]:
            raise ValueError("HECA annotation range differs from its receipt")
        if not path.exists():
            path.write_bytes(chunk)
        parts[start] = chunk
    return SourceRanges(manifest["file_size"], parts)


def read_fields(source):
    import h5py

    fields = {}
    with h5py.File(source, "r") as atlas:
        for key in ("dataset", "genotype", "label_long", "library"):
            group = atlas["obs"][key]
            categories = group["categories"].asstr()[:]
            codes = group["codes"][:]
            if codes.shape != (313527,) or not np.all((codes >= 0) & (codes < len(categories))):
                raise ValueError("Incomplete aligned annotation")
            fields[key] = categories[codes]
    return fields


def coverage(counts, groups, minimum):
    donors_by_library = defaultdict(set)
    for (library, donor), count in counts.items():
        if count >= minimum:
            donors_by_library[library].add(donor)
    mixed = {library: donors for library, donors in donors_by_library.items()
             if {groups[donor] for donor in donors} == {False, True}}
    donors = set().union(*mixed.values())
    return {"minimum_per_donor_per_library": minimum,
            "mixed_libraries": len(mixed), "distinct_cases_in_mixed": sum(groups[d] for d in donors),
            "distinct_controls_in_mixed": sum(not groups[d] for d in donors),
            "gland_cells_in_mixed": sum(counts[lib, d] for lib, members in mixed.items() for d in members)}


def audit(workbook, sdrf, fields):
    clinical, sc_libraries = atlas_design(workbook)
    pin(sdrf, *OXFORD[1:])
    with sdrf.open(encoding="utf-8", newline="") as source:
        oxford = {(r["Characteristics[individual]"], r["Characteristics[library_id]"])
                  for r in csv.DictReader(source, delimiter="\t") if r["Characteristics[10x technology]"] == "scRNAseq"}
    study, donor, state, library = (fields[k] for k in ("dataset", "genotype", "label_long", "library"))
    output = {"all_state_check": verify_state_counts(workbook, study, donor, state)}
    mask = study == "Mareckova"
    if set(zip(donor[mask], library[mask], strict=True)) != oxford or set(library[mask]) != sc_libraries:
        raise ValueError("Cell annotations and Oxford library manifest differ")
    glandular = np.isin(state, [STATE_I, STATE_II])
    for key, cohort in (("Oxford", "Mareckova"), ("Huang", "Huang")):
        groups = {d: case for d, (case, _, eligible) in clinical[key].items() if eligible}
        within_study = study == cohort
        selected = within_study & glandular & np.isin(donor, list(groups))
        counts = Counter(zip(library[selected], donor[selected], strict=True))
        if {d for _, d in counts} != groups.keys():
            raise ValueError("Mismatch between available and eligible donors")
        output[key] = {
            "all_study_libraries": len(set(library[within_study])),
            "eligible_gland_cells": int(selected.sum()),
            "eligible_donors": len(groups),
            "thresholds": [coverage(counts, groups, minimum) for minimum in (1, 20)],
        }
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("sdrf", type=Path)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--parts", type=Path, help="Private cache for the 15 verified public-source ranges")
    choice.add_argument("--h5ad", type=Path, help="Local copy of the published full cell atlas")
    parser.add_argument("--fetch", action="store_true", help="Fetch only missing pinned ranges into --parts")
    args = parser.parse_args()
    if sha256_file(RECEIPT) != RECEIPT_SHA:
        raise ValueError("Unexpected range manifest")
    manifest = json.loads(RECEIPT.read_text())
    if args.h5ad:
        if args.fetch:
            parser.error("--fetch requires --parts")
        if args.h5ad.stat().st_size != manifest["file_size"]:
            raise ValueError("Wrong whole-cell object size")
        with args.h5ad.open("rb") as src:
            for part in manifest["ranges"]:
                src.seek(part["start"])
                if sha256(src.read(part["stop"] - part["start"])).hexdigest() != part["sha256"]:
                    raise ValueError("Unrecognized local HECA annotation bytes")
        fields = read_fields(args.h5ad)
    else:
        fields = read_fields(load_parts(manifest, args.parts, args.fetch))
    print(json.dumps(audit(args.workbook, args.sdrf, fields), indent=2))


if __name__ == "__main__":
    main()
