"""Verify native donor and physical-address locks before any normal-anatomy outcome."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from tempfile import NamedTemporaryFile
from pathlib import Path

from analyze_endometdb_wg6_fosb import RAW_COPY, RAW_FK
from audit_endometdb import HEAD_SHA256, pinned, read_metadata
from audit_endometdb_chip_census import (
    PREFIX_SHA256, RAW_START, RAW_TAG, ZIP_PREFIX_BYTES, classify_addresses, input_bytes,
    pg_table_blocks, raw_dump_bytes, raw_keys,
)
from audit_endometdb_probe_keys import MANIFESTS, manifest
from normal_peritoneum_design import candidate_pairs, frozen_pairs_json, verified_pairs

PRIVATE_ROOT = Path(__file__).resolve().parents[1] / "data" / "private"
CHIPS = (
    ("WG6", 48_701, 6_040_086, "ILMN_1730706"),
    ("HT12_R1", 47_231, 7_160_239, "ILMN_1751607"),
    ("HT12_R2", 47_323, 7_160_239, "ILMN_1751607"),
)


def platform_manifests(directory):
    result = {}
    for (filename, digest), (chip, size, address, probe) in zip(MANIFESTS.items(), CHIPS, strict=True):
        mapping = manifest(directory / filename, digest)
        if len(mapping) != size or mapping.get(address) != probe:
            raise ValueError("Unexpected manufacturer manifest size or identity")
        result[chip] = (set(mapping), address)
    return result


def blocks(directory):
    """Recheck the complete opaque prefix before yielding a PG raw-table block."""
    prefix = directory / "raw-prefix.zip.part"
    if prefix.stat().st_size != ZIP_PREFIX_BYTES:
        raise ValueError("Incorrect native raw-prefix size")
    with prefix.open("rb") as fh:
        if hashlib.file_digest(fh, "sha256").hexdigest() != PREFIX_SHA256:
            raise ValueError("Incorrect native raw-prefix SHA-256 before stream decoding")
    yield from pg_table_blocks(raw_dump_bytes(input_bytes(directory, prefix)), RAW_TAG)


def pair_lock(path, pairs, *, create):
    """Write once for preflight; later consumers must match the private file."""
    expected = frozen_pairs_json(pairs).encode("utf-8")
    if create:
        with NamedTemporaryFile(dir=path.parent, delete=False) as temp:
            pending = Path(temp.name)
            try:
                temp.write(expected)
                temp.flush()
                try:
                    os.link(pending, path)
                except FileExistsError:
                    pass
            finally:
                pending.unlink()
    if not path.is_file() or path.read_bytes() != expected:
        raise ValueError("Missing or changed private paired design lock")


def preflight(directory, *, create_locks=False):
    """Return only the design surviving the complete address-only census.

    CLI preflight creates missing locks. Outcome consumers default to read-only
    verification and must also pass the frozen minimum donor count.
    """
    from inflate64 import Inflater

    directory = Path(directory).resolve()
    if not directory.is_relative_to(PRIVATE_ROOT) or directory == PRIVATE_ROOT:
        raise ValueError("Use an ignored project data/private input directory")
    head = directory / "archive-start.bin"
    header = Inflater().inflate(pinned(head, HEAD_SHA256)[48:])[:RAW_START]
    if RAW_COPY not in header or RAW_FK not in header:
        raise ValueError("Native assay foreign key or signal column order changed")
    clinical = read_metadata(head, directory / "archive-tail.bin")
    candidates, selection = candidate_pairs(clinical)
    pair_lock(directory / "candidate-pairs.json", candidates, create=create_locks)
    ids = {r["id"].encode() for r in clinical["sample_analysis_result"]}
    selected = [a.encode() for pair in candidates for a in pair.peritoneum_assays + pair.uterus_assays]
    statuses, census = classify_addresses(raw_keys(blocks(directory)), selected, ids, platform_manifests(directory))
    pairs, retention = verified_pairs(candidates, {assay.decode(): chip for assay, chip in statuses.items()})
    pair_lock(directory / "verified-pairs.json", pairs, create=create_locks)
    return {"clinical_design": selection, "address_census": census, "chip_gate": retention}, pairs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    report, _ = preflight(args.directory, create_locks=True)
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
