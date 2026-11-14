"""Census native assay probe addresses. Read an offline Figshare v1 ZIP prefix.

Supply a single pinned prefix through byte 512595210 with --prefix, or reuse
the documented fragments, including raw-parts/, under the input directory.
Only native assay IDs and integer probe addresses are decoded from raw rows.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import zlib

from audit_endometdb import HEAD_SHA256, eligible, pinned, read_metadata
from audit_endometdb_probe_keys import MANIFESTS, manifest

ZIP_PREFIX_BYTES = 512_595_211
RAW_START = 38_714
RAW_END = 512_632_179
RAW_TAG = bytes.fromhex("0100d8080000")
PREFIX_SHA256 = "74316094156c43ebb693c060178238b46c0b3f5fa032f8333488063da9e5ce7a"
CHUNK = 256 * 1024


def input_bytes(directory, prefix):
    # Fresh checkouts typically download the single prefix documented by the census.
    if prefix is None and (directory / "raw-prefix.zip.part").is_file():
        prefix = directory / "raw-prefix.zip.part"
    if prefix is not None:
        if prefix.stat().st_size != ZIP_PREFIX_BYTES:
            raise ValueError("Wrong Figshare prefix length")
        spans = [(prefix, ZIP_PREFIX_BYTES)]
    else:
        spans = [(directory / "archive-start.bin", 65_536), (directory / "native-raw-continuation.bin", 2_600_000)]
        for start in range(2_665_536, 510_400_000, 8_000_000):
            end = min(start + 7_999_999, 510_399_999)
            spans.append((directory / "raw-parts" / f"{start:09d}-{end:09d}.bin", end - start + 1))
        spans.append((directory / "native-normalized-head-range.bin", ZIP_PREFIX_BYTES - 510_400_000))
    digest = hashlib.sha256()
    for path, n in spans:
        with path.open("rb") as source:
            while n:
                chunk = source.read(min(CHUNK, n))
                if not chunk:
                    raise ValueError("Short input range")
                n -= len(chunk)
                digest.update(chunk)
                yield chunk
    if digest.hexdigest() != PREFIX_SHA256:
        raise ValueError("Unexpected Figshare raw ZIP prefix SHA-256")


def pg_table_blocks(chunks, header):
    """Parse pg_dump compression frames, including partial headers and payloads."""
    pending = bytearray()
    started = False
    ended = False
    for chunk in chunks:
        pending.extend(chunk)
        if not started:
            if len(pending) < len(header):
                continue
            if not pending.startswith(header):
                raise ValueError("Unexpected pg_dump table tag")
            del pending[:len(header)]
            started = True
        while len(pending) >= 5 and not ended:
            size = int.from_bytes(pending[1:5], "little")
            if pending[0] != 0 or size > 4096:
                raise ValueError("Unexpected pg_dump compression frame")
            if len(pending) < 5 + size:
                break
            if size == 0:
                ended = True
            else:
                yield bytes(pending[5:5 + size])
            del pending[:5 + size]
        if ended and pending:
            raise ValueError("Bytes beyond the raw PG terminator")
    if not ended:
        raise ValueError("Missing raw PG terminator")


def raw_dump_bytes(chunks):
    from inflate64 import Inflater

    decoder = Inflater()
    position = zip_position = 0
    for chunk in chunks:
        # The sole leading 48 bytes are the pinned ZIP member's local header.
        decoded = decoder.inflate(chunk[max(0, 48 - zip_position):])
        zip_position += len(chunk)
        lo, hi = max(0, RAW_START - position), max(0, RAW_END - position)
        if lo < len(decoded) and hi > lo:
            yield decoded[lo:hi]
        position += len(decoded)
    if position != RAW_END or zip_position != ZIP_PREFIX_BYTES:
        raise ValueError("Wrong ZIP or uncompressed PG boundary")


def copy_lines(blocks):
    """Yield complete opaque COPY record bytes and verify both stream endings."""
    decoder = zlib.decompressobj()
    carry = b""
    terminated = False
    for block in blocks:
        decoded = decoder.decompress(block)
        if decoder.unused_data:
            raise ValueError("Unexpected bytes beyond raw zlib stream")
        *lines, carry = (carry + decoded).split(b"\n")
        for line in lines:
            if line == b"\\." and not terminated:
                terminated = True
            elif terminated:
                # pg_dump may append an empty line after its text COPY terminator.
                if line:
                    raise ValueError("Raw record after COPY terminator")
            else:
                yield line
    if not decoder.eof or carry or not terminated:
        raise ValueError("Incomplete raw COPY stream")


def raw_keys(blocks):
    for line in copy_lines(blocks):
        # Leave all quantitative columns in the opaque fourth field.
        _, assay, address, _ = line.split(b"\t", 3)
        yield assay, address


def native_chip(observed, n_duplicates, manifests):
    if not observed:
        return "missing"
    if n_duplicates:
        return "repeated_addresses"
    matches = [name for name, (expected, _) in manifests.items() if observed == expected]
    if not matches:
        return "no_exact_manifest"
    if len(matches) > 1:
        return "ambiguous_manifest"
    chip = matches[0]
    if manifests[chip][1] not in observed:
        raise ValueError("Verified chip lacks FOSB manufacturer address")
    return chip


def census(rows, selected, native_ids, manifests, *, selected_status=None):
    """Report aggregates; optionally retain chip classifications in caller memory."""
    records = {r["analysis_result_id"].encode(): r for r in selected}
    if len(records) != len(selected):
        raise ValueError("Duplicate clinically selected assay")
    unknown = set()
    seen_assays = set()
    addresses = defaultdict(set)
    duplicates = Counter()
    total = 0
    for assay, address in rows:
        total += 1
        seen_assays.add(assay)
        if assay not in native_ids:
            unknown.add(assay)
        if assay in records:
            numeric = int(address)
            if numeric in addresses[assay]:
                duplicates[assay] += 1
            addresses[assay].add(numeric)
    if unknown:
        raise ValueError("Unmapped raw assay outside native clinical namespace")

    by_status = Counter()
    strata = defaultdict(Counter)
    models = {}
    for assay, record in records.items():
        observed = addresses[assay]
        status = native_chip(observed, duplicates[assay], manifests)
        if selected_status is not None:
            selected_status[assay] = status
        by_status[status] += 1
        strata[status, record["analysis_run_id"], record["phase"]][record["group"]] += 1
    for chip in manifests:
        cells = [c for (status, _, _), c in strata.items() if status == chip and c["case"] and c["control"]]
        models[chip] = {"shared_cells": len(cells), **{g: sum(c[g] for c in cells) for g in ("case", "control")}}
    return {
        "native_raw_rows": total, "native_raw_assays": len(seen_assays),
        "native_assay_rows": len(native_ids), "clinically_selected": len(records),
        "selected_assays_with_repeated_addresses": sum(n > 0 for n in duplicates.values()),
        "selected_excess_address_rows": sum(duplicates.values()),
        "selected_by_status": dict(sorted(by_status.items())),
        "overlap_by_complete_manifest": models,
        "selected_by_status_run_phase": [
            {"status": status, "run": run, "phase": phase, "case": c["case"], "control": c["control"]}
            for (status, run, phase), c in sorted(strata.items(), key=lambda item: (item[0][0], int(item[0][1]), item[0][2]))
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Directory containing the pinned clinical segments and manifests")
    parser.add_argument("--prefix", type=Path, help="One ZIP prefix file covering exactly the documented byte range")
    args = parser.parse_args()
    directory = args.directory
    head, tail = directory / "archive-start.bin", directory / "archive-tail.bin"
    from inflate64 import Inflater

    header = Inflater().inflate(pinned(head, HEAD_SHA256)[48:])[:RAW_START]
    fk = b"ADD CONSTRAINT sample_analysis_result_microarray_data_fk FOREIGN KEY (sample_analysis_result_id) REFERENCES public.sample_analysis_result(id)"
    if fk not in header:
        raise ValueError("Missing native raw-to-assay foreign key")
    clinical = read_metadata(head, tail)
    _, selected = eligible(clinical)
    native_ids = {row["id"].encode() for row in clinical["sample_analysis_result"]}
    manifests = {}
    for (filename, digest), (label, fosb_address, fosb_probe) in zip(
        MANIFESTS.items(), [("WG6", 6040086, "ILMN_1730706"), ("HT12_R1", 7160239, "ILMN_1751607"),
                            ("HT12_R2", 7160239, "ILMN_1751607")], strict=True,
    ):
        mapping = manifest(directory / filename, digest)
        if mapping[fosb_address] != fosb_probe:
            raise ValueError("Incorrect FOSB address for manufacturer manifest")
        manifests[label] = (set(mapping), fosb_address)
    data = pg_table_blocks(raw_dump_bytes(input_bytes(directory, args.prefix)), RAW_TAG)
    report = census(raw_keys(data), selected, native_ids, manifests)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
