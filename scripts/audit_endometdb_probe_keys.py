"""Inspect probe annotations and native assay keys in two pinned EndometDB prefixes.

The local inputs, offsets and COPY layouts are specific to Figshare v1. Only
complete COPY lines contribute, with quantitative fields passed over as bytes.
The output contains aggregate counts and manufacturer probe identifiers.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import itertools
import json
from pathlib import Path
import zlib

from audit_endometdb import HEAD_SHA256, eligible, pinned, read_metadata


RAW_RANGE = "90cbe73c4fd51f102f73a7bfdf834597c0bf7f6b50eb9c2f8ec577d2abdd58e6"
NORMALIZED_RANGE = "a62f3374bf39d6b460994089298fec0a90dc8049757c4565504c782c063c129a"
MANIFESTS = {
    "GPL13376-platform.gz": "c6f2c87ab7bc094b22270c907e1f9a6d74e9990b62726b53bf44660b0db700c1",
    "GPL10558-R1.gz": "770fb74bb90e548ee731924310d8f5c48589e63ee530dc11c7686d28a3ba5ba6",
    "GPL10558-R2.gz": "932e2f110bd30ac70ef4d7dbac2f87f8088f982bf9aa596732c7477ab50d43c2",
}
FOSB_ADDRESSES = (6040086, 7160239)


def copy_prefix(buffer, offset, block_tag):
    """Unpack PG's length-prefixed zlib stream, dropping the final incomplete line."""
    if buffer[offset:offset + 6] != block_tag:
        raise ValueError("Unexpected pg_dump table header")
    offset += 6
    compressed = bytearray()
    while offset + 5 <= len(buffer):
        size = int.from_bytes(buffer[offset + 1:offset + 5], "little")
        if buffer[offset] != 0 or not 0 < size <= 4096:
            raise ValueError("Unexpected pg_dump chunk in table prefix")
        offset += 5
        compressed.extend(buffer[offset:offset + size])
        offset += size
    data = zlib.decompressobj().decompress(compressed)
    return data.split(b"\n")[:-1]


def manifest(path, sha256):
    pinned(path, sha256)
    with gzip.open(path, "rt", encoding="utf-8") as source:
        for line in source:
            if "Probe_Id\tArray_Address_Id\t" in line:
                break
        else:
            raise ValueError("Missing manufacturer annotation header")
        addresses = {}
        for row in csv.DictReader(source, fieldnames=line.rstrip().split("\t"), delimiter="\t"):
            address = row.get("Array_Address_Id")
            if address and address.isdecimal():
                address_id = int(address)
                if address_id in addresses:
                    raise ValueError("Nonunique address in manufacturer manifest")
                addresses[address_id] = row["Probe_Id"]
        return addresses


def raw_summary(lines, wg6, ht12_addresses, assay_ids, eligible_ids):
    """Only an assay with a complete unique manifest address set verifies the chip."""
    def assay_and_address(line):
        _, assay, probe, _ = line.split(b"\t", 3)
        return assay.decode(), int(probe)

    groups = itertools.groupby(map(assay_and_address, lines), key=lambda pair: pair[0])
    verified = set()
    seen = set()
    raw_rows = 0
    for assay, group in groups:
        if assay not in assay_ids or assay in seen:
            raise ValueError("Unknown or non-contiguous native assay in raw prefix")
        seen.add(assay)
        addresses = [address for _, address in group]
        distinct = set(addresses)
        verified_chip = len(addresses) == len(distinct) and distinct == wg6.keys()
        if verified_chip:
            verified.add(assay)
        raw_rows += len(addresses)
    return {
        "native_assays_in_prefix": len(seen),
        "complete_wg6_assays": len(verified),
        "complete_wg6_eligible_assays": len(verified & eligible_ids),
        "raw_rows_in_prefix": raw_rows,
        "wg6_distinct_addresses": len(wg6),
        "wg6_addresses_absent_from_ht12": len(wg6.keys() - ht12_addresses),
    }


def normalized_summary(lines):
    counts = Counter()
    for line in lines:
        assay, _, search, symbol, probe, address = line.split(b"\t")
        counts[(assay, search, symbol, probe, address)] += 1
    first_probe = next(iter(counts))[3]
    first_probe_counts = Counter()
    for (assay, _, _, probe, _), n in counts.items():
        if probe == first_probe:
            first_probe_counts[assay] += n
    return {
        "rows": sum(counts.values()),
        "distinct_assays": len({k[0] for k in counts}),
        "distinct_reported_probes": len({k[3] for k in counts}),
        "distinct_keys_excluding_signal": len(counts),
        "zero_search_key_and_address_rows": sum(n for k, n in counts.items() if k[1] == k[4] == b"0"),
        "fosb_symbol_rows": sum(n for k, n in counts.items() if k[2] == b"FOSB"),
        "first_reported_probe": first_probe.decode(),
        "first_probe_rows": sum(first_probe_counts.values()),
        "first_probe_distinct_assays": len(first_probe_counts),
        "first_probe_rows_per_assay_range": [min(first_probe_counts.values()), max(first_probe_counts.values())],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Local directory containing the documented pinned inputs")
    directory = parser.parse_args().directory
    from inflate64 import Inflater

    head = directory / "archive-start.bin"
    clinical = read_metadata(head, directory / "archive-tail.bin")
    _, selected = eligible(clinical)
    native_ids = {row["id"] for row in clinical["sample_analysis_result"]}
    eligible_ids = {row["analysis_result_id"] for row in selected}
    raw = Inflater().inflate(pinned(head, HEAD_SHA256)[48:] +
                             pinned(directory / "native-raw-continuation.bin", RAW_RANGE))
    fk = (b"FOREIGN KEY (sample_analysis_result_id) REFERENCES public.sample_analysis_result(id)")
    if b"ADD CONSTRAINT sample_analysis_result_microarray_data_fk " + fk not in raw[:38714]:
        raise ValueError("Missing native raw-to-assay foreign key")
    normalized_range = pinned(directory / "native-normalized-head-range.bin", NORMALIZED_RANGE)
    normalized = Inflater().inflate(normalized_range[2103069:])
    wg6, ht12_r1, ht12_r2 = [manifest(directory / filename, digest) for filename, digest in MANIFESTS.items()]
    print(json.dumps({
        "address_to_probe": {name: {address: m[address] for address in FOSB_ADDRESSES}
                             for name, m in zip(MANIFESTS, (wg6, ht12_r1, ht12_r2), strict=True)},
        "raw": raw_summary(copy_prefix(raw, 38714, bytes.fromhex("0100d8080000")), wg6,
                           ht12_r1.keys() | ht12_r2.keys(), native_ids, eligible_ids),
        "normalized": normalized_summary(copy_prefix(normalized, 92118, bytes.fromhex("0100da080000"))),
    }, indent=2))


if __name__ == "__main__":
    main()
