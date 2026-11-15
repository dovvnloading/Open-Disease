"""Verify array reuse with CEL intensity measurements across text, XDA and Calvin files.

Inputs are locally cached public study supplement prefixes plus published classifier
tables. Report only aggregate counts. A match requires the float32 mean intensities
of the first 100 physical grid cells in original order; separately verify native
DAT header and grid corners wherever available. No participant linkage is used.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import struct
import zlib

COUNT = 100
SERIES = {"1": "GSE4888", "2": "GSE6364", "3": "GSE51981", "4": "GSE29981"}

INPUT_HASHES = {
    "coldata.tsv": "a2696216f28f6f80a8e7db716f280c45e4d2e4264730723fb3775b94751ebc4a",
    "test.tsv": "35cbd873dabf6f5567f9a01a2c6fad29850a2920597328edeac762cca1f4f273",
}
SERIES_HEADER_HASHES = {
    "GSE4888": "86555a082e3322982b17f3f595e416b56b3cea796225e9bf7bd4b9b080c534b4",
    "GSE6364": "73928e6928400b48b192db0ccf058374f77f36c123f0c4383002cde062da7167",
    "GSE51981": "73d9dc597c10563744e06d02a7800bce7336e302b1f8948d870e8837c4f02909",
    "GSE29981": "f002d1fb7980b735a23a9c166c5faf684304bcc34b088f0137b1279a53ab8933",
}


def hash_files(paths):
    combined = hashlib.sha256()
    for path in sorted(paths):
        combined.update(path.read_bytes())
    return combined.hexdigest()


def counted_string(buf, off, encoding):
    size, = struct.unpack_from(">I", buf, off)
    end = off + 4 + size * 2
    return buf[off + 4:end].decode(encoding), end


def cell_vector(raw):
    if raw.startswith(b"[CEL]"):
        headers = raw.split(b"[INTENSITY]")[0]
        lines = raw.split(b"[INTENSITY]")[1].splitlines()
        j = next(j for j, line in enumerate(lines) if line.startswith(b"CellHeader="))
        if lines[j].split(b"=", 1)[1].split() != [b"X", b"Y", b"MEAN", b"STDV", b"NPIXELS"]:
            raise ValueError("Unexpected text CEL columns")
        means = []
        for ix, row in enumerate(lines[j + 1:j + 1 + COUNT]):
            x, y, intensity, _, _ = row.split()
            if int(x) != ix or int(y) != 0:
                raise ValueError("Unexpected CEL grid order")
            means.append(float(intensity))
        file_format = "text"
    elif raw[:4] == struct.pack("<I", 64):
        _, ver, rows, cols, cells = struct.unpack_from("<5I", raw)
        if ver != 4 or cells != rows * cols or rows < COUNT:
            raise ValueError("Unexpected XDA CEL layout")
        off = 20
        for i in range(3):
            size, = struct.unpack_from("<I", raw, off)
            if i == 0:
                headers = raw[off + 4:off + 4 + size]
            off += 4 + size
        off += 16
        means = [x[0] for x in struct.iter_unpack("<ffh", raw[off:off + COUNT * 10])]
        file_format = "XDA"
    else:
        if raw[:2] != b"\x3b\x01":
            raise ValueError("Unrecognized CEL file")
        group, = struct.unpack_from(">I", raw, 6)
        dataset, = struct.unpack_from(">I", raw, group + 4)
        data_off, = struct.unpack_from(">I", raw, dataset)
        name, k = counted_string(raw, dataset + 8, "utf-16be")
        if name != "Intensity":
            raise ValueError("Expected a Calvin Intensity dataset")
        parameters, columns = struct.unpack_from(">II", raw, k)
        colname, k = counted_string(raw, k + 8, "utf-16be")
        kind, width, rows = struct.unpack_from(">BII", raw, k)
        if (parameters, columns, colname, kind, width) != (0, 1, "Intensity", 6, 4) or rows < COUNT:
            raise ValueError("Calvin channel is not one column of floating-point mean intensities")
        means = struct.unpack_from(f">{COUNT}f", raw, data_off)
        headers = None
        file_format = "Calvin"
    if len(means) != COUNT:
        raise ValueError("Incomplete CEL measurement prefix")
    digest = hashlib.sha256(struct.pack(f">{COUNT}f", *means)).digest()
    key = None
    if headers:
        h = dict(row.split(b"=", 1) for row in headers.splitlines() if b"=" in row)
        key = tuple(h[k] for k in (b"DatHeader", b"GridCornerUL", b"GridCornerUR",
                                    b"GridCornerLR", b"GridCornerLL"))
    return digest, key, file_format


def verify(directory, calvin):
    if any(hash_files([directory / name]) != digest for name, digest in INPUT_HASHES.items()):
        raise ValueError("Published classifier table checksum mismatch")
    if any(hash_files((directory / study).glob("*.part")) != digest
           for study, digest in SERIES_HEADER_HASHES.items()):
        raise ValueError("Native CEL prefix checksum mismatch")
    if hash_files(calvin.glob("*/*.part")) != "5b73bf809a819992c14b6d7b7421875d4ef664cff72b9aa07583843b6b336e3e":
        raise ValueError("Extended Calvin prefix checksum mismatch")

    def read(filename):
        with (directory / filename).open(newline="") as file:
            return list(csv.DictReader(file, delimiter="\t"))
    mapping = read("coldata.tsv")
    predictions = read("test.tsv")
    test = {r["sample_name"]: r for r in predictions}
    if len(mapping) != 216 or len(test) != len(predictions) or len(test) != 37:
        raise ValueError("Unexpected classifier cohort size")
    taxa = {}
    for row in mapping:
        sample = row["sample_name"]
        study = SERIES[row["batch"]]
        if sample in taxa or (sample in test) != (study == "GSE6364"):
            raise ValueError("Published test membership or identifier is ambiguous")
        if sample in test and int(test[sample]["true_label"]) != int(row["condition"] == "endometriosis"):
            raise ValueError("Author diagnosis does not match the published prediction file")
        part = directory / study / (sample + ".part")
        raw = zlib.decompressobj(31).decompress(part.read_bytes())
        if raw[0] == 59:
            longer = (calvin / study / part.name).read_bytes()
            if len(longer) != 8192 or not longer.startswith(part.read_bytes()):
                raise ValueError("Mismatched Calvin CEL compressed prefix")
            raw = zlib.decompressobj(31).decompress(longer)
        fingerprint, key, file_format = cell_vector(raw)
        taxa[sample] = (fingerprint, key, file_format, study, row["phase"])
    train = {i: v for i, v in taxa.items() if i not in test}
    scope = {}
    for group in ["all", "proliferative", "early_secretory", "mid_secretory"]:
        tally = Counter()
        by_source = Counter()
        for i, pred in test.items():
            f, k, fmt, _, phase = taxa[i]
            if group != "all" and phase != group:
                continue
            hit = [v for v in train.values() if v[0] == f]
            header_hit = [v for v in train.values() if k is not None and v[1] == k]
            if Counter(hit) != Counter(header_hit):
                raise ValueError("Intensity and DAT/corner matches disagree")
            tally["tested"] += 1
            tally["correct"] += pred["true_label"] == pred["predicted_label"]
            tally["intensity_matched"] += bool(hit)
            tally["unmatched_prefixes"] += not hit
            tally["matched_training_arrays"] += len(hit)
            for study in {v[3] for v in hit}:
                by_source[study] += 1
        scope[group] = {**tally, "test_matches_by_training_study": dict(sorted(by_source.items()))}
    return dict(n_train=len(train), n_test=len(test), physical_grid_cells=COUNT,
                formats_in_train=dict(Counter(v[2] for v in train.values())), results=scope)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("calvin", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.directory, args.calvin), indent=2))
