"""Read checksum-pinned Pan-UK Biobank BGZF ranges to evaluate the Indian lead variants.

Local mode reads {rsID}.bgzfpart from the cache directory. --fetch permits exact,
bounded byte ranges of the public female phecode-615 aggregate GWAS file.
"""

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import struct
import urllib.request
from pathlib import Path

from india_gwas_portability import FAMILY, PANEL, load_indian_panel

SOURCE = "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/phecode-615-females.tsv.bgz"
HEADER_SHA = "243d0996649b6ddd260256a0ff14ccea4bc42515095f54924b3313459e53618a"


def read_part(path, start, end, digest, fetch):
    if start < 0 or not 0 < end - start + 1 <= 1_000_000:
        raise ValueError("Byte range outside limit")
    if path.exists():
        data = path.read_bytes()
    elif fetch:
        request = urllib.request.Request(
            SOURCE + "?block=" + str(start), headers={"Range": f"bytes={start}-{end}"}
        )
        with urllib.request.urlopen(request, timeout=30) as response:
            expected = f"bytes {start}-{end}/1588773197"
            if (
                response.status != 206
                or response.headers.get("Content-Range") != expected
            ):
                raise ValueError("Server did not honor the pinned byte range")
            data = response.read(end - start + 2)
    else:
        raise FileNotFoundError(path)
    if len(data) != end - start + 1 or hashlib.sha256(data).hexdigest() != digest:
        raise ValueError(f"Range length or checksum mismatch for {path.name}")
    if not path.exists():
        path.write_bytes(data)
    return data


def decompress_blocks(data):
    chunks = []
    pos = 0
    while pos + 18 <= len(data):
        if (
            data[pos : pos + 4] != b"\x1f\x8b\x08\x04"
            or data[pos + 12 : pos + 14] != b"BC"
        ):
            raise ValueError("Unexpected BGZF block header")
        length = struct.unpack_from("<H", data, pos + 16)[0] + 1
        if pos + length > len(data):
            break  # The requested range can end inside a BGZF block.
        chunks.append(gzip.decompress(data[pos : pos + length]))
        pos += length
    return b"".join(chunks).decode()


def read_site(header, text, chrom, pos, alleles):
    rows = [
        fields
        for fields in csv.reader(io.StringIO(text), delimiter="\t")
        if len(fields) == len(header)
        and fields[0] == str(chrom)
        and fields[1].isdigit()
    ]
    positions = [int(r[1]) for r in rows]
    if not positions or not min(positions) <= pos <= max(positions):
        raise ValueError("Pinned range does not cover requested coordinate")
    hits = [
        dict(zip(header, r, strict=True))
        for r in rows
        if int(r[1]) == pos and {r[2], r[3]} == alleles
    ]
    if len(hits) != 1:
        raise ValueError(f"Expected one position-and-allele match for {chrom}:{pos}")
    return hits[0]


def analyze(supplement, ranges, cache, fetch):
    discovery = load_indian_panel(supplement)
    targets = json.loads(ranges.read_text())
    if [t["rsid"] for t in targets] != PANEL:
        raise ValueError("Range manifest differs from the full discovery family")
    cache.mkdir(parents=True, exist_ok=True)
    header_part = read_part(cache / "panukb_header", 0, 65535, HEADER_SHA, fetch)
    header = decompress_blocks(header_part).splitlines()[0].split("\t")
    result = []
    for site in targets:
        name = site["rsid"]
        panel = discovery[name]
        data = read_part(
            cache / (name + ".bgzfpart"),
            site["start"],
            site["end"],
            site["sha256"],
            fetch,
        )
        row = read_site(
            header, decompress_blocks(data), site["chr"], site["pos"], panel["alleles"]
        )
        p = 10 ** (-float(row["neglog10_pval_CSA"]))
        out = dict(
            rsid=name,
            risk_allele=panel["risk"],
            two_sided_p=p,
            bonferroni_p_21=min(1, p * FAMILY),
        )
        sign = 1 if panel["risk"] == row["alt"] else -1
        status = "aligned"
        if row["low_confidence_CSA"] != "false":
            status = "low_confidence"
        elif panel["alleles"] in ({"A", "T"}, {"C", "G"}):
            freq = float(row["af_controls_CSA"])
            aligned_freq = freq if sign == 1 else 1 - freq
            if (
                abs(panel["india_freq"] - 0.5) < 0.1
                or abs(aligned_freq - panel["india_freq"]) > 0.1
            ):
                status = "strand_ambiguous"
        out["status"] = status
        if status == "aligned":
            beta = sign * float(row["beta_CSA"])
            se = float(row["se_CSA"])
            out.update(
                odds_ratio=math.exp(beta),
                lower_95=math.exp(beta - 1.96 * se),
                upper_95=math.exp(beta + 1.96 * se),
            )
        result.append(out)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("supplement", type=Path)
    parser.add_argument("ranges", type=Path)
    parser.add_argument("cache", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    rows = analyze(args.supplement, args.ranges, args.cache, args.fetch)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(dict.fromkeys(k for r in rows for k in r)),
            delimiter="\t",
            lineterminator="\n",
            restval="NA",
        )
        writer.writeheader()
        writer.writerows(rows)
