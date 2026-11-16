"""Read one inclusive interval from a sorted, position-indexed BGZF summary table."""

import gzip
import hashlib
import io
import struct
import urllib.request


def _read_struct(stream, fmt):
    packed = struct.Struct("<" + fmt)
    buf = stream.read(packed.size)
    if len(buf) != packed.size:
        raise ValueError("Truncated TBI structure")
    return packed.unpack(buf)


def _index(index_bytes, chromosome, start, end):
    """Return table columns, coordinate convention and virtual chunk boundaries."""
    if not (1 <= start <= end <= 2**29):
        raise ValueError("Invalid one-based inclusive TBI interval")
    stream = io.BytesIO(gzip.decompress(index_bytes))
    if stream.read(4) != b"TBI\1":
        raise ValueError("Expected tabix TBI index")
    refs, preset, seq, pos, stop, _, _, names_size = _read_struct(stream, "8i")
    names = stream.read(names_size).decode().rstrip("\x00").split("\x00")
    if len(names) != refs or chromosome not in names:
        raise ValueError(f"Chromosome {chromosome!r} absent from index: {names}")
    if preset & 0xFFFF or stop not in (0, pos) or min(seq, pos) < 1:
        raise ValueError("Only generic point-indexed summary tables are supported")
    wanted = {0}
    for offset, shift in ((1, 26), (9, 23), (73, 20), (585, 17), (4681, 14)):
        wanted.update(range(offset + ((start - 1) >> shift), offset + ((end - 1) >> shift) + 1))
    for ref in names:
        chunks = []
        for _ in range(_read_struct(stream, "i")[0]):
            bin_id, count = _read_struct(stream, "Ii")
            if ref == chromosome and bin_id in wanted:
                chunks.extend(_read_struct(stream, "QQ") for _ in range(count))
            else:
                stream.seek(16 * count, io.SEEK_CUR)
        n_linear = _read_struct(stream, "i")[0]
        if ref != chromosome:
            stream.seek(8 * n_linear, io.SEEK_CUR)
            continue
        i = (start - 1) >> 14
        minimum = 0
        if i < n_linear:
            stream.seek(8 * i, io.SEEK_CUR)
            minimum = _read_struct(stream, "Q")[0]
        hits = [(max(first, minimum), last) for first, last in chunks if last > minimum and last > first]
        extent = (min(x[0] for x in hits), max(x[1] for x in hits)) if hits else None
        return seq - 1, pos - 1, bool(preset & 0x10000), extent
    raise ValueError("Index lacks reference bins")


def _virtual_bytes(data, base, begin, end):
    """Decode BGZF blocks intersecting [begin, end) in virtual coordinates."""
    output = bytearray()
    cursor = 0
    while (base + cursor) << 16 < end:
        block = data[cursor:cursor + 18]
        if len(block) != 18 or block[:4] != b"\x1f\x8b\x08\x04" or block[12:16] != b"BC\x02\x00":
            raise ValueError("Invalid BGZF header")
        length = int.from_bytes(block[16:18], "little") + 1
        if cursor + length > len(data):
            raise ValueError("Truncated BGZF block")
        plain = gzip.decompress(data[cursor:cursor + length])
        addr = base + cursor
        low = begin & 65535 if addr == begin >> 16 else 0
        high = end & 65535 if addr == end >> 16 else len(plain)
        if not (0 <= low <= high <= len(plain)):
            raise ValueError("Invalid virtual BGZF subrange")
        output.extend(plain[low:high])
        cursor += length
    return bytes(output)


def regional_rows(url, index_bytes, chrom, start, end, max_bytes):
    """Return header fields, selected rows and receipts for bounded HTTP data ranges.

    max_bytes bounds wire payload for this call, excluding the caller's TBI index.
    """
    seq_col, pos_col, zero_based, extent = _index(index_bytes, chrom, start, end)
    receipts = []

    def fetch(first, last):
        if last - first + 1 + sum(r["bytes"] for r in receipts) > max_bytes:
            raise ValueError(f"BGZF range {first}-{last} exceeds {max_bytes}-byte call budget")
        req = urllib.request.Request(url, headers={"Range": f"bytes={first}-{last}", "Accept-Encoding": "identity"})
        with urllib.request.urlopen(req, timeout=45) as resp:
            content_range = resp.headers.get("Content-Range", "")
            if resp.status != 206 or not content_range.startswith(f"bytes {first}-"):
                raise ValueError(f"Server did not honor BGZF range {first}-{last}: {resp.status} {content_range}")
            span, total = content_range[6:].split("/")
            upper = int(span.split("-")[1])
            if upper != min(last, int(total) - 1):
                raise ValueError(f"Inconsistent BGZF Content-Range {content_range}")
            buf = resp.read(upper - first + 2)
            if len(buf) != upper - first + 1:
                raise ValueError("Truncated HTTP body")
        receipts.append({"range": content_range, "bytes": len(buf), "sha256": hashlib.sha256(buf).hexdigest()})
        return buf

    head = fetch(0, 65535)
    # The header line must fit in the first (possibly shorter) BGZF block.
    block_size = int.from_bytes(head[16:18], "little") + 1
    first_block = gzip.decompress(head[:block_size])
    if b"\n" not in first_block:
        raise ValueError("Missing complete header line in first block")
    header = first_block.split(b"\n", 1)[0].decode().rstrip("\r").lstrip("#").split("\t")
    if extent is None:
        return header, [], receipts
    begin, finish = extent
    first, last_block = begin >> 16, finish >> 16
    last = last_block + 65535 if finish & 65535 else last_block - 1
    buf = fetch(first, last)
    records = _virtual_bytes(buf, first, begin, finish).decode().splitlines()
    rows = []
    for record in records:
        fields = record.split("\t")
        if len(fields) != len(header):
            raise ValueError("Indexed row does not match header width")
        if fields[seq_col] == chrom and start <= int(fields[pos_col]) + int(zero_based) <= end:
            rows.append(fields)
    return header, rows, receipts
