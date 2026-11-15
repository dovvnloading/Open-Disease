"""Restore a pinned seven-record ENA TSV when only portal response order varies."""
import argparse
from hashlib import sha256
from itertools import permutations
from pathlib import Path

from bulk_gse153739 import DIGESTS, ENA


def recover(payload):
    lines = payload.splitlines(keepends=True)
    if len(payload) != 1764 or len(lines) != 8 or not lines[0].startswith(b"run_accession\t"):
        raise ValueError("Unexpected ENA response")
    for rows in permutations(lines[1:]):
        possible = lines[0] + b"".join(rows)
        if sha256(possible).hexdigest() == DIGESTS[ENA]:
            return possible
    raise ValueError("Source records differ beyond ordering")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("response",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    verified=recover(args.response.read_bytes())
    if args.output.exists() and args.output.read_bytes()!=verified:
        raise ValueError("Existing output differs")
    args.output.write_bytes(verified)
