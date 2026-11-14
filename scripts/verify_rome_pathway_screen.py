"""Run an independent base-R check from the public nomination lock and source data.

Usage: python scripts/verify_rome_pathway_screen.py data/private LOCK.json
The public lock supplies the exact membership of every nominated pathway. The separate
R implementation reads the original Affymetrix inputs and emits only aggregate
statistics for these candidates.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import tempfile

from analyze_gse25628_fosb import ANNOTATION_HASH, MATRIX_HASH
from open_disease.geo_microarray import sha256_file


def membership_rows(candidates):
    if not 1 <= len(candidates) <= 10:
        raise ValueError("Expected one to ten locked candidates")
    seen = set()
    for candidate in candidates:
        family, key, direction = (candidate[k] for k in ("family", "id", "direction"))
        if key in seen or direction not in ("higher", "lower"):
            raise ValueError("Duplicate nomination or unexpected direction")
        seen.add(key)
        if family == "gene":
            if not str(key).isdigit() or int(key) <= 0:
                raise ValueError("Nominated gene outside the prospective universe")
            members = (int(key),)
        elif family == "pathway" and str(key).startswith("R-HSA-"):
            members = candidate["entrez_ids"]
            if (len(members) != candidate["members"] or not 15 <= len(members) <= 200
                    or not all(type(gene) is int and gene > 0 for gene in members)
                    or list(members) != sorted(set(members))):
                raise ValueError("Invalid locked pathway membership")
        else:
            raise ValueError("Nominated family or pathway outside the prospective universe")
        if 2354 in members:
            raise ValueError("Previously studied FOSB is excluded")
        for gene in members:
            yield family, key, direction, gene


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", type=Path)
    parser.add_argument("lock", type=Path)
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    cc = json.loads(args.lock.read_text(encoding="utf-8"))["candidates"]
    rome = args.data / "gse25628"
    matrix = rome / "GSE25628_series_matrix.txt.gz"
    annotation = rome / "GPL571.annot.gz"
    if sha256_file(matrix) != MATRIX_HASH or sha256_file(annotation) != ANNOTATION_HASH:
        raise ValueError("Rome input SHA-256 differs from the audited snapshot")
    # The intermediate holds public gene-set definitions, never donor records.
    with tempfile.TemporaryDirectory() as folder:
        projection = Path(folder) / "membership.tsv"
        with projection.open("w", encoding="utf-8", newline="") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(("family", "id", "direction", "gene"))
            writer.writerows(membership_rows(cc))
        subprocess.run([
            args.rscript, str(Path(__file__).parent / "r/verify_rome_pathway_screen.R"),
            str(matrix), str(annotation),
            str(projection),
        ], check=True)


if __name__ == "__main__":
    main()
