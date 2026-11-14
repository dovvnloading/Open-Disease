"""Verify the public epithelial derivative and run its aggregate-only base-R audit."""
from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import subprocess
import tempfile

from open_disease.geo_microarray import sha256_file

DIGEST = "50715d0e2f5a09b8f4075668569a75baff9ce24bf8e62503f259858bb1c7d683"
SIZE = 199_996_542


def verify(path):
    if path.stat().st_size != SIZE or sha256_file(path) != DIGEST:
        raise ValueError("RDS length or SHA-256 differs from the pinned GEO deposit")


WORKBOOK_DIGEST = "21d7996ad209662a550fceb9694b3e9b70eb1a7e9cb44cce26275a60a10891dd"
SOURCE_COLUMNS = ("Patient sample ID", "endometriosis status", "MenstrualCyclePhase",
                  "10x capturing batch", "sequencing batch", "DEG analysis (Figure 3a,b)")


def write_sample_map(workbook, output):
    from openpyxl import load_workbook

    if sha256_file(workbook) != WORKBOOK_DIGEST:
        raise ValueError("Unrecognized Supplementary Data 1")
    book = load_workbook(workbook, read_only=True, data_only=True)
    try:
        records = iter(book.active.values)
        header = next(records)
        if any(header.count(key) != 1 for key in SOURCE_COLUMNS):
            raise ValueError("Unexpected clinical table headers")
        selected = [header.index(key) for key in SOURCE_COLUMNS]
        with output.open("w", encoding="utf-8", newline="") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            writer.writerow(("specimen", "diagnosis", "phase", "capture", "sequencing", "deg"))
            for row in records:
                if row[selected[0]] is not None:
                    writer.writerow([row[i] for i in selected])
    finally:
        book.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("rds", type=Path)
    parser.add_argument("supplement", type=Path, help="Published Supplementary Data 1 XLSX")
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    verify(args.rds)
    env = dict(os.environ, R_DEFAULT_PACKAGES="NULL", OPENBLAS_NUM_THREADS="1")
    script = Path(__file__).parent / "r/audit_gse266265_epithelial.R"
    # Temporary exact joins stay inside this private directory until the R
    # process exits. TemporaryDirectory grants access only to the current user.
    with tempfile.TemporaryDirectory(dir=args.rds.parent) as folder:
        table = Path(folder) / "study-design.tsv"
        write_sample_map(args.supplement, table)
        subprocess.run([args.rscript, str(script), str(args.rds), str(table)], env=env, check=True)
