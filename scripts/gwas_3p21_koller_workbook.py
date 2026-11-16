"""Audit locus-specific rows of Koller supplementary tables 4, 5, 7, and 8."""

import argparse
import hashlib
import json
import re
from pathlib import Path

from openpyxl import load_workbook


SHA256 = "c021b6c20c668644939fd83381307a2c8cf08f18d95156b33e65fb0301997815"
WORKBOOK_URL = (
    "https://media.springernature.com/full/springer-static/esm/art%3A10.1038"
    "%2Fs41588-026-02582-2/MediaObjects/41588_2026_2582_MOESM4_ESM.xlsx"
)
START, END = 48_652_148, 50_652_148  # GRCh37 chr3; frozen design.md
MAPPINGS = {
    "rs9831967": (50_217_234, "C/A/G/T", "2a7d09473dd3e98b1bf81aeaffed6d86a13bf100998a1c6d278a0687355b513f"),
    "rs6804970": (49_804_290, "C/T", "fbafd1e0d2cb0757d7c33813ae7a305bb7644320c93be20021a47f510827dd1a"),
}


def overlaps(chromosome, start, end):
    if not (0 < start <= end):
        raise ValueError("Expected an ascending, positive interval")
    return str(chromosome).removeprefix("chr") == "3" and start <= END and end >= START


def extract_table(sheet, number):
    selected = []
    scanned = 0
    section = None
    # S4 has a separate heading for hits from the non-cross-ancestry GWAS at row 87.
    data_start = 6 if number == 4 else 4
    for row_no, row in enumerate(sheet.iter_rows(min_row=data_start, values_only=True), data_start):
        if row[1] is None:
            continue
        scanned += 1
        if number == 8:
            # S8 reports chromosome and candidate rsID, without variant positions. Retain
            # every chromosome 3 row so a gene-name filter cannot hide other candidates.
            if str(row[4]) == "3":
                selected.append(dict(zip(
                    ["row", "analysis", "tissue", "gene_id", "gene_symbol", "chromosome",
                     "p_gwas_as_labelled", "posterior_probability", "regional_probability",
                     "candidate_snp", "posterior_proportion_explained_by_snp"],
                    [row_no, *row[:10]], strict=True,
                )))
            continue

        chromosome, position = row[1:3]
        in_interval = overlaps(chromosome, position, position)
        if number == 4:
            if in_interval:
                effects = {}
                for analysis, offset in (
                    ("cross_ancestry_combined", 7), ("eur_combined", 27), ("eur_clinical", 32),
                    ("eur_self_reported", 37), ("eur_adenomyosis", 42), ("eur_excluding_adenomyosis", 47),
                ):
                    effects[analysis] = dict(zip(
                        ["beta", "se", "p", "eaf", "heterogeneity_p"], row[offset:offset + 5], strict=True,
                    ))
                selected.append({
                    "row": row_no, "variant": row[0], "position": position, "effect_allele": row[3],
                    "other_allele": row[4], "nearest_gene": row[5], "classification": row[6], "effects": effects,
                })
        else:
            # S5 and S7 print the section label only on its first data row.
            section = row[0] or section
            if not section:
                raise ValueError(f"S{number} analysis label missing before row {row_no}")
            if number == 5 and in_interval:
                selected.append({
                    "row": row_no, "analysis": section, "variant": row[3], "position": position,
                    "beta": row[6], "se": row[7], "p": row[8], "joint_beta": row[9],
                    "joint_se": row[10], "joint_p": row[11],
                })
            if number == 7:
                region = re.fullmatch(r"chr([\w]+):(\d+)-(\d+)", row[8])
                if not region or region[1] != str(chromosome).removeprefix("chr"):
                    raise ValueError(f"S7 locus/chromosome mismatch at row {row_no}")
                locus_overlaps = overlaps(region[1], int(region[2]), int(region[3]))
                if in_interval or locus_overlaps:
                    selected.append({
                        "row": row_no, "analysis": section, "position": position, "variant": row[3],
                        "allele_1": row[4], "allele_2": row[5], "z_score": row[6], "pip": row[7],
                        "source_locus": row[8], "marker_in_interval": in_interval, "locus_overlaps": locus_overlaps,
                    })
    return {"data_rows_scanned": scanned, "rows": selected}


def audit(path):
    with path.open("rb") as src:
        digest = hashlib.file_digest(src, "sha256").hexdigest()
    if digest != SHA256:
        raise ValueError("Koller workbook SHA-256 differs from the final publication pin")
    book = load_workbook(path, read_only=True, data_only=True)
    try:
        tables = {f"S{n}": extract_table(book[f"Supplementary Table {n}"], n) for n in (4, 5, 7, 8)}
    finally:
        book.close()
    return {
        "source": WORKBOOK_URL, "source_sha256": digest, "source_bytes": path.stat().st_size,
        "assembly": "GRCh37", "chromosome": "3", "start_inclusive": START, "end_inclusive": END,
        "position_crosschecks": {
            marker: {
                "url": f"https://grch37.rest.ensembl.org/variation/human/{marker}?content-type=application/json",
                "retrieval": "GET with Accept: application/json, 2026-11-16", "receipt_sha256": sha,
                "chromosome": "3", "grch37_position": pos, "ensembl_allele_string": alleles,
            } for marker, (pos, alleles, sha) in MAPPINGS.items()
        },
        "selection": {
            "S4": "All published lead rows were examined, including the separate non-cross-ancestry section. Report leads in the frozen interval.",
            "S5": "All reported genome-wide significant secondary association rows were examined. Absence from this table does not establish a zero conditional effect.",
            "S7": "Report either a marker in the frozen interval or any published PAINTOR interval overlapping it. The sheet reports selected variants above 0.5, not a complete credible set.",
            "S8": "Report every published positive row on chromosome 3. The SNP coordinates are absent from S8. The p gwas field is transcribed as labelled, without assigning a new meaning or inferring an allele.",
        },
        **tables,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path, help="Local copy of the final Koller workbook")
    print(json.dumps(audit(parser.parse_args().workbook), indent=2))


if __name__ == "__main__":
    main()
