"""Audit published 2023 workbook cells. Takes a local file, never downloads data."""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

LEAD = "rs1352889"
CANDIDATE = "rs13100173"
SHA256 = "01cc0de4f3cb9aeec523ba312a17b931606d627e17e83581700b559741c7585e"
MARKERS = {
    LEAD: (49652148, 49614715),
    CANDIDATE: (50332697, 50295266),
}
# NCBI Gene cross-references for HYAL3 (8372), RBM6 (10180),
# TMEM115 (11070) and GMPPB (29925).
GENES = {
    "ENSG00000186792": "HYAL3",
    "ENSG00000004534": "RBM6",
    "ENSG00000126062": "TMEM115",
    "ENSG00000173540": "GMPPB",
}


def analyze_uk(sheets):
    """Count both coordinate fields independently, preserving the source column labels."""
    hits = {marker: [] for marker in MARKERS}
    coordinate_counts = {marker: Counter() for marker in MARKERS}
    genes = Counter()
    gene_checks = {}
    for sheet, rows in sheets:
        for line, row in enumerate(rows, 1):
            if line < 4 or row[8] != 3:
                continue
            gene, variant_id, start = row[0], str(row[7]), row[9]
            if gene in GENES:
                value = row[11]
                summary = gene_checks.setdefault(GENES[gene], {
                    "n_cis_rows": 0, "column_L_min": value,
                    "column_L_max": value, "flagged_best": [],
                })
                summary["n_cis_rows"] += 1
                summary["column_L_min"] = min(summary["column_L_min"], value)
                summary["column_L_max"] = max(summary["column_L_max"], value)
                if row[13] == 1:
                    summary["flagged_best"].append({
                        "sheet": sheet, "row": line, "id": variant_id, "column_L": value,
                    })
            if variant_id.startswith("3:") and variant_id[2:].isdigit():
                if 48652148 <= int(variant_id[2:]) <= 50652148:
                    genes[gene] += 1
            for marker, (gr37, gr38) in MARKERS.items():
                by_id = variant_id in (f"3:{gr37}", marker)
                by_position = start == gr38
                if not (by_id or by_position):
                    continue
                coordinate_counts[marker]["matched_id"] += int(by_id)
                coordinate_counts[marker]["matched_position"] += int(by_position)
                coordinate_counts[marker]["matched_both"] += int(by_id and by_position)
                hits[marker].append({
                    "gene_id": gene, "sheet": sheet, "row": line,
                    "cis_variant_count_column": row[5],
                    "column_L_as_reported": row[11], "slope_unoriented": row[12],
                })

    marker_summary = {}
    for marker, matches in hits.items():
        ratios = [entry["slope_unoriented"] ** 2 / entry["column_L_as_reported"]
                  for entry in matches if entry["column_L_as_reported"] > 0]
        marker_summary[marker] = {
            k: coordinate_counts[marker][k]
            for k in ("matched_id", "matched_position", "matched_both")
        }
        marker_summary[marker].update({
            "distinct_gene_ids": len({entry["gene_id"] for entry in matches}),
            "slope_squared_divided_by_L_range": [min(ratios), max(ratios)] if ratios else None,
            "prespecified_genes": {
                symbol: [entry for entry in matches if entry["gene_id"] == gene_id]
                for gene_id, symbol in GENES.items()
            },
        })
    return {
        "N_donors": 163,
        "marker_rows": marker_summary,
        "rows_in_audit_interval_by_numeric_variant_id": dict(sorted(genes.items())),
        "all_cis_column_L_checks": gene_checks,
        "caveats": [
            "The numeric variant ID and variant position columns use different coordinates at the BSN lead.",
            "Neither sheet reports REF, ALT, dosage allele, frequency or standard error. No risk-allele effect is inferred.",
            "Column L has the nominal P-value heading. Its proportionality to squared slopes across genes at the lead makes its interpretation as a P-value unsafe.",
        ],
    }


def extract(path, members_out=None):
    from openpyxl import load_workbook

    with path.open("rb") as file:
        digest = hashlib.file_digest(file, "sha256").hexdigest()
    if digest != SHA256:
        raise ValueError("Input differs from the verified 2023 publisher workbook")
    workbook = load_workbook(path, data_only=True, read_only=True)
    try:
        def table(n):
            sheet = workbook.worksheets[n - 1]
            if not sheet.title.startswith(f"Supp{n}."):
                raise ValueError(f"Unexpected sheet at position {n}")
            return sheet.values

        s3 = [r for r in table(3) if r[2] == LEAD]
        s10 = [r for r in table(10) if r[1] == LEAD]
        (lead,) = s3
        (summary,) = s10
        annotated = [r for r in table(11) if r[0] == LEAD]
        members = {r[2] for r in annotated}
        bounds = [min(r[5] for r in annotated), max(r[5] for r in annotated)]
        if len(members) != summary[4] or bounds != list(map(int, summary[6].split("-"))):
            raise ValueError("S10 and S11 disagree on credible-set membership or extent")

        encoded_members = ("\n".join(sorted(members)) + "\n").encode()
        if members_out is not None:
            members_out.write_bytes(encoded_members)
        variants = {f"chr3:{positions[0]}" for positions in MARKERS.values()}
        return {
            "doi": "10.1038/s41588-023-01323-z",
            "workbook_sha256": digest, "bytes": path.stat().st_size,
            "S3_overall": {"lead": lead[2], "hg19": lead[3], "risk_allele_EUR_JPT_frequency": lead[4],
                           "odds_ratio_CI": lead[5], "p_as_printed": lead[6], "distinct_signals": lead[13]},
            "S8_BSN_rows": sum(r[1] == LEAD or r[0] == "BSN/3p21.31" for r in table(8)),
            "S10_S11_99pct_credible_set": {
                "index": LEAD, "unique_members": len(members), "hg19_min_max": bounds,
                "sorted_member_ids_lf_sha256": hashlib.sha256(encoded_members).hexdigest(),
                "lead_is_member": LEAD in members, "HYAL3_candidate_is_member": CANDIDATE in members,
                "named_cross_study_markers": {
                    marker: marker in members
                    for marker in ("rs9831967", "rs6804970", "rs6446298", "rs6792892")
                },
                "bp_candidate_beyond_largest_member": MARKERS[CANDIDATE][0] - bounds[1],
                "individual_posteriors_in_S11": False,
            },
            "S12_coding_members": [{"variant": r[2], "gene": r[11], "consequence": r[8]}
                                   for r in table(12) if r[0] == LEAD],
            "S13_high_confidence_BSN_rows": sum(r[1] == LEAD for r in table(13)),
            "S31_chr3_selected_endometrial_eqtl_rows": [
                {"variant": r[0], "position_as_printed": r[2], "allele1": r[3],
                 "allele2": r[4], "allele1_frequency": r[5], "beta": r[7],
                 "se": r[8], "p": r[9], "gene": r[10]} for r in table(31) if r[1] == 3
            ],
            "S33_selected_top10000_rows_at_two_markers": [
                {"variant_id": r[0], "effect_allele": r[3], "other_allele": r[4],
                 "frequency": r[5], "beta": r[6], "se": r[7], "p": r[8], "sample_size_column": r[9]}
                for r in table(33) if r[0] in variants
            ],
            "S34_S35_UK_endometrial_QTL": analyze_uk((n, table(n)) for n in (34, 35)),
        }
    finally:
        workbook.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("--members-out", type=Path, help="Optional local text file for sorted S11 member IDs")
    args = parser.parse_args()
    print(json.dumps(extract(args.workbook, args.members_out), indent=2, allow_nan=False))
