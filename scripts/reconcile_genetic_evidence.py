"""Cross-check frozen nominations against public GWAS and ABO aggregate tables."""
import argparse
import csv
import hashlib
import io
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZipFile

from openpyxl import load_workbook


def read_table(workbook, name, header_row):
    rows = workbook[name].values
    for _ in range(header_row - 1):
        next(rows)
    header = next(rows)
    return [dict(zip(header, row)) for row in rows if row[0] is not None]


def ensembl(value):
    return value.split(".")[0]


def genetic_lookup(workbook, lock, crosswalk):
    with crosswalk.open() as source:
        mapping = {row["entrez_id"]: row for row in csv.DictReader(source, delimiter="\t")}
    families = {}
    for candidate in json.loads(lock.read_text())["candidates"]:
        ids = candidate.get("entrez_ids", [candidate["id"]])
        label = candidate.get("name") or mapping[candidate["id"]]["symbol"]
        families[label] = {mapping[str(i)]["ensembl_gene_id"] for i in ids}
    locked = set.union(*families.values())
    if len(locked) != len(mapping) or "" in locked:
        raise ValueError("Locked Ensembl mapping is incomplete or ambiguous")

    table8 = read_table(workbook, "Supplementary Table 8", 3)
    if not all(row["Posterior probability"] > .7
               and row["Posterior probability explained by SNP"] > .5 for row in table8):
        raise ValueError("Colocalization table contains rows below published thresholds")
    reported = {ensembl(r["Ensembl gene ID"]) for r in table8}

    table13 = read_table(workbook, "Supplementary Table 13", 3)
    twas = {ensembl(r["Gene transcript"]): r for r in table13}
    if len(twas) != len(table13):
        raise ValueError("MultiXcan table contains repeated Ensembl IDs")

    family_summary = {}
    for label, ids in families.items():
        observed = [twas[i] for i in sorted(ids & twas.keys())]
        passing = sorted(r["Gene name"] for r in observed if r["S MultiXcan p"] < 3.72e-6)
        family_summary[label] = {
            "locked_genes": len(ids), "multi_tissue_models": len(observed),
            "uterus_models": sum(isinstance(r["Uterus p-value"], (int, float)) for r in observed),
            "colocalized_ids": sorted(ids & reported), "multi_tissue_hits": passing,
            "minimum_multi_tissue_p": min((r["S MultiXcan p"] for r in observed), default=None),
        }
    abo = twas[mapping["28"]["ensembl_gene_id"]]
    return {
        "unique_locked_genes": len(locked), "colocalization_rows": len(table8),
        "colocalization_unique_genes": len(reported),
        "colocalization_unique_snps": len({r["Candidate SNP"] for r in table8}),
        "locked_genes_colocalized": sorted(reported & locked),
        "multi_tissue_genes_total": len(twas),
        "multi_tissue_hits_genomewide": sum(r["S MultiXcan p"] < 3.72e-6 for r in twas.values()),
        "multi_tissue_locked_genes_tested": len(locked & twas.keys()),
        "multi_tissue_locked_hits": sorted(twas[i]["Gene name"] for i in locked & twas.keys()
                                          if twas[i]["S MultiXcan p"] < 3.72e-6),
        "families": family_summary,
        "ABO_selected_tissue_models": {
            t: {"z": abo[f"{t} z-score"], "p": abo[f"{t} p-value"]}
            for t in ("Uterus", "Whole blood", "Breast mammary tissue", "Ovary", "Vagina",
                      "Sigmoid colon", "Transverse colon")
        },
    }


def blood_groups(data):
    with ZipFile(data / "elife-83116-table2-data1.zip") as archive:
        text = archive.read("table2_summary.tds").decode()
    rows = [r for r in csv.DictReader(io.StringIO(text), delimiter="\t")
            if r["phecode"] == "615.0" and r["phenotype"] == "Endometriosis"]
    if len(rows) != 5 or len({(r["pyrs"], r["cases"]) for r in rows}) != 1:
        raise ValueError("Danish endpoint/denominator mismatch")
    with ZipFile(data / "elife-83116-supp7.docx") as archive:
        document = ET.fromstring(archive.read("word/document.xml"))
    endpoint = []
    for row in document.findall(".//{*}tr"):
        cells = ["".join(c.itertext()).strip() for c in row.findall("{*}tc")]
        if cells[:2] == ["615", "Endometriosis"]:
            endpoint.append(cells)
    if len(endpoint) != 1 or int(endpoint[0][2]) != int(rows[0]["cases"]):
        raise ValueError("Danish A-versus-O row does not reconcile")
    fixed_reference = {}
    for group, i in zip(("A", "B", "AB", "RhD"), (4, 7, 10, 13)):
        match = re.fullmatch(r"([\d.]+) \(([\d.]+), ([\d.]+)\)", endpoint[0][i + 1])
        if match is None:
            raise ValueError(f"Invalid source interval for {group}")
        irr, lower, upper = map(float, match.groups())
        fixed_reference[group] = {"group_events": int(endpoint[0][i]), "irr": irr,
                                  "adjusted_95ci": [lower, upper],
                                  "adjusted_p": float(endpoint[0][i + 2])}

    uk = load_workbook(data / "uk.xlsx", read_only=True, data_only=True)
    rows_n80 = [r for r in uk["Table S2"].values if r[0] == "N80"]
    if len(rows_n80) != 1 or "Endometriosis" not in rows_n80[0][2]:
        raise ValueError("UK ICD-10 N80 row not uniquely identified")
    row = rows_n80[0]
    return {
        "denmark_female_phecode_615": {
            "total_events": int(rows[0]["cases"]), "person_years": float(rows[0]["pyrs"]),
            "ref_O_or_Rh_negative": fixed_reference,
            "ref_other_ABO_groups_or_Rh_negative": {
                r["blood_group"]: {"irr": float(r["estimate"]), "raw_p": float(r["p.value"]),
                                   "adjusted_p": float(r["p.value_FDR"])} for r in rows
            },
        },
        "ukb_mixed_sex_main_diagnosis_N80": {
            "cases": row[8], "controls": row[9],
            "probability_change_per_tag_allele": {
                name: {"beta": row[i], "p": row[i+1], "se": row[i+2]}
                for name, i in zip(("A1_rs507666_A", "A2_rs8176704_A", "B_rs8176746_T",
                                    "O_rs687289_G"), (10, 13, 16, 19))
            }
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.crosswalk.parent / "genetics_sources.json").read_text())
    for path in [args.lock, args.crosswalk] + [args.source_dir / s for s in manifest["download"]]:
        with path.open("rb") as source:
            digest = hashlib.file_digest(source, "sha256").hexdigest()
        if digest != manifest["sha256"][path.name]:
            raise ValueError(f"Input digest does not match: {path.name}")
    book = load_workbook(args.source_dir / "koller.xlsx", read_only=True, data_only=True)
    result = genetic_lookup(book, args.lock, args.crosswalk)
    result.update(blood_groups(args.source_dir))
    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
