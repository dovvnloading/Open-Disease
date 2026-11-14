"""Independently check the HECA workbook census by relational joins.

Use: python scripts/verify_heca_glandular.py data/private/heca/41588_2024_1873_MOESM3_ESM.xlsx
Reads Tables 1–3 and prints only aggregate results. No expression inputs.
"""
from __future__ import annotations

import argparse
import csv
from itertools import combinations
from collections import Counter, defaultdict
from hashlib import file_digest
import json
from pathlib import Path
import sqlite3

from openpyxl import load_workbook

WORKBOOK_SHA256 = "f70030404af3f6070c10f22a064e20188b0fb3251c1042605563c0aea00d2ed3"


def unique_clinical_key(donor, study):
    """Use exact published naming conventions for the Table 1 to Table 3 join."""
    study = {"Huangb": "Huang", "Garcia-Alonso": "GarciaAlonso",
             "Mareckova-cells": "Mareckova", "Mareckova-cells & nuclei": "Mareckova"}.get(study, study)
    if study == "Lai":
        donor = donor.removesuffix("a")  # footnote a is attached directly to one cell
    return study, donor if study == "Fonseca" else f"{study}_{donor}"


def verify_anchor_libraries(path, sc_libraries, clinical, eligible):
    """Demultiplexed donor assignments must retain one diagnosis across repeated captures."""
    with path.open("rb") as source:
        if file_digest(source, "sha256").hexdigest() != "cdcd25f5519b3ddac7b68a41140ae577b929a856c35732e598fc95fcb0bbf36b":
            raise ValueError("Unrecognized E-MTAB-14039 SDRF checksum")
    pairs = set()
    library_donors = defaultdict(set)
    with path.open(encoding="utf-8", newline="") as handle:
        for record in csv.DictReader(handle, delimiter="\t"):
            lib = record["Characteristics[library_id]"]
            if lib not in sc_libraries:
                continue
            key = "Mareckova_" + record["Characteristics[individual]"]
            group = {"endometriosis": "case", "normal": "control"}[record["Characteristics[disease]"]]
            if (key, lib) in pairs or clinical[key] != (group, record["Characteristics[menstrual_cycle_stage]"]):
                raise ValueError("Repeated sc donor-library link or inconsistent phenotype")
            pairs.add((key, lib))
            if key in eligible:
                library_donors[lib].add(key)
    if {lib for _, lib in pairs} != sc_libraries or {key for key, _ in pairs} != clinical.keys():
        raise ValueError("Clinical anchor donors or sc library IDs do not reconcile")
    cases = {key for key, group in eligible.items() if group == "case"}
    mixed = [members for members in library_donors.values() if members & cases and members - cases]
    total = 0
    for choice in combinations(eligible, len(cases)):
        assigned = set(choice)
        if all(len(assigned & members) == len(cases & members) for members in library_donors.values()):
            total += 1
    return {"anchor_sc_links": len(pairs), "anchor_sc_libraries": len(sc_libraries),
            "anchor_sc_donors": len(clinical), "eligible_anchor_libraries": len(library_donors),
            "mixed_anchor_libraries": len(mixed), "unique_donors_in_mixed_libraries": len(set().union(*mixed)),
            "label_allocations_preserving_all_eligible_library_case_totals": total}


def verify_workbook(path, sdrf=None):
    with path.open("rb") as source:
        if file_digest(source, "sha256").hexdigest() != WORKBOOK_SHA256:
            raise ValueError("Unrecognized HECA workbook checksum")
    workbook = load_workbook(path, read_only=True, data_only=True)
    try:
        clinical = list(workbook["Table1_donors"].values)
        libraries = list(workbook["Table2_CellRanger"].values)
        counts = list(workbook["Table3_cell_type_counts"].values)
    finally:
        workbook.close()
    if (list(clinical[2][3:8]) != ["Locationd", "Binary Stagee", "Stagef", "Hormonal treatmentg ",
                                    "Endometriosis_stageh"]
            or list(counts[2][:4]) != ["Dataset_donorID", "1 | SOX9 basalis", "2 | SOX9 functionalis I",
                                      "3 | SOX9 functionalis II"]
            or libraries[2][0] != "Library ID"):
        raise ValueError("HECA metadata headers differ")

    con = sqlite3.connect(":memory:")
    try:
        con.executescript("""
            CREATE TABLE donors (
                id TEXT PRIMARY KEY, study TEXT, location TEXT, phase TEXT, hormones TEXT,
                disease TEXT, ascertainment TEXT, pathology TEXT);
            CREATE TABLE abundance (id TEXT PRIMARY KEY, cells INTEGER, glands INTEGER);
        """)
        studies = Counter()
        for row in clinical[3:]:
            if len(row) < 11 or not row[1]:
                continue
            donor, dataset, _, tissue, stage, _, hormones, diagnosis, source, pathology = row[:10]
            studies[dataset] += 1
            if dataset == "Mareckova-nuclei":
                continue
            if diagnosis not in (0, 1, 2, 3, 4, "TBC"):
                raise ValueError("Unexpected diagnosis convention")
            # Table 1 footnote h documents confirmed disease with unknown stage as TBC.
            group = "control" if diagnosis == 0 else "case"
            study, key = unique_clinical_key(donor, dataset)
            con.execute("INSERT INTO donors VALUES (?,?,?,?,?,?,?,?)",
                        (key, study, tissue, stage, hormones, group, source, pathology))
        for row in counts[3:]:
            if not row or not row[0]:
                continue
            values = row[1:37]
            if len(values) != 36 or any(int(v) != v or v < 0 for v in values):
                raise ValueError("Expected 36 nonnegative integer cell counts")
            con.execute("INSERT INTO abundance VALUES (?,?,?)", (row[0], sum(values), row[2] + row[3]))
        con.executescript("""
            CREATE TEMP VIEW joined AS SELECT * FROM donors INNER JOIN abundance USING (id);
            CREATE TEMP VIEW proliferative AS SELECT * FROM joined WHERE phase='Proliferative'
                AND hormones='NA';
            CREATE TEMP VIEW eligible AS SELECT * FROM proliferative WHERE location='endometrium'
                AND glands>=20;
        """)
        sc = con.execute("SELECT COUNT(*), SUM(cells) FROM joined").fetchone()
        unmatched = con.execute("""SELECT (SELECT COUNT(*) FROM donors) + (SELECT COUNT(*) FROM abundance)
                                   - 2*(SELECT COUNT(*) FROM joined)""").fetchone()[0]
        if unmatched:
            raise ValueError("Unmatched clinical or count keys in scRNA component")
        groups = {study: {"case": 0, "control": 0} for study, in con.execute("SELECT DISTINCT study FROM donors")}
        for study, disease, n in con.execute("SELECT study,disease,COUNT(*) FROM eligible GROUP BY study,disease"):
            groups[study][disease] = n
        mixed = [study for study in groups if min(groups[study].values()) > 0]
        report = {
            "total_donors": sum(studies.values()), "sc_donors": sc[0], "sc_cells": sc[1],
            "sn_donors": studies["Mareckova-nuclei"] + studies["Mareckova-cells & nuclei"],
            "shared_sc_sn_donors": studies["Mareckova-cells & nuclei"],
            "sc_natural_proliferative": con.execute("SELECT COUNT(*) FROM proliferative").fetchone()[0],
            "sc_natural_proliferative_superficial": con.execute(
                "SELECT COUNT(*) FROM proliferative WHERE location='endometrium'").fetchone()[0],
            "eligible_sc_donors": con.execute("SELECT COUNT(*) FROM eligible").fetchone()[0],
            "eligible_gland_cells": con.execute("SELECT SUM(glands) FROM eligible").fetchone()[0],
            "eligible_by_study": dict(sorted(groups.items())),
            "mixed_study_sc_donors": sum(sum(groups[study].values()) for study in mixed),
            "eligible_controls_with_other_pathology": con.execute(
                "SELECT COUNT(*) FROM eligible WHERE disease='control' AND pathology!='C'").fetchone()[0],
            "qc_library_rows": dict(Counter(r[1] for r in libraries[3:] if len(r) > 1 and r[1])),
        }
        if sdrf is not None:
            clinical_anchor = {key: (group, phase) for key, group, phase in con.execute(
                "SELECT id,disease,phase FROM donors WHERE study='Mareckova'")}
            eligible_anchor = dict(con.execute("SELECT id,disease FROM eligible WHERE study='Mareckova'"))
            sc_libraries = {r[0] for r in libraries[3:] if len(r) > 1 and r[1] == "scRNAseq"}
            report["sdrf"] = verify_anchor_libraries(sdrf, sc_libraries, clinical_anchor, eligible_anchor)
        return report
    finally:
        con.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("--sdrf", type=Path, help="Optional original E-MTAB-14039 SDRF")
    args = parser.parse_args()
    print(json.dumps(verify_workbook(args.workbook, args.sdrf), indent=2))
