"""Audit HECA Supplement tables 1–3 from a local workbook, reporting aggregate counts only."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import re

from openpyxl import load_workbook
from open_disease.geo_microarray import sha256_file

SIZE = 400034
SHA256 = "f70030404af3f6070c10f22a064e20188b0fb3251c1042605563c0aea00d2ed3"
SOURCE = "https://pmc-oa-opendata.s3.amazonaws.com/PMC11387200.1/41588_2024_1873_MOESM3_ESM.xlsx"
COHORT = {"Mareckova-cells": "Mareckova", "Mareckova-cells & nuclei": "Mareckova",
          "Huangb": "Huang", "Lai": "Lai", "Tan": "Tan", "Wang": "Wang",
          "Garcia-Alonso": "GarciaAlonso", "Fonseca": "Fonseca"}
STATE_I = "2 | SOX9 functionalis I"
STATE_II = "3 | SOX9 functionalis II"
STAGES = ("all_sc_donors", "proliferative", "untreated_superficial_proliferative", "at_least_20_gland_cells")


def check(ok, message):
    if not ok:
        raise ValueError(message)


def records(sheet, identifier, required):
    """The printed headers follow the title and blank row. Table 1 has an unnamed key."""
    raw = iter(sheet.values)
    next(raw)
    next(raw)
    names = list(next(raw))
    if sheet.title == "Table1_donors":
        check(names[:2] == [None, "Dataset"], "Unexpected clinical key and adjacent dataset header")
        names[0] = identifier
    names = [(v.strip() if isinstance(v, str) else v) for v in names]
    present = [v for v in names if v is not None]
    check(len(set(present)) == len(present) and set(required) <= set(present), "Missing or duplicate table field")
    rows, notes = [], []
    for values in raw:
        row = {name: value for name, value in zip(names, values) if name is not None}
        if not any(v is not None for v in row.values()):
            continue
        if sheet.title == "Table1_donors" and row["Dataset"] is None:
            check(isinstance(row[identifier], str) and sum(v is not None for v in row.values()) == 1,
                  "Malformed clinical footnote")
            notes.append(row[identifier])
        else:
            check(isinstance(row[identifier], str) and bool(row[identifier]), "Missing source key")
            rows.append(row)
    check(len({r[identifier] for r in rows}) == len(rows), "Duplicate source key")
    return rows, notes


def is_case(row):
    stage = row["Endometriosis_stageh"]
    check(stage in (0, 1, 2, 3, 4, "TBC"), "Unknown disease-stage code")
    pathology = {x.strip() for x in row["Endometrial pathologyi"].split("+")}
    check((stage != 0) == ("E" in pathology), "Clinical stage and pathology disagree")
    return stage != 0


def sc_count_key(study, donor, footnoted_lai):
    # Footnote a identifies exactly one Lai source label with a phase reassignment.
    if study == "Lai" and donor == footnoted_lai + "a":
        donor = footnoted_lai
    return donor if study == "Fonseca" else f"{study}_{donor}"


def count_pair(counter):
    return {"cases": counter[True], "controls": counter[False]}


def summarize(clinical, cell_counts, qc, notes):
    footnoted = [m[1] for text in notes if (m := re.match(r"^a for donor (\S+) the original publication", text))]
    check(len(footnoted) == 1, "Expected one explicit Lai phase footnote")
    state_columns = [k for k in cell_counts[0] if re.match(r"^\d+ \| ", k)]
    check(len(state_columns) == 36 and {STATE_I, STATE_II} <= set(state_columns), "Unexpected cell-state universe")
    indexed = {row["Dataset_donorID"]: row for row in cell_counts}
    total_cells = 0
    for row in cell_counts:
        for name in state_columns:
            n = row[name]
            check(isinstance(n, (int, float)) and 0 <= n < 1e7 and int(n) == n, "Invalid cell count")
            total_cells += int(n)
    check(len(indexed) == 63 and total_cells == 313527, "Source scRNA denominator mismatch")
    cohorts = defaultdict(lambda: {stage: Counter() for stage in STAGES})
    cell_yield, verified, other_pathologies = Counter(), defaultdict(Counter), Counter()
    used = set()
    clinical_status = Counter()
    for row in clinical:
        case = is_case(row)
        clinical_status[case] += 1
        dataset = row["Dataset"]
        if dataset == "Mareckova-nuclei":
            continue
        check(dataset in COHORT, "Unrecognized contributing dataset")
        study = COHORT[dataset]
        label = sc_count_key(study, row["source_key"], footnoted[0])
        check(label in indexed and label not in used, "Clinical-to-abundance join missing or ambiguous")
        used.add(label)
        cohort = cohorts[study]
        cohort[STAGES[0]][case] += 1
        if row["Binary Stagee"] != "Proliferative":
            continue
        cohort[STAGES[1]][case] += 1
        if row["Locationd"] != "endometrium" or row["Hormonal treatmentg"] != "NA":
            continue
        cohort[STAGES[2]][case] += 1
        n = int(indexed[label][STATE_I] + indexed[label][STATE_II])
        if n < 20:
            continue
        cohort[STAGES[3]][case] += 1
        cell_yield[study] += n
        verified[study][row["Endometriosis_stage assignment"]] += 1
        other_pathologies[study] += row["Endometrial pathologyi"] not in {"C", "E"}
    check(used == set(indexed), "Incomplete scRNA donor coverage")
    mixed = [s for s in sorted(cohorts) if all(cohorts[s][STAGES[-1]][case] for case in (True, False))]
    qc_counts = Counter(r["snRNAseq"] for r in qc)
    check(qc_counts == {"scRNAseq": 22, "snRNAseq": 41, "spatial": 8}, "Unexpected QC modality census")
    check(set(qc[0]) == {"Library ID", "snRNAseq", "Number of Cells", "Reads per Cell", "Genes per Cell",
                         "Fraction Reads in Spots Under Tissue", "Mean Reads per Spot", "Mean Reads Under Tissue per Spot",
                         "Median UMI Counts per Spot", "Median Genes per Spot", "Genes Detected"}, "QC schema changed")
    check(not ({r["Library ID"] for r in qc} & {r["source_key"] for r in clinical}),
          "Unexpected shared key between clinical and QC tables")
    by_study = {s: {t: count_pair(counts[t]) for t in STAGES} for s, counts in sorted(cohorts.items())}
    return {"source": SOURCE, "bytes": SIZE, "sha256": SHA256,
            "clinical_donors_across_assays": count_pair(clinical_status),
            "donors_in_both_assays": sum(r["Dataset"] == "Mareckova-cells & nuclei" for r in clinical),
            "sc_cells_all_states": total_cells,
            "sc_attrition": {t: count_pair(sum((c[t] for c in cohorts.values()), Counter())) for t in STAGES},
            "sc_attrition_by_study": by_study, "mixed_studies": mixed,
            "donors_in_mixed_studies": count_pair(sum((cohorts[s][STAGES[-1]] for s in mixed), Counter())),
            "selected_gland_cells_in_mixed_studies": sum(cell_yield[s] for s in mixed),
            "diagnosis_ascertainment_in_mixed_studies": dict(sum((verified[s] for s in mixed), Counter())),
            "additional_endometrial_pathologies_in_mixed_studies": sum(other_pathologies[s] for s in mixed),
            "library_qc_records_by_modality": dict(sorted(qc_counts.items())),
            "supplement_contains_donor_to_capture_crosswalk": False}


def audit(workbook):
    check(workbook.stat().st_size == SIZE and sha256_file(workbook) == SHA256, "Unrecognized HECA workbook")
    book = load_workbook(workbook, read_only=True, data_only=True)
    try:
        clinical, notes = records(book["Table1_donors"], "source_key", ("Dataset", "Locationd", "Binary Stagee",
                                  "Hormonal treatmentg", "Endometriosis_stageh", "Endometriosis_stage assignment",
                                  "Endometrial pathologyi"))
        counts, _ = records(book["Table3_cell_type_counts"], "Dataset_donorID", ("Dataset_donorID", STATE_I, STATE_II))
        qc, _ = records(book["Table2_CellRanger"], "Library ID", ("Library ID", "snRNAseq"))
        check(len(clinical) == 121, "Unexpected clinical donor count")
        return summarize(clinical, counts, qc, notes)
    finally:
        book.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    print(json.dumps(audit(parser.parse_args().workbook), indent=2, allow_nan=False))
