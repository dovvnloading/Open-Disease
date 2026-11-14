"""Verify the locked Melbourne genes from underlying XLSX XML and exact decimal p-values."""
from __future__ import annotations

import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re
from xml.etree import ElementTree as ET
from zipfile import ZipFile

NS = {
    "s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
WORKBOOK = "b83a3d32a5882f750d978f4de1a6473fdd5d7627dbccf7d7c1774162e104e8e3"
LOCK = "ba140f2bb026b0416e1c705ab28f616f430cc4295f94d174feaf5d15947c1e31"
CANDIDATES = {"ABO": 28, "KIF3C": 3797, "ZNF835": 90485, "GLTP": 51228, "F11": 2160}
DIRECTIONS = ["higher", "lower", "lower", "higher", "lower"]
HEADER = ["PROBE_ID", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "FC", "ILMN_GENE", "CHROMOSOME"]
TITLE = ("Table S20. Differential expression of genes expressed in >90% of samples. "
         "Differentially expressed genes in Endometriosis Cases and Controls, ranked by p-value")


def sheet_rows(workbook):
    """Resolve the S20 relationship and stream values with their explicit sheet column."""
    with ZipFile(workbook) as archive:
        doc = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet = [s for s in doc.findall("s:sheets/s:sheet", NS) if s.get("name") == "TableS20"]
        if len(sheet) != 1:
            raise ValueError("Expected exactly one TableS20")
        ref = sheet[0].get("{%s}id" % NS["r"])
        rel = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        filenames = [r.get("Target") for r in rel if r.get("Id") == ref]
        if len(filenames) != 1 or not re.fullmatch(r"worksheets/sheet[0-9]+\.xml", filenames[0]):
            raise ValueError("Unrecognized TableS20 reference")
        strings = ["".join(si.itertext()) for si in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        with archive.open("xl/" + filenames[0]) as source:
            for _, row in ET.iterparse(source, events=["end"]):
                if row.tag != "{%s}row" % NS["s"]:
                    continue
                index = int(row.get("r"))
                cells = {}
                for cell in row.findall("s:c", NS):
                    address = re.fullmatch(r"([A-Z]+)([1-9][0-9]*)", cell.get("r"))
                    if not address or int(address[2]) != index or address[1] in cells:
                        raise ValueError("Inconsistent XLSX coordinates")
                    text = cell.findtext("s:v", "", NS)
                    cells[address[1]] = strings[int(text)] if text and cell.get("t") == "s" else text
                yield index, cells
                row.clear()


def verify(workbook, lock):
    if hashlib.sha256(workbook.read_bytes()).hexdigest() != WORKBOOK or hashlib.sha256(lock.read_bytes()).hexdigest() != LOCK:
        raise ValueError("Workbook or candidate lock differs from frozen inputs")
    nominees = [(int(x["id"]), x["direction"]) for x in json.loads(lock.read_text())["candidates"] if x["family"] == "gene"]
    if nominees != list(zip(CANDIDATES.values(), DIRECTIONS, strict=True)):
        raise ValueError("Unexpected candidate family")
    records = sheet_rows(workbook)
    first, second = next(records), next(records)
    if first[0] != 1 or first[1]["A"] != TITLE or second[0] != 2 or [second[1][k] for k in "ABCDEFGHIJ"] != HEADER:
        raise ValueError("Unexpected source heading")
    probes = set()
    matched = {gene: [] for gene in CANDIDATES}
    for expected, (actual, row) in enumerate(records, start=3):
        probe = row["A"]
        if actual != expected or probe in probes or not re.fullmatch(r"ILMN_[0-9]+", probe):
            raise ValueError("Nonunique or unrecognized source probe row")
        probes.add(probe)
        gene = row.get("I", "").strip()
        if gene in matched:
            coef, p = Decimal(row["B"]), Decimal(row["E"])
            if not coef.is_finite() or not p.is_finite() or not 0 <= p <= 1:
                raise ValueError("Invalid published coefficient or significance")
            matched[gene].append((probe, coef, p))
    if len(probes) != 15_262:
        raise ValueError("Unexpected source panel size")

    one = Decimal(1)
    raw = {gene: min(one, len(rows) * min(r[2] for r in rows)) if rows else one for gene, rows in matched.items()}
    # Holm compares the minimum of k remaining p-values with its k-fold threshold.
    ordered = sorted(raw, key=raw.get)
    adjusted = {
        g: min(one, max((len(ordered) - k) * raw[ordered[k]] for k in range(i + 1)))
        for i, g in enumerate(ordered)
    }
    return {
        "source_probes": len(probes), "gene_family_size": len(raw),
        "genes": [{
            "gene": g, "entrez": CANDIDATES[g], "tested": bool(matched[g]), "matched_probes": len(matched[g]),
            "probes": [{"id": probe, "logFC_reported": float(coef), "nominal_p": float(p)} for probe, coef, p in matched[g]],
            "p_within_gene": float(raw[g]) if matched[g] else None, "holm_p_for_bookkeeping": float(adjusted[g]),
        } for g in CANDIDATES],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("nomination_lock", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.workbook, args.nomination_lock), indent=2, allow_nan=False))
