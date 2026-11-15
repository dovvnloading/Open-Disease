"""Cross-check reference parsing and all seven exact tests on the audited cohort loaders.

The verifier reads the workbook directly with OOXML, independently maps symbols,
computes assignment standard deviations over the full enumerated population and
uses chunked matrix multiplication for every permuted gene contrast.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import gzip
from itertools import combinations, product
import json
from pathlib import Path
import re
from statistics import mean
from xml.etree import ElementTree as ET
import zipfile

import numpy as np

from analyze_melbourne_genes import SOURCE_HASH
from genome_endometdb import load_finland
from genome_rome import load_rome
from open_disease.geo_microarray import sha256_file


def direct_reference(workbook, platform, finnish_ids):
    by_name = defaultdict(set)
    by_id = defaultdict(set)
    with gzip.open(platform, "rt", newline="") as stream:
        for line in stream:
            if line.strip() == "!platform_table_begin":
                break
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["ID"] == "!platform_table_end":
                break
            gene, symbol = row["Gene ID"], row["Gene symbol"]
            if (row["ID"].startswith("AFFX") or not re.fullmatch("[1-9][0-9]*", gene)
                    or not re.fullmatch("[A-Za-z0-9][A-Za-z0-9_.-]*", symbol)):
                continue
            by_name[symbol].add(int(gene))
            by_id[int(gene)].add(symbol)
    bridge = {s: g for g, symbols in by_id.items() for s in symbols
              if len(symbols) == 1 and len(by_name[s]) == 1 and g in finnish_ids}
    if sha256_file(workbook) != SOURCE_HASH:
        raise AssertionError("Workbook hash mismatch")
    out = {}
    with zipfile.ZipFile(workbook) as archive:
        ns = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
        rel = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
        strings = ["".join(t.text or "" for t in si.iter(ns + "t"))
                   for si in ET.fromstring(archive.read("xl/sharedStrings.xml"))]
        relationships = {x.attrib["Id"]: x.attrib["Target"] for x in
                         ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))}
        for sheet in ET.fromstring(archive.read("xl/workbook.xml")).iter(ns + "sheet"):
            label = sheet.attrib["name"]
            if label not in {"TableS3", "TableS20"}:
                continue
            parts = defaultdict(list)
            target = relationships[sheet.attrib[rel]]
            target = target.lstrip("/") if target.startswith("/") else "xl/" + target
            for _, row in ET.iterparse(archive.open(target)):
                if row.tag != ns + "row":
                    continue
                if int(row.attrib["r"]) > 2:
                    cells = {}
                    for cell in row.iter(ns + "c"):
                        col = re.match(r"[A-Z]+", cell.attrib["r"])[0]
                        val = cell.findtext(ns + "v")
                        if val is not None:
                            cells[col] = strings[int(val)] if cell.get("t") == "s" else val
                    symbol = cells.get("I", "").strip()
                    if symbol in bridge:
                        parts[bridge[symbol]].append(float(cells["D"]))
                row.clear()
            ids = sorted(parts)
            out[label] = np.array(ids), np.array([mean(parts[g]) for g in ids])
    return out


def direct_test(scores, cases, strata, reference, block=128):
    options = []
    observed = []
    for s in sorted(set(strata)):
        ids = np.where(strata == s)[0]
        k = sum(cases[ids])
        local = []
        for choice in combinations(range(len(ids)), k):
            weights = np.zeros(len(cases))
            picked = np.zeros(len(ids), dtype=bool)
            picked[list(choice)] = True
            weights[ids[picked]] = 1 / k
            weights[ids[~picked]] = -1 / (len(ids) - k)
            local.append(weights)
            if np.array_equal(picked, cases[ids]):
                observed.append(weights)
        options.append(local)
    obs_weights = np.mean(observed, axis=0)
    all_weights = np.array([np.mean(combo, axis=0) for combo in product(*options)])
    # Each block is a set of genes. Accumulate the correlation's sufficient statistics
    # without holding the assignments-by-transcriptome array in memory.
    n = len(reference)
    y = reference - mean(reference)
    sum_x = np.zeros(len(all_weights))
    sum_x2 = np.zeros(len(all_weights))
    product_xy = np.zeros(len(all_weights))
    obs_x = np.empty(n)
    for start in range(0, n, block):
        sl = slice(start, start + block)
        m = scores[:, sl] - scores[0, sl]
        x = all_weights @ m
        sd = x.std(axis=0)
        x /= sd
        sum_x += x.sum(axis=1)
        sum_x2 += np.einsum("ij,ij->i", x, x)
        product_xy += x @ y[sl]
        obs_x[sl] = (obs_weights @ m) / sd
    rho = product_xy / np.sqrt((sum_x2 - sum_x**2 / n) * (y @ y))
    actual = np.corrcoef(obs_x, reference)[0, 1]
    return float(actual), int(np.count_nonzero(np.abs(rho) >= abs(actual) - 1e-12)), len(all_weights)


def verify(args, native):
    with args.finnish.open() as stream:
        allowed = {int(row["id"]) for row in csv.DictReader(stream, delimiter="\t")}
    refs = direct_reference(args.workbook, args.rome / "GPL571.annot.gz", allowed)
    rome = load_rome(args.rome)
    fin = json.loads(args.finnish_results.read_text())
    rome_expected = {r["sheet"]: r for r in json.loads(args.rome_results.read_text())}
    checks = []

    def check(cohort, which, sel, case, strata, row, field, tail):
        genes, scores, _, _ = cohort
        ids, ref = refs[which]
        col = np.searchsorted(genes, ids)
        np.testing.assert_array_equal(genes[col], ids)
        r, k, n = direct_test(scores[sel][:, col], case, strata, ref)
        if (len(ids), k, n) != (row["genes"], row[tail], row["assignments"]) or abs(r - row[field]) > 1e-10:
            raise AssertionError("Independent reconstruction differs")
        checks.append(dict(genes=len(ids), r=r, extreme=k, assignments=n))

    fg, fs, fc, strata = native
    for which in ("TableS20", "TableS3"):
        check(rome, which, slice(None), rome[2], rome[3], rome_expected[which], "rome_r", "rome_extreme")
        check(native, which, slice(None), fc, strata, fin["disease_contrast"][which], "r", "extreme")
    phases = np.array([s.split("-")[-1] for s in strata])
    for phase in ("proliferative", "secretory"):
        use = phases == phase
        check(native, "TableS3", use, fc[use], strata[use], fin["disease_contrast_by_phase"][phase], "r", "extreme")
    st = np.array([s.split("-")[0] + str(c) for s, c in zip(strata, fc, strict=True)])
    keep = np.isin(st, [s for s in np.unique(st) if len(set(phases[st == s])) == 2])
    check(native, "TableS3", keep, phases[keep] == "secretory", st[keep],
          fin["positive_control_secretory_vs_proliferative"], "r", "extreme")
    return checks


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for flag in ("rome", "finland", "workbook", "finnish", "rome-results", "finnish-results"):
        parser.add_argument("--" + flag, required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args, load_finland(args.finland)), indent=2))
