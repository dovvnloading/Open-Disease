"""Verify the ovarian RNA-seq panel independently from OOXML cells and exact decimal sums."""
import argparse
import json
from collections import defaultdict
from decimal import Decimal
from itertools import product
from pathlib import Path
from xml.etree import ElementTree as ET
from zipfile import ZipFile

from open_disease.geo_microarray import sha256_file


NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"


def verify(book, receipt):
    if sha256_file(book) != receipt["input_sha256"][book.name]:
        raise ValueError("Input workbook changed")
    totals = defaultdict(lambda: [Decimal(0)] * 10)
    with ZipFile(book) as archive:
        strings = ["".join(n.text or "" for n in item.iter(NS + "t")) for item in
                   ET.parse(archive.open("xl/sharedStrings.xml")).getroot()]
        header = None
        with archive.open("xl/worksheets/sheet1.xml") as stream:
            for _, item in ET.iterparse(stream, events=("end",)):
                if item.tag != NS + "row":
                    continue
                cells = [c.findtext(NS + "v") for c in item]
                if header is None:
                    header = [strings[int(i)] for i in cells[4:]]
                else:
                    vector = totals[strings[int(cells[1])]]
                    for i, text in enumerate(cells[4:]):
                        vector[i] += Decimal(text)
                item.clear()
    anchors = receipt["panel_genes"]
    scores = [0.0] * 10
    observed = {}
    for gene, r in anchors.items():
        target = totals[r["ensembl"]]
        pts = []
        for i, v in enumerate(target):
            below = sum(a[i] < v for a in totals.values())
            ties = sum(a[i] == v for a in totals.values())
            pts.append(100 * (below + ties / 2) / len(totals))
            scores[i] += r["direction"] * pts[-1] / len(anchors)
        observed[gene] = pts
    # Pair by worksheet aliases, independently of native GEO row ordering.
    pairs = [(header.index("OSA_" + key[3:]), i) for i, key in enumerate(header) if key.startswith("EM_")]
    diffs = [scores[l] - scores[u] for l, u in pairs]
    effect = sum(diffs)/len(diffs)
    if abs(effect - receipt["primary_pp"]["effect_pp"]) > 1e-8:
        raise ValueError("Decimal verification disagrees with panel percentile shift")
    for name, pts in observed.items():
        delta = sum(pts[l] - pts[u] for l, u in pairs) / len(pairs)
        if abs(delta - anchors[name]["difference_pp"]) > 1e-8:
            raise ValueError("Decimal verification disagrees with a gene percentile shift")
    upper = sum(sum(d*s for d, s in zip(diffs, signs)) >= sum(diffs) - 1e-9
                for signs in product((-1, 1), repeat=len(diffs)))
    if upper != receipt["primary_pp"]["upper"]:
        raise ValueError("Decimal verification disagrees with exact test")
    return {"verified_genes": len(anchors), "reference_genes": len(totals), "pairs": len(pairs),
            "panel_effect_pp": effect, "upper_tail_assignments": upper}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("book", type=Path)
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.book, json.loads(args.receipt.read_text())), indent=2))
