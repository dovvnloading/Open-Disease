"""Render the aggregate phase-specific physical scan overlap without extra dependencies."""
import json
from pathlib import Path
from xml.etree import ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
DIR = ROOT / "evidence/matrisome-holdout-2026-11-15"


def draw(overlap):
    svg = ET.Element("svg", xmlns="http://www.w3.org/2000/svg", viewBox="0 0 760 260", role="img")
    ET.SubElement(svg, "title").text = "Training scan reuse in the 37-array endometriosis test set"
    ET.SubElement(svg, "desc").text = "Shared training scans by phase: proliferative 9 of 11, early secretory 8 of 9, mid-secretory 17 of 17."
    ET.SubElement(svg, "rect", width="760", height="260", fill="#fff")
    def label(x, y, value, size=16, color="#183449", weight="normal"):
        el = ET.SubElement(svg, "text", x=str(x), y=str(y), fill=color,
                           attrib={"font-family": "system-ui, sans-serif", "font-size": str(size), "font-weight": weight})
        el.text = value
    label(25, 37, f"{overlap['matched_test_scans']} of {overlap['test_records']} test scans also appeared in training", 25, weight="bold")
    colors = ("#11666e", "#d3dce2")
    for i, (phase, display) in enumerate((("proliferative", "Proliferative"), ("early_secretory", "Early secretory"), ("mid_secretory", "Mid-secretory"))):
        counts = overlap["phase"][phase]
        matched = counts["matched"]
        total = matched + counts.get("unmatched", 0)
        y = 68 + 47*i
        label(25, y+19, display)
        for j in range(total):
            ET.SubElement(svg, "rect", x=str(186 + j*25), y=str(y), width="21", height="24",
                          rx="2", fill=colors[int(j >= matched)])
        label(634, y+19, f"{matched} / {total}", weight="bold")
    for x, color, display in ((25, colors[0], "Scan found in training"), (300, colors[1], "No training scan match found")):
        ET.SubElement(svg, "rect", x=str(x), y="222", width="17", height="17", rx="2", fill=color)
        label(x+26, 236, display, 15)
    ET.indent(svg)
    return ET.tostring(svg, encoding="unicode") + "\n"


if __name__ == "__main__":
    (DIR / "scan-overlap.svg").write_text(draw(json.loads((DIR / "overlap.json").read_text())))
