"""Plot the two independent case–control estimates for every locked nominee.

Usage: python scripts/plot_cross_cohort_validation.py
Requires matplotlib. The defaults use only aggregate files tracked in Git.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


EVIDENCE = Path(__file__).resolve().parents[1] / "evidence/cross-cohort-2026-11-14"
# Display labels checked against the pinned GPL571 Entrez/symbol annotation.
GENE_LABELS = {"28": "ABO", "3797": "KIF3C", "90485": "ZNF835", "51228": "GLTP", "2160": "F11"}
FINLAND, ROME = "#135E80", "#B44E23"
INK, SECONDARY = "#172B35", "#40535E"
ALT = ("Case-minus-control relative-expression ranks in Finland and Rome for five genes "
       "and five pathways. ABO is higher in both cohorts. The other nine estimates "
       "reverse direction, including all five pathways. Each cohort uses its own array reference.")


def read_effects(lock_path, results_path):
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    expected = [(c["family"], c["id"]) for c in lock["candidates"]]
    with results_path.open(encoding="utf-8", newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    by_id = {(r["family"], r["id"]): r for r in rows}
    if set(expected) != by_id.keys() or len(rows) != len(set(expected)) or len(expected) != 10:
        raise ValueError("Expected all ten unique locked candidates and exactly one comparison per candidate")
    output = {"gene": [], "pathway": []}
    for c in lock["candidates"]:
        family, key = c["family"], c["id"]
        row = by_id[family, key]
        fin, rome = float(row["finnish_effect_pp"]), float(row["rome_effect_pp"])
        if (not math.isfinite(fin) or not math.isfinite(rome) or fin != c["effect_pp"]
                or row["direction"] != c["direction"]):
            raise ValueError("Reported effect or direction disagrees with the lock")
        label = GENE_LABELS[key] if family == "gene" else c["name"]
        output[family].append((label, fin, rome))
    if any(len(values) != 5 for values in output.values()):
        raise ValueError("Expected five gene nominees and five pathway nominees")
    return output


def draw(groups, output):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "text.color": INK,
                         "axes.labelcolor": INK, "xtick.color": SECONDARY, "ytick.color": INK,
                         "svg.fonttype": "none", "svg.hashsalt": "open-disease-2026-11-14"})
    fig, axes = plt.subplots(2, 1, figsize=(8.6, 6.45))
    fig.subplots_adjust(left=0.355, right=0.966, top=0.756, bottom=0.19, hspace=0.59)
    fig.text(0.039, 0.954, "Endometrial expression across cohorts",
             fontsize=17.5, fontweight="bold", va="top")
    fig.text(0.039, 0.896, "Case minus control, in within-array percentile points", fontsize=11.2)
    handles = [Line2D([], [], linestyle="none", marker=marker, markersize=7, color=color)
               for marker, color in (("o", FINLAND), ("D", ROME))]
    fig.legend(handles, ["Finland · 12 cases, 13 controls", "Rome · 8 cases, 6 controls"],
               loc="center left", bbox_to_anchor=(0.027, 0.846), ncol=2, frameon=False,
               handletextpad=0.2, columnspacing=1.9, fontsize=9.6)
    for ax, (family, section, step) in zip(axes, [("gene", "Genes", 10), ("pathway", "Pathways", 2)]):
        values = groups[family]
        labels, fin, rome = zip(*values)
        limit = step * math.ceil(1.07 * max(abs(x) for x in (*fin, *rome)) / step)
        ax.set_xlim(-limit, limit)
        ax.set_xticks(range(-limit, limit + 1, step))
        ax.set_ylim(4.52, -0.52)
        ax.set_yticks(range(5), labels)
        for y, a, b in zip(range(5), fin, rome):
            ax.plot([a, b], [y, y], color="#ABB9C0", linewidth=1.25, zorder=2, solid_capstyle="round")
        ax.scatter(fin, range(5), marker="o", s=55, facecolor=FINLAND, edgecolor="white",
                   linewidth=0.8, zorder=3)
        ax.scatter(rome, range(5), marker="D", s=44, facecolor=ROME, edgecolor="white",
                   linewidth=0.8, zorder=4)
        ax.axvline(0, color="#536672", linewidth=1.05, zorder=1)
        ax.set_axisbelow(True)
        ax.xaxis.grid(True, linewidth=0.6, color="#D6DFE3", linestyle=(0, (2, 3)))
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.spines["bottom"].set_color("#B4C0C7")
        ax.tick_params(axis="y", length=0, pad=7, labelsize=9.2)
        ax.tick_params(axis="x", length=3, labelsize=8.8, pad=4)
        fig.text(0.039, ax.get_position().y1 + 0.013, section, fontsize=11.5, fontweight="bold")
    fig.text(0.039, 0.092, "Positive values indicate higher relative expression in cases.",
             fontsize=9.3, color=SECONDARY)
    fig.text(0.039, 0.057,
             "Panel scales differ. Percentiles use each platform’s own probe set.",
             fontsize=9.3, color=SECONDARY)
    fig.savefig(output, format="svg", metadata={"Title": "Endometrial expression across cohorts",
                                                "Description": ALT, "Date": None})
    plt.close(fig)
    # Matplotlib adds trailing spaces in multiline SVG path attributes.
    output.write_text("\n".join(line.rstrip() for line in output.read_text().splitlines()) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=EVIDENCE / "nominations.json")
    parser.add_argument("--results", type=Path, default=EVIDENCE / "validation.tsv")
    parser.add_argument("--output", type=Path, default=EVIDENCE / "validation.svg")
    args = parser.parse_args()
    draw(read_effects(args.lock, args.results), args.output)
