"""Draw the matched-membership comparison using only published aggregate JSON.

From the repository root, run python scripts/plot_swiss_glandular.py.
The optional plotting dependency is matplotlib.
"""
from __future__ import annotations

import argparse
import json
from math import isclose, isfinite
from pathlib import Path
from statistics import fmean

import matplotlib

matplotlib.use("Agg")
from matplotlib import pyplot as plt

EVIDENCE = Path(__file__).resolve().parents[1] / "evidence/swiss-glandular-2026-11-14"
PATHWAYS = {
    "R-HSA-9008059": "IL-37 signaling",
    "R-HSA-5654741": "Signaling by FGFR3",
    "R-HSA-163560": "Triglyceride catabolism",
    "R-HSA-5654710": "PI-3K cascade downstream\nof FGFR3",
    "R-HSA-432722": "Golgi-associated vesicle\nbiogenesis",
}
COLUMNS = (
    ("Finland", "Whole-tissue array", "#155e80", "o"),
    ("Rome", "Whole-tissue array", "#a54820", "D"),
    ("Switzerland", "Donor-first ranks", "#186e61", "s"),
    ("Switzerland", "Cell-first ranks", "#186e61", "s"),
)
INK, SECONDARY = "#182d36", "#475b64"
TITLE = "Endometrial pathway scores across three cohorts"


def read_comparison(prior_file, swiss_file):
    prior = json.loads(prior_file.read_text(encoding="utf-8"))
    swiss = json.loads(swiss_file.read_text(encoding="utf-8"))
    endpoints = swiss["endpoints"]
    if prior.keys() != PATHWAYS.keys() or endpoints.keys() != {"GLOBAL", *PATHWAYS}:
        raise ValueError("Aggregate inputs must contain all five fixed pathways and the Swiss joint test")
    rows = []
    for key in PATHWAYS:
        values = (prior[key]["finland_matched_pp"], prior[key]["rome_matched_pp"],
                  endpoints[key]["difference"], endpoints[key]["cell_first"])
        if not all(isfinite(x) and abs(x) <= 100 for x in values):
            raise ValueError("Expected finite percentile-point differences")
        rows.append(values)
    joint = tuple(fmean(column) for column in zip(*rows, strict=True))
    for value, field in zip(joint[-2:], ("difference", "cell_first"), strict=True):
        if not isclose(value, endpoints["GLOBAL"][field], rel_tol=0, abs_tol=1e-10):
            raise ValueError("Swiss joint estimate differs from the equal-weight five-pathway mean")
    return [joint, *rows], endpoints["GLOBAL"]["extreme"], swiss["assignments"]


def draw(rows, extreme, assignments, output):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "svg.fonttype": "none", "svg.hashsalt": "swiss-glandular-2026-11-14",
                         "text.color": INK, "axes.labelcolor": INK, "xtick.color": SECONDARY})
    fig = plt.figure(figsize=(12.0, 6.55), facecolor="white")
    fig.text(.032, .949, TITLE, size=20, weight="bold", va="top")
    fig.text(.032, .868,
             "Case-minus-control difference in relative expression, in within-cohort percentile points.",
             size=12.0)
    left, right, bottom, top = .315, .980, .238, .720
    gap = .017
    width = (right - left - 3 * gap) / 4
    labels = ["Five-pathway mean", *PATHWAYS.values()]
    y_min, y_max = -.56, len(labels) - .41
    # Map one figure-level label to each row across the four aligned, equally scaled panels.
    for i, label in enumerate(labels):
        fraction = (i - y_min) / (y_max - y_min)
        fig.text(.032, top - fraction * (top - bottom), label, ha="left", va="center",
                 fontsize=12.5 if i == 0 else 11.3, fontweight="bold" if i == 0 else "normal",
                 linespacing=1.06)
    fig.text(.032, .756, "Endpoint", size=10.8, color=SECONDARY, weight="bold")
    for j, (cohort, method, color, marker) in enumerate(COLUMNS):
        x0 = left + j * (width + gap)
        ax = fig.add_axes([x0, bottom, width, top - bottom])
        ax.set_xlim(-5.65, 5.65)
        ax.set_ylim(y_max, y_min)
        ax.set_yticks([])
        ax.set_xticks([-5, 0, 5], ["−5", "0", "+5"])
        ax.tick_params(axis="x", length=0, pad=7, labelsize=9.2)
        ax.spines[:].set_visible(False)
        ax.axhspan(y_min, .5, color="#edf2f3", zorder=0)
        for i in range(len(labels) - 1):
            ax.axhline(i + .5, lw=.55, color="#d7e0e2", zorder=0)
        for tick in (-5, 5):
            ax.axvline(tick, color="#d6dee0", lw=.55, linestyle=(0, (2, 3)), zorder=1)
        ax.axvline(0, color="#86989e", lw=.9, zorder=1)
        fig.text(x0 + width / 2, .801, cohort, ha="center", color=color, weight="bold", size=13)
        fig.text(x0 + width / 2, .759, method, ha="center", size=10.2, color=SECONDARY)
        for i, row in enumerate(rows):
            value = row[j]
            ax.plot([0, value], [i + .095, i + .095], color=color, lw=1.35, zorder=2)
            ax.plot(value, i + .095, marker=marker, markerfacecolor="white" if j == 3 else color,
                    markeredgecolor=color, markersize=5.8, markeredgewidth=1.2, zorder=3)
            label = f"{value:+.2f}".replace("-", "−")
            ax.annotate(label, (value, i + .095), xytext=(0, 6), textcoords="offset points",
                        va="bottom", ha="center", size=9.5, color=INK,
                        weight="bold" if i == 0 else "normal")

    fig.text(.032, .148,
             "Each pathway uses the same measured genes across cohorts. Positive values indicate higher scores in cases.",
             size=10.4)
    fig.text(.032, .103,
             "Finland and Rome rank features on their respective arrays. Switzerland ranks 29,034 glandular epithelial features.",
             size=10.4, color=SECONDARY)
    fig.text(.032, .058,
             f"The Swiss donor-first joint test gave an exact two-sided p of {extreme / assignments:.3f} "
             f"({extreme} of {assignments} allocations within capture batches).",
             size=10.4, color=SECONDARY)
    description = ("Case minus control in within-cohort percentile points. The five-pathway mean is "
                   + ", ".join(f"{cohort} ({method}) {value:+.2f}"
                               for (cohort, method, _, _), value in zip(COLUMNS, rows[0], strict=True))
                   + ". Each of the five individual pathways is also shown. The ranking universe differs by assay.")
    fig.savefig(output, format="svg", metadata={"Title": TITLE, "Description": description, "Date": None})
    plt.close(fig)
    svg = output.read_text(encoding="utf-8")
    output.write_text("\n".join(line.rstrip() for line in svg.splitlines()) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prior", type=Path, default=EVIDENCE / "prior-cohorts.json")
    parser.add_argument("--swiss", type=Path, default=EVIDENCE / "results.json")
    parser.add_argument("--output", type=Path, default=EVIDENCE / "comparison.svg")
    options = parser.parse_args()
    draw(*read_comparison(options.prior, options.swiss), options.output)
