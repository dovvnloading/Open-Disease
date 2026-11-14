"""Plot the five published cohorts on their common frozen pathway membership."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np

from plot_swiss_glandular import PATHWAYS

ROOT = Path(__file__).resolve().parents[1] / "evidence"
SWISS = ROOT / "swiss-glandular-2026-11-14"
EXTERNAL = ROOT / "external-bulk-2026-11-14"


def matched_effects():
    def read(path):
        return json.loads(path.read_text(encoding="utf-8"))

    prior = read(SWISS / "prior-cohorts.json")
    swiss = read(SWISS / "results.json")["endpoints"]
    coverage = read(SWISS / "mapping.json")["coverage"]
    new = read(EXTERNAL / "results.json")
    values = []
    for key in PATHWAYS:
        # The new restricted set is a subset of the Swiss set by construction.
        # Equal cardinality therefore establishes identity with the earlier set.
        if new["membership"][key]["swiss"] != coverage[key]["mapped"]:
            raise ValueError("Five-cohort memberships no longer coincide")
        values.append([prior[key]["finland_matched_pp"], prior[key]["rome_matched_pp"],
                       swiss[key]["difference"], new["GSE6364"][key]["swiss_members"],
                       new["GSE153739"][key]["swiss_members"]])
    mean = np.mean(values, axis=0)
    stored = [swiss["GLOBAL"]["difference"], new["GSE6364"]["GLOBAL"]["swiss_members"],
              new["GSE153739"]["GLOBAL"]["swiss_members"]]
    if not np.allclose(mean[2:], stored, rtol=0, atol=1e-10):
        raise ValueError("Global estimates disagree with frozen equal weights")
    return np.vstack([mean, values])


def draw(output):
    effects = matched_effects()
    plt.rcParams.update({"font.family": "DejaVu Sans", "svg.fonttype": "none", "svg.hashsalt": "open-disease-bulk-2026"})
    cmap = LinearSegmentedColormap.from_list("relative", ["#9169a8", "#d9c7e3", "#f6f6f3", "#a8d5ca", "#448978"])
    fig, ax = plt.subplots(figsize=(11.4, 5.5))
    fig.subplots_adjust(top=.71, bottom=.24, left=.31, right=.97)
    ax.imshow(effects, cmap=cmap, norm=TwoSlopeNorm(vcenter=0, vmin=-4.7, vmax=4.7), aspect="auto")
    ax.set_xticks(range(5), ["Finland\nwhole tissue", "Rome\nwhole tissue", "Switzerland\nglandular cells",
                             "United States\nwhole tissue¹", "Mumbai\nwhole tissue"])
    ax.xaxis.tick_top()
    ax.set_yticks(range(6), ["Five-pathway mean", *PATHWAYS.values()])
    ax.tick_params(axis="both", which="both", length=0, pad=8, labelsize=10)
    ax.set_xticks(np.arange(-.5, 5), minor=True)
    ax.set_yticks(np.arange(-.5, 6), minor=True)
    ax.grid(which="minor", color="white", linewidth=2)
    ax.spines[:].set_visible(False)
    ax.get_yticklabels()[0].set_fontweight("bold")
    for y, row in enumerate(effects):
        for x, value in enumerate(row):
            ax.text(x, y, f"{value:+.2f}".replace("-", "−"), ha="center", va="center",
                    color="white" if value < -3.6 else "#162b2b", weight="bold" if y == 0 else "normal")
    fig.text(.035, .94, "The same pathway genes across five endometriosis cohorts", fontsize=17, weight="bold")
    fig.text(.035, .875, "Case minus control, in within-sample percentile points", fontsize=12)
    fig.text(.035, .141, "Positive values indicate higher ranks in cases. Platforms and reference features differ.", fontsize=10)
    fig.text(.035, .086, "¹ No U.S. scan date included both groups. Mumbai flowcell sharing is unknown.", fontsize=10)
    fig.savefig(output, format="svg", metadata={"Date": None, "Title": "Matched pathway genes across five cohorts"})
    plt.close(fig)
    output.write_text("\n".join(line.rstrip() for line in output.read_text().splitlines()) + "\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, default=EXTERNAL / "comparison.svg")
    draw(ap.parse_args().output)
