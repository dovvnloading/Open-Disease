"""Render the five frozen genes in three external paired studies from aggregate JSON."""

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgb
from matplotlib.patches import Rectangle

EVIDENCE = Path(__file__).resolve().parents[1] / "evidence" / "lesion-transfer"
SYMBOLS = {
    55026: "TMEM255A",
    9077: "DIRAS3",
    2568: "GABRP",
    6299: "SALL1",
    3512: "JCHAIN",
}
TITLES = ("Hull", "CodeLink", "Peking")
ASSAYS = ("Affymetrix array", "CodeLink array", "RNA-seq")
ACCESSIONS = ("GSE11691", "GSE5108", "GSE282532")
INK, MUTED = "#172F38", "#4A616B"
CMAP = LinearSegmentedColormap.from_list(
    "paired_direction", ["#075975", "#81B9CC", "#F7F7F4", "#E5B58E", "#A04B25"], N=257
)


def aggregate_inputs(directory):
    def read(filename):
        return json.loads((directory / filename).read_text())

    lock = {g["entrez"]: g["direction"] for g in read("panel.json")}
    hull = read("results.json")["primary"]
    code = read("codelink.json")
    peking = read("peking-rnaseq.json")
    hull_by_id = {g["entrez"]: g for g in hull["panel_genes"]}
    if (
        set(lock) != set(SYMBOLS)
        or set(hull_by_id) != set(lock)
        or set(peking["panel_genes"]) != set(SYMBOLS.values())
    ):
        raise ValueError("Aggregate panel membership mismatch")
    effects = {}
    for key, symbol in SYMBOLS.items():
        a, b, c = (
            hull_by_id[key],
            code["genes"][str(key)],
            peking["panel_genes"][symbol],
        )
        if a["direction"] != lock[key] or c["direction"] != lock[key]:
            raise ValueError("Aggregate gene direction changed")
        effects[key] = (a["validation_pp"], b["effect_pp"], c["difference_pp"])
    summaries = (hull["panel_pp"], code["panel"], peking["primary_pp"])
    for i, cohort in enumerate(summaries):
        signed = sum(lock[key] * row[i] for key, row in effects.items()) / len(effects)
        if not math.isclose(signed, cohort["effect_pp"], abs_tol=1e-9):
            raise ValueError("Panel score differs from mean signed gene estimate")
    return lock, effects, summaries


def contrast_ink(color):
    srgb = to_rgb(color)
    linear = [
        (c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4) for c in srgb
    ]
    luminance = sum(c * w for c, w in zip(linear, (0.2126, 0.7152, 0.0722)))
    return "#FFFFFF" if luminance < 0.21 else "#15272E"


def figure(directory, output):
    lock, effects, summaries = aggregate_inputs(directory)
    normal = Normalize(-70, 70)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "svg.fonttype": "none",
            "svg.hashsalt": "open-disease-lesion-panel-v1",
        }
    )
    fig = plt.figure(figsize=(10.8, 7.15), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set(xlim=(0, 1080), ylim=(715, 0))
    ax.axis("off")

    def label(x, y, s, size=12, color=INK, weight="normal", ha="left"):
        return ax.text(
            x, y, s, ha=ha, va="center", fontsize=size, color=color, weight=weight
        )

    label(
        48,
        41,
        "Five-gene expression in lesions and matched uterine tissue",
        19,
        weight="bold",
    )
    label(
        48,
        80,
        "Mean within-person difference in expression rank, in percentile points",
        13,
        MUTED,
    )
    label(48, 127, "Gene", 11, weight="bold")
    label(184, 127, "Frozen direction", 10, MUTED)
    xs = [412, 663, 914]
    for i, (x, s) in enumerate(zip(xs, summaries)):
        label(x, 120, TITLES[i], 15, weight="bold", ha="center")
        label(x, 145, f"{ASSAYS[i]}  ·  {s['pairs']} pairs", 10.3, MUTED, ha="center")
    ax.plot((48, 1033), (163, 163), color="#B9C8CD", lw=0.8)
    # Group genes by the frozen direction, sorting symbols within each group.
    for j, key in enumerate(sorted(SYMBOLS, key=lambda k: (-lock[k], SYMBOLS[k]))):
        y = 194 + 50 * j
        label(48, y, SYMBOLS[key], 13, weight="bold")
        label(184, y, "Higher" if lock[key] == 1 else "Lower", 11, MUTED)
        for x, value in zip(xs, effects[key]):
            fill = CMAP(normal(value))
            ax.add_patch(
                Rectangle((x - 119, y - 23), 238, 46, facecolor=fill, edgecolor="none")
            )
            label(
                x,
                y,
                f"{value:+.1f}".replace("-", "−"),
                14,
                contrast_ink(fill),
                weight="bold",
                ha="center",
            )
    label(48, 438, "Lesion minus uterine rank", 10, MUTED)
    for j in range(70):
        ax.add_patch(
            Rectangle(
                (328 + 5 * j, 428), 5.1, 13, facecolor=CMAP(j / 69), edgecolor="none"
            )
        )
    for x, s in [(328, "−70"), (503, "0"), (678, "+70")]:
        label(x, 456, s, 9, MUTED, ha="center")
    ax.plot((48, 1033), (480, 480), color="#B9C8CD", lw=0.8)
    label(48, 510, "Signed panel score", 14, weight="bold")
    label(48, 542, "Pairs with higher score", 11, MUTED)
    for x, s in zip(xs, summaries):
        label(x, 510, f"{s['effect_pp']:+.1f}", 17, weight="bold", ha="center")
        label(x, 542, f"{s['higher']} of {s['pairs']}", 12, MUTED, ha="center")
    ax.plot((48, 1033), (567, 567), color="#B9C8CD", lw=0.8)
    label(
        48,
        593,
        "Cell values retain each gene’s measured direction. The panel score reverses the sign for genes",
        11,
        MUTED,
    )
    label(
        48,
        615,
        "frozen as lower, then averages all five changes. Ranks use each assay’s own feature set.",
        11,
        MUTED,
    )
    label(48, 660, "GEO series, left to right", 10, weight="bold")
    label(258, 660, "  ·  ".join(ACCESSIONS), 10, MUTED)
    fig.savefig(
        output,
        format="svg",
        metadata={
            "Date": None,
            "Title": "Five-gene expression in lesions and matched uterine tissue",
            "Description": "Descriptive aggregate gene rank shifts in three external paired cohorts. "
            "Cell color encodes lesion minus uterine difference in percentile points.",
        },
    )
    plt.close(fig)
    lines = output.read_text(encoding="utf-8").splitlines()
    output.write_text(
        "\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=EVIDENCE)
    parser.add_argument("--out", type=Path, default=EVIDENCE / "panel-comparison.svg")
    args = parser.parse_args()
    figure(args.data, args.out)
