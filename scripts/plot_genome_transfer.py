"""Plot the conditional genome-score contrasts and eligible Rome scan order."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

EVIDENCE = Path(__file__).resolve().parents[1] / "evidence/genome-transfer-2026-11-15"


def save_figure(fig, filename):
    path = EVIDENCE / filename
    fig.savefig(path, format="svg", metadata={"Date": None})
    # Matplotlib indents multiline SVG paths with trailing spaces.
    path.write_text("\n".join(line.rstrip() for line in path.read_text().splitlines()) + "\n")
    plt.close(fig)


def draw():
    result = json.loads((EVIDENCE / "results.json").read_text())
    swiss = json.loads((EVIDENCE / "swiss.json").read_text())
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "svg.fonttype": "none",
                         "svg.hashsalt": "open-disease-genome-transfer"})
    fig, (forest, order) = plt.subplots(1, 2, figsize=(11.8, 4.5), gridspec_kw={"width_ratios": [1.27, 1]})
    fig.subplots_adjust(left=.21, right=.98, bottom=.28, top=.71, wspace=.39)
    fig.text(.035, .91, "The Finnish expression profile in three external comparisons", fontsize=17,
             weight="bold", color="#182f36")
    fig.text(.035, .835, "Case minus control, using each cohort’s measured genes and the original Finnish weights",
             fontsize=10.5, color="#43575c")
    colors = ["#b45225", "#136c72", "#136c72"]
    rows = [result["rome_within_mixed_scan_dates"], result["mumbai"], swiss]
    forest.axvline(0, color="#829094", lw=1, zorder=0)
    for i, (row, color) in enumerate(zip(rows, colors, strict=True)):
        val = row.get("score_difference_pp", row.get("primary_score_difference_pp"))
        lo, hi = row["constant_shift_95_interval"]
        y = 2 - i
        forest.plot([lo, hi], [y, y], color=color, lw=2.25, solid_capstyle="round")
        forest.plot(val, y, "o", ms=7, color=color)
        forest.text(val, y+.25, f"{val:+.2f}", ha="center", fontsize=9, color=color)
    forest.set(yticks=[2, 1, 0], yticklabels=["Rome whole tissue\n5 cases, 6 controls",
               "Mumbai whole tissue\n4 cases, 3 controls", "Swiss glandular cells\n6 cases, 7 controls"],
               xlim=(-6.1, 6.1), ylim=(-.5, 2.6), xlabel="Weighted percentile points (95% constant-shift interval)")
    forest.set_xticks([-6, -4, -2, 0, 2, 4, 6])
    order.set_title("Rome scanner order within mixed dates", fontsize=10.5, loc="left", pad=12, weight="bold")
    dates = [v for v in result["rome_scan_metadata"]["by_date_chronological"] if v["case"] and v["control"]]
    gaps = result["rome_scan_metadata"]["earliest_case_minus_latest_control_seconds_by_mixed_date"]
    for i, (d, gap) in enumerate(zip(dates, gaps, strict=True)):
        if gap <= 0:
            raise ValueError("The plotted scan separation requires all controls first")
        y = 1 - i
        for x in range(d["control"] + d["case"]):
            case = x >= d["control"]
            order.scatter(x+1, y, color=colors[0] if case else "#536c77",
                          marker="s" if case else "o", s=48)
        minutes, seconds = divmod(gap, 60)
        order.text(.85, y-.26, f"{minutes} min {seconds} s between groups", fontsize=8, color="#55666c")
    order.set(yticks=[1, 0], yticklabels=["Day 1", "Day 2"], xticks=range(1, 7), xlim=(.55, 6.45),
              ylim=(-.53, 1.5), xlabel="Acquisition order among eligible biopsies")
    handles = [Line2D([], [], marker=m, color=c, ls="", label=s)
               for m,c,s in [("o","#536c77","Control"), ("s",colors[0],"Case")]]
    order.legend(handles=handles, loc="upper right", frameon=False, ncol=2, fontsize=8.5,
                 handletextpad=.3, columnspacing=1.0, borderaxespad=0)
    for ax in (forest, order):
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#a6b1b5")
        ax.tick_params(axis="y", length=0, pad=7)
    fig.text(.035, .075, "Rome’s interval conditions on scan date. Diagnosis is ordered within both dates. "
             "Each cohort has its own rank reference.", fontsize=9.1, color="#43575c")
    save_figure(fig, "comparison.svg")
    phase = json.loads((EVIDENCE / "melbourne-phase.json").read_text())
    fig, ax = plt.subplots(figsize=(9.2, 4.6))
    fig.subplots_adjust(left=.10, right=.97, bottom=.23, top=.74)
    fig.text(.045, .91, "Menstrual phase tracks the Finnish score in Melbourne", fontsize=16,
             weight="bold", color="#182f36")
    fig.text(.045, .835, "Mean score by histologic stage and sequencing batch, 257 biopsies",
             color="#43575c", fontsize=10)
    colors = ["#27696b", "#4b4c93", "#b45225"]
    shapes = ["o", "s", "D"]
    for i, (batch, measures) in enumerate(phase["by_batch"].items()):
        cells = list(measures["cells"].values())
        ax.scatter([k+1+(i-1)*.15 for k, v in enumerate(cells) if v["n"]],
                   [v["mean_score"] for v in cells if v["n"]], s=40, marker=shapes[i],
                   color=colors[i], label=batch.capitalize(), zorder=2)
    labels = ["Menstrual", "Early\nproliferative", "Mid\nproliferative", "Late\nproliferative",
              "Early\nsecretory", "Mid\nsecretory", "Late\nsecretory"]
    ax.set(xticks=list(range(1,8)), xticklabels=labels, ylabel="Weighted percentile score")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#a6b1b5")
    ax.set_axisbelow(True)
    ax.grid(axis="y", color="#e1e7e6")
    ax.tick_params(length=0, pad=5)
    ax.legend(loc="lower right", frameon=False, ncol=3, fontsize=9)
    fig.text(.045, .063, "Mid-secretory minus mid-proliferative: +1.94 points (95% CI +1.79 to +2.09), adjusted for batch.",
             color="#43575c", fontsize=9.5)
    save_figure(fig, "melbourne-phase.svg")


if __name__ == "__main__":
    draw()
