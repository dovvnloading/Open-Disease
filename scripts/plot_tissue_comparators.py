"""Plot the fixed panel's six aggregate tissue contrasts on one percentile-point scale."""
import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt

from plot_lesion_panel import EVIDENCE, INK, MUTED, aggregate_inputs


def figure(evidence, adjacent_json, output):
    _, _, upper = aggregate_inputs(evidence)
    source = json.loads(adjacent_json.read_text())
    lower = [source[k] for k in ("core_minus_uterus","adjacent_minus_uterus","core_minus_adjacent")]
    if not math.isclose(lower[0]["score_pp"] - lower[1]["score_pp"],lower[2]["score_pp"],abs_tol=1e-9):
        raise ValueError("Same-cohort tissue contrasts do not reconcile")
    effects = [x["effect_pp"] for x in upper] + [x["score_pp"] for x in lower]
    higher = [x["higher"] for x in upper] + [x["positive"] for x in lower]
    totals = [x["pairs"] for x in (*upper,*lower)]
    low, high = lower[2]["exact_95pct_constant_shift_pp"]

    plt.rcParams.update({"font.family":"DejaVu Sans","svg.fonttype":"none",
                         "svg.hashsalt":"open-disease-tissue-comparators-v1"})
    fig = plt.figure(figsize=(10.5,7.0),facecolor="white")
    ax = fig.add_axes((0,0,1,1))
    ax.set(xlim=(0,1050),ylim=(700,0))
    ax.axis("off")

    def text(x,y,label,size=11,color=INK,weight="normal",align="left"):
        ax.text(x,y,label,fontsize=size,color=color,fontweight=weight,va="center",ha=align)

    def point(value):
        return 382 + 8.1*value

    def signed(value):
        return f"{value:+.2f}".replace("-","−")

    text(42,39,"The five-gene score depends on the tissue comparator",19,weight="bold")
    text(42,76,"Mean within-patient change in signed expression rank, in percentile points",12,MUTED)
    text(42,132,"Lesion minus uterine tissue",14,weight="bold")
    text(42,156,"Three separately recruited cohorts",10.5,MUTED)
    text(850,141,"Mean change",10,MUTED,align="center")
    text(971,141,"Higher score",10,MUTED,align="center")

    rows = ("Hull  ·  Affymetrix","CodeLink  ·  microarray","Peking  ·  RNA-seq",
            "Lesion core minus uterus","Adjacent tissue minus uterus","Lesion core minus adjacent")
    heights = (198,243,288,417,462,507)
    for start,stop in ((176,309),(396,530)):
        for tick in (-10,0,10,20,30,40,50):
            ax.plot((point(tick),point(tick)),(start,stop),lw=1.2 if tick==0 else .6,
                    color="#71858C" if tick==0 else "#DFE6E8",zorder=0)
    ax.plot((42,1010),(329,329),lw=.8,color="#BBC8CB")
    text(42,355,"Peritoneal lesion, adjacent tissue and uterus",14,weight="bold")
    text(42,379,"Three contrasts in the same six patients  ·  single-cell RNA-seq aggregated by biopsy",10.5,MUTED)
    for i,(name,y,effect,npos,n) in enumerate(zip(rows,heights,effects,higher,totals,strict=True)):
        text(42,y,name,11,weight="bold" if i==5 else "normal")
        text(850,y,signed(effect),12,weight="bold",align="center")
        text(971,y,f"{npos} of {n}",11,MUTED,align="center")
        color = "#255F77" if i<3 else "#945531"
        if i==5:
            ax.plot((point(low),point(high)),(y,y),lw=2,color=color,zorder=2)
            for endpoint in (low,high):
                ax.plot((point(endpoint),point(endpoint)),(y-5,y+5),lw=1.4,color=color)
        ax.scatter(point(effect),y,s=60,facecolor=color,edgecolor="white",linewidth=.9,zorder=3)
    for tick in (-10,0,10,20,30,40,50):
        text(point(tick),552,str(tick).replace("-","−"),9.5,MUTED,align="center")
    text(543,580,"Difference in score (percentile points)",10.5,MUTED,align="center")
    ax.plot((42,1010),(605,605),lw=.8,color="#BBC8CB")
    text(42,625,"Whisker: exact 95% constant-shift interval for lesion core minus adjacent tissue "
         f"({signed(low)} to {signed(high)}).",10,MUTED)
    text(42,648,"All other dots show point estimates alone. “Higher score” counts patients with a positive paired difference.",10,MUTED)
    text(42,675,"GEO: GSE11691, GSE5108, GSE282532 and GSE179640.",9.5,MUTED)
    fig.savefig(output,format="svg",metadata={"Date":None,"Title":"Paired five-gene panel scores by tissue comparator",
                "Description":"Six separate within-person comparisons on a common scale. The interval "
                "is available for the direct peritoneal lesion-core versus adjacent-tissue comparison."})
    plt.close(fig)
    lines = output.read_text(encoding="utf-8").splitlines()
    output.write_text("\n".join(line.rstrip() for line in lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data",type=Path,default=EVIDENCE)
    parser.add_argument("--adjacent",type=Path,default=EVIDENCE.parent/"adjacent-peritoneum"/"results.json")
    parser.add_argument("--out",type=Path,default=EVIDENCE/"tissue-comparators.svg")
    args = parser.parse_args()
    figure(args.data,args.adjacent,args.out)
