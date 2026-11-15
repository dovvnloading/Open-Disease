"""Exploratory sensitivity for nearby reported markers when an Indian lead is unavailable."""

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path
from india_gwas_portability import (
    PANEL,
    file_hash,
    PDF_SHA,
    load_indian_panel,
    load_validation,
)


def secondary(pdf, gwas, primary_table):
    if file_hash(pdf, "sha256") != PDF_SHA:
        raise ValueError("Unexpected discovery PDF")
    primary = load_indian_panel(pdf)
    with primary_table.open() as stream:
        missing = [
            r["rsid"]
            for r in csv.DictReader(stream, delimiter="\t")
            if r["status"] == "absent"
        ]
    text = subprocess.check_output(
        ["pdftotext", "-layout", str(pdf), "-"], text=True
    ).replace("·", ".")
    s2 = text.split("Supplementary Table S2: List of variants", 1)[1].split(
        "Supplementary Table S3:", 1
    )[0]
    selected = {}
    for line in s2.splitlines():
        fields = line.split()
        if len(fields) != 9 or not fields[0].startswith("rs") or fields[0] in PANEL:
            continue
        rsid, chrom, pos, a1, a2, _, odds, _, p = fields
        nearby = [
            lead
            for lead in missing
            if primary[lead]["chromosome"] == chrom
            and abs(int(primary[lead]["pos"]) - int(pos)) <= 250_000
        ]
        if nearby:
            if len(nearby) != 1 or rsid in selected:
                raise ValueError("Unexpected nonunique secondary marker")
            selected[rsid] = (nearby[0], chrom, pos, a1, a2, float(odds), float(p))
    validation, _ = load_validation(gwas, selected)
    output = []
    for rsid, (lead, chrom, pos, a1, a2, odds, p) in selected.items():
        row = validation.get(rsid)
        out = dict(
            lead=lead,
            rsid=rsid,
            discovery_risk_allele=a1 if odds > 1 else a2,
            india_risk_or=max(odds, 1 / odds),
            india_p=p,
            status="absent",
        )
        if row:
            if (row["chromosome"], row["base_pair_location"]) != (chrom, pos) or {
                a1,
                a2,
            } != {row["effect_allele"], row["other_allele"]}:
                raise ValueError("Regional marker mismatch")
            out["external_p"] = float(row["p_value"])
            out["external_bonferroni_p"] = min(1, len(selected) * out["external_p"])
            if {a1, a2} in ({"A", "T"}, {"G", "C"}):
                out["status"] = "strand_ambiguous"
            else:
                sign = 1 if row["effect_allele"] == out["discovery_risk_allele"] else -1
                b = sign * float(row["beta"])
                se = float(row["standard_error"])
                out.update(
                    status="aligned",
                    external_or=math.exp(b),
                    lower_95=math.exp(b - 1.96 * se),
                    upper_95=math.exp(b + 1.96 * se),
                )
        output.append(out)
    return output


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("supplement", "gwas", "primary", "output"):
        p.add_argument(name, type=Path)
    args = p.parse_args()
    results = secondary(args.supplement, args.gwas, args.primary)
    args.output.write_text(json.dumps(results, indent=2) + "\n")
