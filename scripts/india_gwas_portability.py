"""Test the 21 exact lead variants from Anand et al. against GWAS Catalog GCST90205183.

Requires two locally supplied aggregate files and pdftotext on PATH. No network access.
Only analysis-derived variant-level summaries and diagnostic counts are written.
"""

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
import subprocess
from statistics import NormalDist
from pathlib import Path

PANEL = [
    "rs872595",
    "rs145379000",
    "rs114511502",
    "rs61829563",
    "rs6710746",
    "rs1861768",
    "rs2679859",
    "rs7563305",
    "rs12631821",
    "rs1005940",
    "rs58794955",
    "rs9401193",
    "rs56991723",
    "rs11017546",
    "rs1599469",
    "rs76975203",
    "rs35236165",
    "rs12867495",
    "rs28698234",
    "rs140109794",
    "rs2838083",
]
CONTROLS = {"rs10917151": "A", "rs10122243": "T"}
PDF_SHA = "f3b072f47d0a6ebd8588fdc015b4d56bf64c817b6cbfb320a999663156431469"
GWAS_MD5 = "b5a0ce654a6e738b056048aee4296a4a"
FAMILY = len(PANEL)
SIMULTANEOUS_Z = NormalDist().inv_cdf(1 - 0.05 / (2 * FAMILY))


def file_hash(path, method):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, method).hexdigest()


def load_indian_panel(pdf):
    if file_hash(pdf, "sha256") != PDF_SHA:
        raise ValueError("Unexpected discovery supplement SHA-256")
    raw = subprocess.check_output(["pdftotext", "-layout", str(pdf), "-"], text=True)
    raw = raw.replace("·", ".").replace("chr7:7432100 ", "rs56991723 ")
    s2 = raw.split("Supplementary Table S2: List of variants", 1)[1].split(
        "Supplementary Table S3:", 1
    )[0]
    s4 = raw.split("Supplementary Table S4: Allele frequency variation", 1)[1].split(
        "Supplementary Table S5:", 1
    )[0]
    roles = {}
    for line in s2.splitlines():
        fields = line.split()
        if len(fields) != 9 or fields[0] not in PANEL:
            continue
        rsid, chrom, pos, a1, a2, freq, odds, se, p = fields
        odds, freq = float(odds), float(freq)
        if rsid in roles:
            raise ValueError(f"Duplicate discovery estimate: {rsid}")
        roles[rsid] = dict(
            chromosome=chrom,
            pos=pos,
            risk=(a1 if odds > 1 else a2),
            alleles={a1, a2},
            india_freq=(freq if odds > 1 else 1 - freq),
            india_or=max(odds, 1 / odds),
            india_se=float(se),
            india_p=float(p),
        )
    for line in s4.replace("chr7:74321005:G:A", "rs56991723").splitlines():
        fields = line.split()
        if fields and fields[0] in PANEL:
            # Nine trailing numerical fields: ECGRI, AFR, AMR, EAS, EUR, SAS, min, max, range.
            freqs = [float(v) for v in fields[-9:]]
            r = roles[fields[0]]
            if abs(freqs[0] - r["india_freq"]) > 0.0011:
                raise ValueError(f"S2/S4 risk orientation disagreement: {fields[0]}")
            r["europe_risk_freq"] = freqs[4]
    if set(roles) != set(PANEL) or any(
        "europe_risk_freq" not in r for r in roles.values()
    ):
        raise ValueError("Discovery family identity mismatch")
    return roles


def load_validation(gwas, additional=()):
    if file_hash(gwas, "md5") != GWAS_MD5:
        raise ValueError("Unexpected validation GWAS MD5")
    selected = set(PANEL) | CONTROLS.keys() | set(additional)
    wanted = {name.encode() for name in selected}
    records = {}
    with gzip.open(gwas, "rb") as stream:
        columns = stream.readline().decode().split()
        if columns[-1] != "rsid":
            raise ValueError("Unexpected GWAS header")
        for n, line in enumerate(stream, start=1):
            if line.rstrip().rsplit(b"\t", 1)[-1] not in wanted:
                continue
            row = dict(zip(columns, line.decode().split(), strict=True))
            rsid = row["rsid"]
            if rsid in records:
                raise ValueError(f"Multiple validation rows: {rsid}")
            records[rsid] = row
    return records, n


def orient(panel, row):
    ea, oa = row["effect_allele"], row["other_allele"]
    if (row["chromosome"], row["base_pair_location"]) != (
        panel["chromosome"],
        panel["pos"],
    ):
        return "coordinate_mismatch", None
    if {ea, oa} != panel["alleles"]:
        return "allele_mismatch", None
    sign = 1 if panel["risk"] == ea else -1
    if {ea, oa} in ({"A", "T"}, {"C", "G"}):
        freq = float(row["effect_allele_frequency"])
        risk_freq = freq if sign == 1 else 1 - freq
        india, eur = panel["india_freq"], panel["europe_risk_freq"]
        if (
            abs(india - 0.5) < 0.1
            or abs(eur - 0.5) < 0.1
            or (india - 0.5) * (eur - 0.5) < 0
            or abs(risk_freq - eur) > 0.1
        ):
            return "strand_ambiguous", None
    return "aligned", sign


def analyze(pdf, gwas):
    discovery = load_indian_panel(pdf)
    records, n_checked = load_validation(gwas)
    table = []
    for rsid in PANEL:
        panel, row = discovery[rsid], records.get(rsid)
        output = {
            "rsid": rsid,
            "risk_allele": panel["risk"],
            "india_or": round(panel["india_or"], 5),
            "india_p": panel["india_p"],
            "status": "absent",
        }
        if row:
            status, sign = orient(panel, row)
            output["status"] = status
            if status in {"aligned", "strand_ambiguous"}:
                output["external_two_sided_p"] = float(row["p_value"])
                output["external_bonferroni_p_21"] = min(
                    1, FAMILY * float(row["p_value"])
                )
            if sign is not None:
                beta, se = sign * float(row["beta"]), float(row["standard_error"])
                output.update(
                    external_or=math.exp(beta),
                    lower_95=math.exp(beta - 1.96 * se),
                    upper_95=math.exp(beta + 1.96 * se),
                    external_log_or=beta,
                    external_se=se,
                    lower_family_95=math.exp(beta - SIMULTANEOUS_Z * se),
                    upper_family_95=math.exp(beta + SIMULTANEOUS_Z * se),
                )
        table.append(output)
    aligned = [r for r in table if r["status"] == "aligned"]
    positive = sum(r["external_log_or"] > 0 for r in aligned)
    nonzero = sum(r["external_log_or"] != 0 for r in aligned)
    if nonzero:
        sign_p = (
            sum(math.comb(nonzero, i) for i in range(positive, nonzero + 1))
            / 2**nonzero
        )
    else:
        sign_p = None
    controls = {}
    for rsid, risk in CONTROLS.items():
        row = records[rsid]
        if risk not in {row["effect_allele"], row["other_allele"]}:
            raise ValueError(f"Positive control allele mismatch: {rsid}")
        beta = float(row["beta"]) * (1 if row["effect_allele"] == risk else -1)
        controls[rsid] = {
            "alleles": [row["effect_allele"], row["other_allele"]],
            "risk_allele": risk,
            "external_or": math.exp(beta),
            "external_p": float(row["p_value"]),
        }
    summary = dict(
        tested_family=FAMILY,
        validation_rows_scanned=n_checked,
        validation_present=sum(r["status"] != "absent" for r in table),
        aligned=len(aligned),
        concordant=positive,
        nonzero_log_or=nonzero,
        binomial_one_sided=sign_p,
        same_direction_bonferroni_pass=sum(
            r["external_log_or"] > 0 and r["external_bonferroni_p_21"] <= 0.05
            for r in aligned
        ),
        family_ci_within_10_percent=sum(
            r["lower_family_95"] > 1 / 1.1 and r["upper_family_95"] < 1.1
            for r in aligned
        ),
        family_upper_below_india_point=sum(
            r["upper_family_95"] < r["india_or"] for r in aligned
        ),
        positive_controls=controls,
    )
    return table, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("supplement_pdf", type=Path)
    parser.add_argument("gwas_gz", type=Path)
    parser.add_argument("output_directory", type=Path)
    args = parser.parse_args()
    rows, summary = analyze(args.supplement_pdf, args.gwas_gz)
    args.output_directory.mkdir(parents=True, exist_ok=True)
    with (args.output_directory / "variants.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(dict.fromkeys(k for r in rows for k in r)),
            delimiter="\t",
            lineterminator="\n",
            restval="NA",
        )
        writer.writeheader()
        writer.writerows(rows)
    (args.output_directory / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
