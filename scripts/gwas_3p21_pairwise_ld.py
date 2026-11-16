"""Recalculate two-marker reference LD from the archived LDpair aggregate receipt."""

import argparse
import hashlib
import json
import math
import re
from pathlib import Path


INPUT = Path(__file__).resolve().parents[1] / "evidence/gwas-3p21-colocalization-audit/pairwise_ld_input.json"


def checked_receipt(source):
    data = source["receipt_text"].encode("utf-8")
    digest = hashlib.sha256(data).hexdigest()
    if len(data) != source["receipt_bytes"] or digest != source["receipt_sha256"]:
        raise ValueError("Aggregate receipt length or SHA-256 differs from the source pin")
    return source["receipt_text"]


def calculate(counts, row_alleles, col_alleles, counted):
    """Use haplotypes (not diploid dosage), with signed LD for the designated alleles."""
    expected = {(a, b) for a in row_alleles for b in col_alleles}
    if len(set(row_alleles)) != 2 or len(set(col_alleles)) != 2 or set(counts) != expected or counted not in expected:
        raise ValueError("Expected all four cells of a biallelic haplotype table")
    if any(type(n) is not int or n < 0 for n in counts.values()):
        raise ValueError("Haplotype counts must be nonnegative integers")
    total = sum(counts.values())
    n1 = sum(n for (a, _), n in counts.items() if a == counted[0])
    n2 = sum(n for (_, b), n in counts.items() if b == counted[1])
    if not (0 < n1 < total and 0 < n2 < total):
        raise ValueError("LD is undefined when either allele frequency is zero or one")

    # delta / total**2 is D. Working from integer margins avoids successive rounding.
    delta = counts[counted] * total - n1 * n2
    var_product = n1 * (total - n1) * n2 * (total - n2)
    if delta >= 0:
        d_bound = min(n1 * (total - n2), (total - n1) * n2)
    else:
        d_bound = min(n1 * n2, (total - n1) * (total - n2))
    r_squared = delta**2 / var_product
    chi_square = total * r_squared
    return {
        "phased_chromosomes": total,
        "counted_allele_1_frequency": n1 / total,
        "counted_allele_2_frequency": n2 / total,
        "d": delta / total**2,
        "signed_r": delta / math.sqrt(var_product),
        "r_squared": r_squared,
        "signed_d_prime": delta / d_bound,
        "pearson_chi_squared_df1_uncorrected": chi_square,
        "pearson_asymptotic_p_value": math.erfc(math.sqrt(chi_square / 2)),
    }


def verify(source):
    receipt = checked_receipt(source)
    first, second = source["marker_order"]
    markers = [first["rsid"], second["rsid"]]
    if source["input_assembly"] != "GRCh37" or source["population"] != "EUR":
        raise ValueError("Expected the GRCh37 European reference")
    if markers != ["rs1352889", "rs13100173"]:
        raise ValueError("Unexpected marker order")
    label = "\n".join(f"{m['rsid']} (chr{m['chromosome']}:{m['position']})" for m in (first, second))
    if f"Query SNPs:\n{label}\n\nEUR Haplotypes:" not in receipt:
        raise ValueError("Mismatch between source markers, coordinates and receipt")
    rows = re.findall(r"^\s*([ACGT])_([ACGT]): (\d+) \([\d.]+\)$", receipt, re.MULTILINE)
    if len(rows) != 4:
        raise ValueError("Expected four LDpair haplotype observations")
    counts = {(a, b): int(n) for a, b, n in rows}
    counted = (first["counted_allele"], second["counted_allele"])
    result = calculate(counts, first["alleles"], second["alleles"], counted)
    if result["phased_chromosomes"] != 2 * source["panel_donors"]:
        raise ValueError("Haplotype count disagrees with the diploid reference-panel size")
    printed_stats = {
        "r_squared": "R2",
        "signed_d_prime": "D'",
        "pearson_chi_squared_df1_uncorrected": "Chi-sq",
    }
    for key, printed in printed_stats.items():
        matched = re.search(rf"^\s*{re.escape(printed)}: ([\d.]+)$", receipt, re.MULTILINE)
        if not matched or round(abs(result[key]), 4) != float(matched[1]):
            raise ValueError(f"{key} fails to reproduce the four-decimal LDpair result")

    crosscheck, = json.loads(checked_receipt(source["ensembl"]))
    ensembl_ids = (crosscheck["variation1"], crosscheck["variation2"], crosscheck["population_name"])
    if ensembl_ids != (*markers, "1000GENOMES:phase_3:EUR"):
        raise ValueError("Ensembl response refers to another pair or population")
    return {
        "assembly": source["input_assembly"],
        "panel": source["reference_panel"],
        "population": source["population"],
        "reference_donors": source["panel_donors"],
        "counted_alleles": dict(zip(markers, counted)),
        "counts_by_haplotype": {"_".join(k): v for k, v in sorted(counts.items())},
        "ldpair_count_calculation": result,
        "ldpair_receipt_sha256": source["receipt_sha256"],
        "ensembl_same_panel_crosscheck": {
            "r_squared": float(crosscheck["r2"]),
            "d_prime": float(crosscheck["d_prime"]),
            "receipt_sha256": source["ensembl"]["receipt_sha256"],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    source = json.loads(parser.parse_args().input.read_text(encoding="utf-8"))
    print(json.dumps(verify(source), indent=2))


if __name__ == "__main__":
    main()
