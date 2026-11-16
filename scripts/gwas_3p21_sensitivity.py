"""Check the frozen single-signal gates using the alternative 2023 disease input."""

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

from gwas_3p21_disease_slice import INTERVAL, SLICE_SHA256
from gwas_3p21_rahmioglu_workbook import GENES

START, END = INTERVAL[1:]
# Two GRCh38 insertion bases have no one-to-one GRCh37 source base.
UNMAPPED = {50_214_403, 50_269_326}
TISSUES = ("Uterus", "Colon_Sigmoid")
DISEASE_FILE = "clinical-metadata/gwas_region.tsv"
REFERENCE_FILE = "GRCh38.fa.txt"
CSET_FILE = "rahmioglu_s11_bsn_member_rsids.txt"
DIGESTS = {
    DISEASE_FILE: SLICE_SHA256,
    REFERENCE_FILE: "973fa89f78cecb010e127f09eed36c995be0217bf81eff4da500f0ffca42059c",
    CSET_FILE: "468d2dffb084718971627a063fe28ad511d047ca01a9d514ae68e533a045c15b",
    "Uterus.tsv": "ab140f050bf5186e75e27d857068a22fa4e7ecb0585b5eda9e00c5c96790f990",
    "Colon_Sigmoid.tsv": "239fe0f800b9002681a68243b49ef86728e7c4c6eac636794cb4a207387a2153",
}
REQUIRED = ("rs1352889", "rs13100173", "rs6804970")
BASES = set("ACGT")


def table(path):
    with path.open() as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def valid_numbers(effect, se, p, frequency):
    return (math.isfinite(float(effect)) and math.isfinite(float(se)) and float(se) > 0
            and 0 <= float(p) <= 1 and 0 < float(frequency) < 1)


def canonical_snv(pos, a1, a2, reference):
    """Resolve REF/ALT against the upper-strand reference sequence; no strand guessing."""
    if pos < START or pos > END or pos in UNMAPPED:
        return None, "unmapped"
    alleles = {a1, a2}
    if len(alleles) != 2 or not alleles <= BASES:
        return None, "not_biallelic_snv"
    ref = reference[pos - START]
    if ref not in alleles:
        return None, "reference_allele_mismatch"
    alt = a2 if ref == a1 else a1
    return (pos, ref, alt), None


def disease_table(rows, reference):
    positions = Counter(int(r["base_pair_location"]) for r in rows)
    filtered, markers, counts = {}, {}, Counter()
    for r in rows:
        pos = int(r["base_pair_location"])
        assert r["chromosome"] == "3"
        key, reason = canonical_snv(pos, r["effect_allele"], r["other_allele"], reference)
        if reason is None and positions[pos] != 1:
            reason = "duplicate_or_multiallelic_position"
        if reason is None and {key[1], key[2]} in ({"A", "T"}, {"C", "G"}):
            reason = "unresolved_palindrome"
        if reason is None and not valid_numbers(r["beta"], r["standard_error"], r["p_value"], r["effect_allele_frequency"]):
            reason = "invalid_effect_se_p_or_frequency"
        counts[reason or "eligible"] += 1
        if reason is not None:
            continue
        beta = float(r["beta"])
        filtered[key] = (r["rsid"], beta if key[2] == r["effect_allele"] else -beta)
        markers[r["rsid"]] = key
    return filtered, markers, counts


def joint_row(r, reference, positions, disease):
    """Return at most one QC exit, the normalized key and the aligned disease beta."""
    pos = int(r["position"])
    key, reason = canonical_snv(pos, r["ref"], r["alt"], reference)
    if reason is None and key[1:] != (r["ref"], r["alt"]):
        reason = "qtl_ref_mismatch"
    if reason is None and positions[pos] != 1:
        reason = "duplicate_or_multiallelic_gene_position"
    if reason is None and {key[1], key[2]} in ({"A", "T"}, {"C", "G"}):
        reason = "unresolved_palindrome"
    if reason is None and (not valid_numbers(r["beta"], r["se"], r["pvalue"], r["maf"])
                           or float(r["maf"]) > 0.5 or int(r["an"]) <= 0):
        reason = "invalid_effect_se_p_maf_or_an"
    if reason is None and key not in disease:
        reason = "no_eligible_disease_allele_match"
    if reason:
        return reason, None
    return "retained", key


def evaluate_gene(records, reference, disease, markers, watch, credible):
    positions = Counter(int(r["position"]) for r in records)
    counts = Counter()
    retained, covered_cset, observed_cset = set(), set(), set()
    # The strongest measured cis-eQTL is identified before intersection with the GWAS,
    # including indels. Restricting its definition to the intersection hides peak loss.
    best_p = min(float(r["pvalue"]) for r in records)
    peaks, all_peaks_retained, focus_rows = [], True, {}
    for r in records:
        reason, key = joint_row(r, reference, positions, disease)
        counts[reason] += 1
        if float(r["pvalue"]) == best_p:
            peaks.append({"variant": r["variant"], "rsid": r["rsid"], "status": reason})
            all_peaks_retained &= reason == "retained"
        if r["rsid"] in credible:
            observed_cset.add(r["rsid"])
        if key is None:
            continue
        retained.add(key)
        rsid, aligned_disease_beta = disease[key]
        if rsid in credible:
            covered_cset.add(rsid)
        if rsid in watch and r["gene_id"] in GENES:
            focus_rows[rsid] = {"ref": key[1], "alt_effect_allele": key[2],
                                "disease_beta_alt": aligned_disease_beta, "qtl_beta_alt": float(r["beta"]),
                                "qtl_se": float(r["se"]), "qtl_p": float(r["pvalue"]), "qtl_maf": float(r["maf"])}
    missing = [rsid for rsid in watch if markers.get(rsid) not in retained]
    return {
        "source_rows": len(records), "qc_rows": dict(sorted(counts.items())),
        "credible_set_observed_by_qtl_rsid": len(observed_cset),
        "credible_set_retained_by_disease_rsid": len(covered_cset),
        "missing_required_markers": missing, "minimum_qtl_p_before_intersection": best_p,
        "qtl_peak": peaks, "all_tied_qtl_peaks_retained": bool(all_peaks_retained),
        "eligible": not missing and bool(all_peaks_retained),
        **({"focus_marker_effects":focus_rows} if focus_rows else {}),
    }


def evaluate(root):
    for filename, digest in DIGESTS.items():
        with (root / filename).open("rb") as f:
            if hashlib.file_digest(f, "sha256").hexdigest() != digest:
                raise ValueError(f"Digest mismatch: {filename}")
    ref = (root / REFERENCE_FILE).read_text()
    if len(ref) != END - START + 1:
        raise ValueError("Reference does not cover the full inclusive interval")
    credible = set((root / CSET_FILE).read_text().splitlines())
    rows = list(table(root / DISEASE_FILE))
    min_p = min(float(r["p_value"]) for r in rows)
    disease_peaks = sorted({r["rsid"] for r in rows if float(r["p_value"]) == min_p})
    watch = list(dict.fromkeys(disease_peaks + list(REQUIRED)))
    disease, markers, counts = disease_table(rows, ref)
    tissues = {}
    for t in TISSUES:
        by_gene = defaultdict(list)
        for r in table(root / (t + ".tsv")):
            if r["gene_id"] != r["molecular_trait_id"] or r["chromosome"] != "3":
                raise ValueError("Unexpected trait or chromosome")
            by_gene[r["gene_id"]].append(r)
        genes = {g: evaluate_gene(records, ref, disease, markers, watch, credible)
                 for g, records in sorted(by_gene.items())}
        counts_marker = {rs: sum(rs not in g["missing_required_markers"] for g in genes.values()) for rs in watch}
        missing_focus = [symbol for ensg, symbol in GENES.items() if ensg not in genes]
        num_eligible = sum(g["eligible"] for g in genes.values())
        if num_eligible:
            raise ValueError("A gene cleared the gates. Run the frozen statistical comparison before closing the audit.")
        tissues[t] = {"source_rows": sum(g["source_rows"] for g in genes.values()), "genes_assayed": len(genes),
                      "eligible_genes": num_eligible, "prespecified_genes_absent": missing_focus,
                      "genes_retaining_required_marker":counts_marker,
                      "genes_retaining_all_tied_qtl_peaks":sum(g["all_tied_qtl_peaks_retained"] for g in genes.values()),
                      "genes_with_zero_joint_rows":sum(g["qc_rows"].get("retained", 0)==0 for g in genes.values()),
                      "genes": genes}
    return {
        "analysis":"Eligibility for one-causal-variant sensitivity, GCST90205183 vs original GTEx v8",
        "input_sha256": DIGESTS, "GRCh38_region_chr3_inclusive":[START,END],
        "GWAS_raw_rows":len(rows), "GWAS_QC_rows":dict(sorted(counts.items())),
        "published_2023_BSN_99pct_credible_set":{"total":len(credible),
             "rsids_in_raw_2023_changed_GWAS":len(credible & {r["rsid"] for r in rows}),
             "rsids_after_GWAS_QC":len(credible & set(markers))},
        "GWAS_unthresholded_regional_peak_rsids":disease_peaks, "regional_peak_p":min_p,
        "required_markers":watch, "prespecified_genes":GENES,
        "tissues":tissues, "models_fitted":0, "prior_sensitivity_runs":0,
        "H0_H1_H2_H3_H4":None, "sigmoid_colon_HYAL3_control":"Failed the same coverage gate",
        "VEZT_control":"Independent endometrial source does not supply harmonizable full signed statistics",
        "interpretation":"All genes failed the frozen coverage rule. No colocalization probabilities were calculated. Relative effect directions at jointly measured individual variants do not estimate sharing.",
    }


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("local_inputs", type=Path)
    print(json.dumps(evaluate(cli.parse_args().local_inputs), indent=2, allow_nan=False))
