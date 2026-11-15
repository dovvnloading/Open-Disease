"""Test an externally defined uterine TWAS direction in endometrial cohorts."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from fractions import Fraction
import json
from pathlib import Path

import numpy as np
from openpyxl import load_workbook

from bulk_gse153739 import load_mumbai
from cross_cohort_discovery import assignment_design, holm
from genome_endometdb import load_finland
from genome_rome import load_rome
from genome_signature_transfer import DISCOVERY, DISCOVERY_SHA
from open_disease.geo_microarray import sha256_file
from rome_scan_design import scan_design
from swiss_epithelial_mapping import HGNC_FILE, HGNC_SHA256
from swiss_global_interval import allocation_contrasts, shift_limits
from swiss_genome_signature import load_glandular_profiles

BOOK_SHA = "c021b6c20c668644939fd83381307a2c8cf08f18d95156b33e65fb0301997815"
PLAN = "ba37b1ae3518d0750a644de31f8fd779fd1434d8"


def selected_genes(book, hgnc):
    if sha256_file(book) != BOOK_SHA or sha256_file(hgnc) != HGNC_SHA256:
        raise ValueError("Genetic or nomenclature source differs from its pinned release")
    wb = load_workbook(book, read_only=True, data_only=True)
    try:
        table = iter(wb["Supplementary Table 13"].values)
        next(table)
        next(table)
        fields = next(table)
        models = []
        for values in table:
            row = dict(zip(fields, values, strict=True))
            p, z = row["Uterus p-value"], row["Uterus z-score"]
            if (p, z) == ("NA", "NA"):
                continue
            if not (0 < p <= 1 and np.isfinite(z)):
                raise ValueError("Invalid uterine prediction statistic")
            models.append((p, z, row["Gene name"], row["Gene transcript"].split(".")[0]))
    finally:
        wb.close()
    if len(models) != 2533 or len({gene for _, _, _, gene in models}) != len(models):
        raise ValueError("Uterine model universe changed")
    with hgnc.open(newline="", encoding="utf-8") as stream:
        approved = {r["symbol"]: r for r in csv.DictReader(stream, delimiter="\t") if r["status"] == "Approved"}
    selected = []
    for p, z, symbol, stable in sorted(models):
        if p > .05 / len(models):
            continue
        entry = approved[symbol]
        if stable != entry["ensembl_gene_id"] or not entry["entrez_id"].isdigit() or not entry["location"]:
            raise ValueError(f"Unresolved gene mapping for {symbol}")
        selected.append(dict(symbol=symbol, entrez=int(entry["entrez_id"]), band=entry["location"],
                             uterus_z=z, uterus_p=p, direction=int(np.sign(z))))
    return len(models), selected


def score_cohort(candidates, genes, matrix, cases, strata):
    """Give each represented band equal weight, then each gene within it equal weight."""
    positions = {int(gene): j for j, gene in enumerate(genes)}
    if len(positions) != len(genes):
        raise ValueError("Repeated cohort gene ID")
    cols = [positions[row["entrez"]] for row in candidates]
    raw = matrix[:, cols]
    if not candidates or not np.isfinite(raw).all():
        raise ValueError("Missing or nonfinite score inputs")
    signed = raw * np.array([row["direction"] for row in candidates])
    bands = [row["band"] for row in candidates]
    distinct = sorted(set(bands))
    grouped = np.column_stack([signed[:, np.array(bands) == band].mean(axis=1) for band in distinct])
    score = grouped.mean(axis=1)
    observed, allocations, _ = assignment_design(cases, strata)
    effect = float(observed @ score)
    reference = allocations @ score
    extreme = int((reference >= effect - 1e-10).sum())
    both = int((abs(reference) >= abs(effect) - 1e-10).sum())
    # Second calculation enumerates stratum combinations on rational donor scores.
    direct, reference_exact = allocation_contrasts(
        score.tolist(), cases.astype(int).tolist(), strata.tolist())
    tol = Fraction(1, 10**10)
    if (abs(effect - float(direct)) > 1e-10
            or sum(a + tol >= direct for a, _ in reference_exact) != extreme
            or sum(abs(a) + tol >= abs(direct) for a, _ in reference_exact) != both):
        raise ValueError("Score or allocation tail failed independent enumeration")
    interval = shift_limits(direct, reference_exact)[:2]
    native = observed @ raw
    band_effects = observed @ grouped
    if not np.isclose(effect, band_effects.mean(), atol=1e-10, rtol=0):
        raise ValueError("Donor-first and gene-first estimates disagree")
    return dict(cases=int(cases.sum()), controls=int((~cases).sum()), strata=len(set(strata)),
                case_minus_control_pp=effect, assignments=len(reference), positive_tail=extreme,
                positive_p=extreme / len(reference), two_sided_tail=both, two_sided_p=both / len(reference),
                constant_shift_95_interval=[None if v is None else float(v) for v in interval],
                per_gene_native_pp=dict(zip([r["symbol"] for r in candidates], native.tolist(), strict=True)),
                per_band_oriented_pp=dict(zip(distinct, band_effects.tolist(), strict=True)),
                omit_one_band_pp={b: float(np.delete(band_effects, j).mean())
                                  for j, b in enumerate(distinct) if len(distinct) > 1})


def analyze(book, data, swiss_files=None, rscript="Rscript"):
    hgnc = data / "gse266265" / HGNC_FILE
    count, selected = selected_genes(book, hgnc)
    fin = load_finland(data / "gse141549")
    rg, rx, rc, _ = load_rome(data / "gse25628")
    scan_cases, dates, _ = scan_design(data / "gse25628")
    if not np.array_equal(scan_cases, rc):
        raise ValueError("Rome series and scanner groups disagree")
    mixed = np.isin(dates, np.intersect1d(dates[rc], dates[~rc]))
    rome = rg, rx[mixed], rc[mixed], dates[mixed]
    mg, mx, mc, _ = load_mumbai(data / "gse153739", hgnc)
    mumbai = mg, mx, mc, np.repeat("proliferative", len(mc))
    membership = {"Finland": set(fin[0]), "Rome": set(rg), "Mumbai": set(mg)}
    for row in selected:
        row["missing_from"] = [s for s, present in membership.items() if row["entrez"] not in present]
    common = [r for r in selected if not r["missing_from"]]
    cohorts = {study: score_cohort(common, *values)
               for study, values in [("Finland", fin), ("Mumbai", mumbai), ("Rome", rome)]}
    adjusted = holm(np.array([cohorts[s]["positive_p"] for s in ("Finland", "Mumbai")]))
    for study, value in zip(("Finland", "Mumbai"), adjusted, strict=True):
        cohorts[study]["holm_positive_p"] = float(value)
    fg, fx, fc, fs = fin
    prolif = np.char.endswith(fs, "-proliferative")
    cohorts["Finland_proliferative"] = score_cohort(common, fg, fx[prolif], fc[prolif], fs[prolif])
    if sha256_file(DISCOVERY) != DISCOVERY_SHA:
        raise ValueError("Finnish reference table differs")
    with DISCOVERY.open(newline="") as stream:
        reference = {int(r["id"]):float(r["effect_pp"]) for r in csv.DictReader(stream, delimiter="\t")}
    for row in common:
        if abs(cohorts["Finland"]["per_gene_native_pp"][row["symbol"]] - reference[row["entrez"]]) > 1e-10:
            raise ValueError("Prior Finnish gene contrast failed to reproduce")
    # Secondary comparisons were specified after seeing the primary coverage and results.
    complete = {study: score_cohort(selected, *design)
                for study, design in (("Finland", fin), ("Mumbai", mumbai))}
    adjust_complete = holm(np.array([value["positive_p"] for value in complete.values()]))
    for result, adjusted_p in zip(complete.values(), adjust_complete, strict=True):
        result["holm_positive_p"] = float(adjusted_p)
    secondary = {"complete_panel_bulk": complete}
    if swiss_files:
        profiles = load_glandular_profiles(*swiss_files, hgnc, rscript)
        available = set(profiles[0])
        swiss = {}
        for name, panel in (("common_panel", common), ("complete_panel", selected)):
            missing = [r["symbol"] for r in panel if r["entrez"] not in available]
            swiss[name] = {"missing_genes": missing}
            if not missing:
                swiss[name]["result"] = score_cohort(panel, *profiles)
        secondary["Swiss_glandular"] = swiss
    return dict(plan=PLAN, source_sha256=BOOK_SHA, models=count, selection_threshold=.05/count,
                selected=selected, measured_bands=dict(Counter(r["band"] for r in common)),
                primary_replication=bool(all(v <= .05 for v in adjusted)
                    and all(cohorts[k]["case_minus_control_pp"] > 0 for k in ("Finland", "Mumbai"))),
                cohorts=cohorts, secondary=secondary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("data", type=Path)
    parser.add_argument("--swiss", nargs=2, type=Path, metavar=("RDS", "CLINICAL_WORKBOOK"))
    parser.add_argument("--rscript", default="Rscript")
    args = parser.parse_args()
    print(json.dumps(analyze(args.workbook, args.data, args.swiss, args.rscript), indent=2, allow_nan=False))
