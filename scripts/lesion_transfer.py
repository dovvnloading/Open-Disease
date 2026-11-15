"""Evaluate a frozen ovarian-to-extraovarian tissue contrast using paired arrays."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from itertools import product
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, spearmanr

from genetic_anchor_tissue import (
    GENES, SERIES, checked, holm, paired_summary, read_annotation, read_paired_arrays,
)


DISCOVERY, VALIDATION = "GSE7305", "GSE11691"


def feature_sets(directory):
    mappings = [read_annotation(checked(directory / f"{platform}.annot.gz"))
                for platform, _, _ in SERIES.values()]
    a, b = mappings
    excluded = {0, *GENES.values()}
    ids = sorted((set(a.values()) & set(b.values())) - excluded)
    shared_probes = {p for p in a.keys() & b.keys() if a[p] == b[p] and a[p] not in excluded}
    shared_ids = sorted({a[p] for p in shared_probes})
    return np.array(ids), shared_probes, np.array(shared_ids)


def gene_differences(arrays, ids, reference=None):
    probes, values, lesion, uterus, mapping = arrays
    if reference is not None:
        take = np.array([p in reference for p in probes])
        values = values[take]
        probes = [p for p, keep in zip(probes, take, strict=True) if keep]
    ranks = 100 * (rankdata(values, axis=0) - 0.5) / len(probes)
    gene_index = {g: i for i, g in enumerate(ids)}
    groups = np.array([gene_index.get(mapping[p], -1) for p in probes])
    keep = groups >= 0
    counts = np.bincount(groups[keep], minlength=len(ids))
    if not np.all(counts):
        raise ValueError("A supposedly measured common gene has no probes")
    deltas = np.zeros((len(ids), len(lesion)))
    np.add.at(deltas, groups[keep], (ranks[:, lesion] - ranks[:, uterus])[keep])
    return deltas / counts[:, None]


def panel_indices(effect, ids):
    high = [i for i in np.lexsort((ids, -effect)) if effect[i] > 0][:2]
    low = [i for i in np.lexsort((ids, effect)) if effect[i] < 0][:2]
    if len(high) != 2 or len(low) != 2:
        raise ValueError("Need at least two nonzero genes in each direction")
    chosen = high + low
    fifth = next(i for i in np.lexsort((ids, -abs(effect))) if i not in chosen)
    return chosen + [fifth]


def frozen_panel(ids, disc):
    effect = disc.mean(axis=1)
    return [{"entrez": int(ids[i]), "discovery_pp": float(effect[i]),
             "direction": int(np.sign(effect[i]))} for i in panel_indices(effect, ids)]


def projection(weights, validation):
    weight_sum = np.abs(weights).sum()
    if weight_sum == 0:
        raise ValueError("No discovery variation")
    return (weights / weight_sum) @ validation


def summarize(scores):
    summary = paired_summary(scores)
    # The column weights are frozen before examining validation tissues.
    return {**summary, "one_sided_p": summary["upper"] / summary["assignments"]}


def gene_tails(diffs):
    n = diffs.shape[1]
    signs = np.array(list(product((-1., 1.), repeat=n)))
    upper, lower = [], []
    for start in range(0, len(diffs), 256):
        block = diffs[start:start + 256]
        observed = block.sum(axis=1)[:, None]
        null = block @ signs.T
        upper.extend((null >= observed - 1e-9).mean(axis=1))
        lower.extend((null <= observed + 1e-9).mean(axis=1))
    return np.array(upper), np.array(lower)


def gene_results(ids, discovery, validation):
    disc_up, disc_down = gene_tails(discovery)
    val_up, val_down = gene_tails(validation)
    conjunction = np.minimum(1, 2 * np.minimum(np.maximum(disc_up, val_up),
                                               np.maximum(disc_down, val_down)))
    corrected = np.asarray(holm(conjunction))
    return [(int(g), float(a), float(b), float(p), float(q))
            for g, a, b, p, q in zip(ids, discovery.mean(axis=1), validation.mean(axis=1),
                                     conjunction, corrected, strict=True)]


def transfer(ids, discovery, validation, panel):
    effect = discovery.mean(axis=1)
    test_effect = validation.mean(axis=1)
    selected = {r["entrez"]: r["direction"] for r in panel}
    present = np.flatnonzero(np.isin(ids, list(selected)))
    directions = np.array([selected[k] for k in ids[present]])
    if not len(directions):
        raise ValueError("No projected panel member")
    panel_result = summarize(projection(directions, validation[present]))
    panel_result["genes"] = len(directions)
    singles = []
    for i, direction in zip(present, directions, strict=True):
        summary = summarize(direction * validation[i])
        singles.append(dict(entrez=int(ids[i]), direction=int(direction),
                            validation_pp=float(test_effect[i]), directional_p=summary["one_sided_p"],
                            directional_pairs=summary["higher"]))
    adj = holm([g["directional_p"] for g in singles] + [1.] * (len(panel) - len(singles)))
    for gene, q in zip(singles, adj):
        gene["directional_holm"] = float(q)
    missing = sorted(set(selected) - set(ids[present]))
    leave = [projection((discovery.sum(axis=1) - discovery[:, i]) / (discovery.shape[1] - 1),
                        validation).mean() for i in range(discovery.shape[1])]
    weights = effect.copy()
    weights[np.lexsort((ids, -abs(effect)))[:int(np.ceil(len(ids) / 100))]] = 0
    return {
        "genes": len(ids),
        "global_pp": summarize(projection(effect, validation)),
        "panel_pp": panel_result, "panel_genes": singles, "missing_panel_entrez": missing,
        "pearson": float(np.corrcoef(effect, test_effect)[0, 1]),
        "spearman": float(spearmanr(effect, test_effect).statistic),
        "discovery_omit_one_refit_effect_pp": [float(min(leave)), float(max(leave))],
        "top_one_percent_removed": summarize(projection(weights, validation)),
        "validation_gene_sign_only": summarize(projection(effect, np.sign(validation))),
        "direction_agreement": int(np.sum(effect * test_effect > 0)),
    }


def read_gene_names(directory):
    names = {}
    with gzip.open(directory / "GPL570.annot.gz", "rt", newline="") as f:
        rows = csv.reader(f, delimiter="\t")
        for row in rows:
            if row == ["!platform_table_begin"]:
                break
        header = next(rows)
        g, s = header.index("Gene ID"), header.index("Gene symbol")
        for row in rows:
            if row == ["!platform_table_end"]:
                break
            if row[g].isdigit() and row[s] and "///" not in row[s]:
                names.setdefault(int(row[g]), set()).add(row[s])
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["lock", "validate"])
    parser.add_argument("data", type=Path)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--genes", type=Path, help="Optional aggregate gene results TSV")
    args = parser.parse_args()
    ids, probes, common_ids = feature_sets(args.data)
    discovery_arrays = read_paired_arrays(args.data, DISCOVERY)
    discovery = gene_differences(discovery_arrays, ids)
    lock = frozen_panel(ids, discovery)
    if args.mode == "lock":
        if args.lock.exists():
            raise FileExistsError("Request a new path for the panel lock")
        args.lock.write_text(json.dumps(lock, indent=2) + "\n")
        return
    frozen = json.loads(args.lock.read_text())
    if len(frozen) != len(lock) or any(
        a["entrez"] != b["entrez"] or a["direction"] != b["direction"]
        or not np.isclose(a["discovery_pp"], b["discovery_pp"], rtol=0, atol=1e-9)
        for a, b in zip(lock, frozen, strict=True)
    ):
        raise ValueError("Discovery panel does not reproduce")
    validation_arrays = read_paired_arrays(args.data, VALIDATION)
    validation = gene_differences(validation_arrays, ids)
    shared_discovery = gene_differences(discovery_arrays, common_ids, probes)
    shared_validation = gene_differences(validation_arrays, common_ids, probes)
    rows = gene_results(ids, discovery, validation)
    if args.genes:
        with args.genes.open("w", newline="") as f:
            writer = csv.writer(f, delimiter="\t", lineterminator="\n")
            writer.writerow(["entrez", "discovery_effect_pp", "validation_effect_pp",
                             "conjunction_p", "holm_p"])
            writer.writerows(rows)
    names = read_gene_names(args.data)
    summary = {
        "discovery": DISCOVERY, "validation": VALIDATION,
        "excluded_known_anchor_entrez": sorted(GENES.values()),
        "primary": transfer(ids, discovery, validation, frozen),
        "shared_probe_ids": len(probes),
        "shared_probe_sensitivity": transfer(common_ids, shared_discovery, shared_validation, frozen),
        "lowest_gene_conjunction_p": min(r[3] for r in rows),
        "holm_genes_lt_0.05": sum(r[4] < 0.05 for r in rows),
        "panel_symbols": {str(r["entrez"]): sorted(names[r["entrez"]]) for r in frozen},
    }
    print(json.dumps(summary, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
