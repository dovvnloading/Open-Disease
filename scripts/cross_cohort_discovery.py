"""Screen Finnish genome-wide relative expression, then test a locked list in Rome.

All array-level matrices stay in memory. Only cross-donor aggregate contrasts,
permutation counts and feature identifiers are written to the output directory.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
from itertools import combinations, product
import json
from math import prod
from pathlib import Path

import numpy as np

from genome_endometdb import finland_annotation_ids, load_finland
from genome_rome import load_rome, rome_annotation_ids
from pathway_atlas_reactome import load_pathways

TOL = 1e-10
PLAN = "94dc9a2965cea7640d37be11286c9df962650f12"
LOCK_HASH = "ba140f2bb026b0416e1c705ab28f616f430cc4295f94d174feaf5d15947c1e31"
GENE_HASH = "eb7a4192c14d0ab20494afdd16c341a0cf5ed1e1ebec5f20818ae5cea82e4695"
PATHWAY_HASH = "96143a21a5e0b313b74b1a7b37d431c486b0c06d21191fa88e4cf2a477e3e494"
DISCOVERY_FIELDS = ("family", "id", "members", "effect_pp", "statistic", "raw_extreme",
                    "p_raw", "maxT_extreme", "p_maxT")


def digest(value):
    value = json.dumps(value, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(value.encode()).hexdigest()


def universe(data):
    gene_ids = np.intersect1d(finland_annotation_ids(data / "gse141549"),
                              rome_annotation_ids(data / "gse25628"))
    gene_ids = gene_ids[gene_ids != 2354]
    pathways, qc = load_pathways(data / "pathway-atlas", gene_ids)
    serialized = [[key, name, list(members)] for key, (name, members) in pathways.items()]
    if (len(gene_ids) != 12044 or digest(gene_ids.tolist()) != GENE_HASH
            or len(pathways) != 1049 or digest(serialized) != PATHWAY_HASH):
        raise ValueError("The common genes or projected pathway universe changed")
    pins = {"genes": digest(gene_ids.tolist()), "pathways": digest(serialized)}
    return gene_ids, pathways, pins, qc


def align(genes, scores, requested):
    positions = np.searchsorted(genes, requested)
    if np.any(positions >= len(genes)) or not np.array_equal(genes[positions], requested):
        raise ValueError("A locked gene is absent from the measured matrix")
    return scores[:, positions]


def assignment_design(case, strata):
    """One contrast-weight row per possible allocation, including the observed."""
    case = np.asarray(case)
    strata = np.asarray(strata)
    if case.ndim != 1 or case.dtype != np.bool_ or strata.shape != case.shape:
        raise ValueError("Misaligned clinical design")
    cells = [np.flatnonzero(strata == label) for label in np.unique(strata)]
    observed = np.zeros(len(case))
    options, factors = [], []
    for ids in cells:
        ncase = int(case[ids].sum())
        ncontrol = len(ids) - ncase
        if min(ncase, ncontrol) < 1:
            raise ValueError("Stratum lacks cases or controls")
        weights = []
        for combo in combinations(range(len(ids)), ncase):
            w = np.full(len(ids), -1 / ncontrol)
            w[list(combo)] = 1 / ncase
            weights.append(w / len(cells))
        options.append(np.asarray(weights))
        observed[ids] = np.where(case[ids], 1 / ncase, -1 / ncontrol) / len(cells)
        factors.append((ids, (1 / ncase + 1 / ncontrol) / len(cells)**2))
    n = prod(len(w) for w in options)
    choices = np.asarray(list(product(*(range(len(w)) for w in options))), dtype=int)
    weights = np.empty((n, len(case)))
    for i, (ids, _) in enumerate(factors):
        weights[:, ids] = options[i][choices[:, i]]
    if np.count_nonzero(np.all(weights == observed, axis=1)) != 1:
        raise ValueError("Observed assignment missing or repeated")
    return observed, weights, factors


def test_family(scores, design, directions=None, chunk_size=128):
    """Apply the synchronized exact permutation distribution to one feature family."""
    observed, weights, factors = design
    x = np.asarray(scores, dtype=np.float64)
    if x.ndim != 2 or x.shape[0] != len(observed) or not np.isfinite(x).all() or not x.shape[1]:
        raise ValueError("Invalid feature matrix")
    # Removing the stratum's fixed mean leaves every contrast intact and makes
    # exactly constant columns zero before standardization.
    centered = np.empty_like(x)
    variance = np.zeros(x.shape[1])
    for ids, factor in factors:
        values = x[ids] - x[ids[0]]
        centered[ids] = values - values.mean(axis=0)
        variance += values.var(axis=0, ddof=1) * factor
    scale = np.sqrt(variance)
    variable = scale > 0
    effect = observed @ centered
    norm = np.divide(centered, scale, out=np.zeros_like(centered), where=variable)
    z = observed @ norm
    one_sided = directions is not None
    if one_sided:
        directions = np.asarray(directions)
        if directions.shape != z.shape or not np.isin(directions, (-1, 1)).all():
            raise ValueError("Expected one locked sign per feature")
    null_max = np.zeros(weights.shape[0])
    raw = np.empty(len(z), dtype=np.int64)
    for start in range(0, len(z), chunk_size):
        sl = slice(start, min(start + chunk_size, len(z)))
        perm = weights @ norm[:, sl]
        if one_sided:
            perm *= directions[sl]
            target = z[sl] * directions[sl]
        else:
            np.abs(perm, out=perm)
            np.maximum(null_max, perm.max(axis=1), out=null_max)
            target = np.abs(z[sl])
        raw[sl] = np.count_nonzero(perm >= target - TOL, axis=0)
    if (raw < 1).any():
        raise ValueError("Permutation tail omitted the observed assignment")
    # searchsorted includes ties in the right tail, identically to the raw test.
    adj = None
    if not one_sided:
        adj = len(weights) - np.searchsorted(np.sort(null_max), np.abs(z) - TOL)
    return effect, z, raw, adj, variable


def holm(ps):
    ps = np.asarray(ps, dtype=float)
    order = np.argsort(ps, kind="stable")
    adjusted = np.empty_like(ps)
    adjusted[order] = np.minimum(1, np.maximum.accumulate(ps[order] * np.arange(len(ps), 0, -1)))
    return adjusted


def nominate(rows, member_map):
    chosen = []
    for row in sorted(rows, key=lambda r: (-round(abs(r["statistic"]), 10),
                                            int(r["id"].split("-")[-1]))):
        if row["statistic"] == 0:
            continue
        if member_map is not None:
            group = set(member_map[row["id"]][1])
            if any(2 * len(group & other) >= len(group | other) for _, other in chosen):
                continue
        else:
            group = set()
        chosen.append((row, group))
        if len(chosen) == 5:
            break
    return [dict(row, direction="higher" if row["effect_pp"] > 0 else "lower")
            for row, _ in chosen]


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, delimiter="\t", lineterminator="\n", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def discovery(data, out):
    shared, paths, pins, qc = universe(data)
    genes, scores, case, strata = load_finland(data / "gse141549")
    design = assignment_design(case, strata)
    n = len(design[1])
    if n != 12096 or len(case) != 25 or int(case.sum()) != 12:
        raise ValueError("Unexpected Finnish design")
    g = align(genes, scores, shared)
    path_scores = np.column_stack([align(shared, g, members).mean(axis=1)
                                   for _, members in paths.values()])
    candidates, summary = [], {}
    for family, ids, x in (("gene", list(map(str, shared)), g), ("pathway", list(paths), path_scores)):
        effect, z, raw, adj, variable = test_family(x, design)
        rows = [dict(family=family, id=key, members=1 if family == "gene" else len(paths[key][1]),
                     effect_pp=float(effect[i]), statistic=float(z[i]), raw_extreme=int(raw[i]),
                     p_raw=float(raw[i] / n), maxT_extreme=int(adj[i]), p_maxT=float(adj[i] / n))
                for i, key in enumerate(ids)]
        write_tsv(out / f"discovery_{family}.tsv", DISCOVERY_FIELDS, rows)
        candidates += nominate(rows, paths if family == "pathway" else None)
        summary[family] = dict(tested=len(rows), variable=int(variable.sum()),
                               threshold_hits=int((adj / n <= 0.025).sum()))
    for row in candidates:
        if row["family"] == "pathway":
            row["name"], row["entrez_ids"] = paths[row["id"]]
    lock = dict(protocol=PLAN, universe=pins, assignment_count=n, summary=summary,
                pathway_qc=qc, candidates=candidates)
    (out / "nominations.json").write_text(json.dumps(lock, ensure_ascii=False, indent=2) + "\n")
    return lock


def validation(data, lock_path, out):
    payload = lock_path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != LOCK_HASH:
        raise ValueError("Validation requires the public nomination snapshot from commit f2a3755")
    lock = json.loads(payload)
    shared, paths, pins, _ = universe(data)
    if lock["protocol"] != PLAN or lock["universe"] != pins:
        raise ValueError("Nomination or universe mismatch")
    candidates = lock["candidates"]
    if Counter(r["family"] for r in candidates) != {"gene": 5, "pathway": 5}:
        raise ValueError("Unexpected nomination family")
    if len({(r["family"], r["id"]) for r in candidates}) != len(candidates):
        raise ValueError("Duplicate nomination")
    selected = []
    for r in candidates:
        if r["direction"] not in {"higher", "lower"} or (r["effect_pp"] > 0) != (r["direction"] == "higher"):
            raise ValueError("Locked direction mismatch")
        members = [int(r["id"])] if r["family"] == "gene" else list(paths[r["id"]][1])
        if r.get("entrez_ids", members) != members or not np.isin(members, shared).all():
            raise ValueError("Locked gene membership differs")
        selected.append(members)
    genes, scores, case, strata = load_rome(data / "gse25628")
    family = np.column_stack([align(genes, scores, members).mean(axis=1) for members in selected])
    signs = np.array([1 if r["direction"] == "higher" else -1 for r in candidates])
    design = assignment_design(case, strata)
    n = len(design[1])
    if n != 3003 or len(case) != 14 or int(case.sum()) != 8:
        raise ValueError("Unexpected Rome design")
    effect, _, raw, _, _ = test_family(family, design, signs)
    adj = holm(raw / n)
    deleted = []
    for i in range(len(case)):
        keep = np.arange(len(case)) != i
        deleted.append(family[keep & case].mean(axis=0) - family[keep & ~case].mean(axis=0))
    lo, hi = np.min(deleted, axis=0), np.max(deleted, axis=0)
    result = []
    for i, r in enumerate(candidates):
        passed = bool(effect[i] * signs[i] > 0 and adj[i] <= .05)
        result.append(dict(family=r["family"], id=r["id"], direction=r["direction"],
                           finnish_effect_pp=r["effect_pp"], finnish_p_maxT=r["p_maxT"],
                           rome_effect_pp=float(effect[i]), rome_extreme=int(raw[i]),
                           rome_p=float(raw[i]/n), rome_p_holm=float(adj[i]),
                           rome_loo_low=float(lo[i]), rome_loo_high=float(hi[i]),
                           validation_pass=passed, both_pass=bool(passed and r["p_maxT"] <= .025)))
    write_tsv(out / "validation.tsv", result[0].keys(), result)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=["discover", "validate"])
    parser.add_argument("--data", type=Path, default=Path("data/private"))
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--lock", type=Path)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.stage == "discover":
        record = discovery(args.data, args.out)
        print(json.dumps(record["summary"]))
    elif args.lock is None:
        parser.error("validate requires --lock")
    else:
        rows = validation(args.data, args.lock, args.out)
        print(json.dumps({"tested": len(rows), "validated": sum(r["validation_pass"] for r in rows),
                          "supported": sum(r["both_pass"] for r in rows)}))
