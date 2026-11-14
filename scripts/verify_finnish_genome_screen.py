"""Independently verify the Finnish discovery screen using per-stratum combination sums.

The verification uses the same clinical gates and source annotation as the
primary screen. It recomputes ranks, group scores, exhaustive null standard
deviations, and multiplicity tails; it never opens the held-out Rome matrix.
Only cohort-level feature statistics or fingerprints are written.
"""
from __future__ import annotations

import argparse
import csv
from itertools import combinations
import json
from pathlib import Path

import numpy as np

from analyze_endometdb_wg6_fosb import N_FEATURES, verified_wg6_source
from genome_endometdb import load_wg6_annotation, selected_signals
from genome_rome import rome_annotation_ids
from pathway_atlas_reactome import load_pathways


def sorted_gene_scores(signals, entrez):
    """Independently locate both ends of each tied feature in the sorted array."""
    genes = np.unique(entrez[entrez > 0])
    mapped = entrez > 0
    indexes = np.searchsorted(genes, entrez[mapped])
    counts = np.bincount(indexes, minlength=len(genes))
    scores = np.zeros((len(signals), len(genes)))
    for donor, x in enumerate(signals):
        ordered = np.sort(x)
        intensity = x[mapped]
        # At positions [left,right), 100*(one-based midrank-.5)/N
        # equals 50*(left+right)/N, including all rank-reference features.
        ties = np.searchsorted(ordered, intensity, "left") + np.searchsorted(ordered, intensity, "right")
        gene_sums = np.bincount(indexes, weights=ties, minlength=len(genes))
        scores[donor] = gene_sums * 50 / (counts * len(ordered))
    return genes, scores


def recompute(matrix, cases, strata, block=48):
    """Build the full sharp null from Cartesian sums; obtain variance from all its values."""
    slices = [np.flatnonzero(strata == name) for name in sorted(set(strata))]
    choices = [np.array(list(combinations(range(len(ix)), int(cases[ix].sum())))) for ix in slices]
    nperm = int(np.prod([len(combinations_) for combinations_ in choices]))
    nfeatures = matrix.shape[1]
    raw = np.empty(nfeatures, dtype=int)
    observed = np.zeros(nfeatures)
    effect = np.zeros(nfeatures)
    maxima = np.zeros(nperm)
    for start in range(0, nfeatures, block):
        section = slice(start, min(start + block, nfeatures))
        selected = matrix[:, section]
        null_sum = np.zeros((1, selected.shape[1]))
        observed_sum = np.zeros(selected.shape[1])
        varies = np.zeros(selected.shape[1], dtype=bool)
        for ix, choose in zip(slices, choices):
            values = selected[ix, :]
            case_flags = cases[ix]
            ncase = int(case_flags.sum())
            ncontrol = len(ix) - ncase
            if min(ncase, ncontrol) < 1:
                raise ValueError("Stratum lacks a case or control")
            varies |= values.max(axis=0) != values.min(axis=0)
            case_sums = values[choose, :].sum(axis=1)
            group_delta = case_sums / ncase - (values.sum(axis=0) - case_sums) / ncontrol
            observed_sum += values[case_flags, :].mean(axis=0) - values[~case_flags, :].mean(axis=0)
            null_sum = (null_sum[:, None, :] + group_delta[None, :, :]).reshape((-1, selected.shape[1]))
        sd = np.std(null_sum, axis=0, ddof=0)
        if (sd[varies] <= 0).any():
            raise ValueError("Null has zero dispersion for a nonconstant feature")
        z = np.divide(null_sum, sd, out=np.zeros_like(null_sum), where=varies)
        zobs = np.divide(observed_sum, sd, out=np.zeros_like(observed_sum), where=varies)
        raw[section] = np.count_nonzero(np.abs(z) >= np.abs(zobs) - 1e-10, axis=0)
        maxima = np.maximum(maxima, np.max(np.abs(z), axis=1))
        observed[section] = zobs
        effect[section] = observed_sum / len(slices)
    adjusted = nperm - np.searchsorted(np.sort(maxima), np.abs(observed) - 1e-10)
    return effect, observed, raw, adjusted, nperm


def write_family(path, ids, outcomes, memberships=None):
    effects, z, raw, adjusted, permutations = outcomes
    # Reconstruct candidate order independently of table order.
    order = np.lexsort((np.array([int(str(i).split("-")[-1]) for i in ids]), -np.round(np.abs(z), 10)))
    selected = []
    for i in order:
        if z[i] == 0:
            continue
        fid = ids[i]
        if memberships:
            terms = set(memberships[fid][1])
            if any(2 * len(terms & set(memberships[j][1])) >= len(terms | set(memberships[j][1]))
                   for j in selected):
                continue
        selected.append(fid)
        if len(selected) == 5:
            break
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(["id", "effect", "z", "raw_extreme", "maxT_extreme", "permutations"])
        for i, identifier in enumerate(ids):
            writer.writerow([identifier, effects[i], z[i], raw[i], adjusted[i], permutations])
    return list(map(str, selected))


def compare_family(independent_path, primary_path):
    """Fail on a changed candidate universe, tail count, or materially different statistic."""
    with primary_path.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    primary = {row["id"]: row for row in rows}
    with independent_path.open(newline="") as fh:
        verified = list(csv.DictReader(fh, delimiter="\t"))
    if (len(primary) != len(rows) or len(primary) != len(verified)
            or primary.keys() != {row["id"] for row in verified}):
        raise ValueError("Discovery families disagree")
    for row in verified:
        published = primary[row["id"]]
        for field in ("raw_extreme", "maxT_extreme"):
            if int(row[field]) != int(published[field]):
                raise ValueError("Discovery tail counts disagree")
        if (abs(float(row["effect"]) - float(published["effect_pp"])) > 1e-9
                or abs(float(row["z"]) - float(published["statistic"])) > 1e-9):
            raise ValueError("Discovery estimate or standardized statistic disagrees")
    return len(verified)


def run(finland, rome, pathways, output, compare=None):
    output.mkdir(parents=True, exist_ok=True)
    annotation = load_wg6_annotation(finland)
    design, reference, source = verified_wg6_source(finland)
    if reference != annotation.keys():
        raise ValueError("Finnish probe census mismatch")
    physical, probe_ids, case, strata = selected_signals(source, design, annotation)
    genes, scores = sorted_gene_scores(physical, probe_ids)
    # Historical positive sentinel, never included in the screened families.
    fosb_pos = np.flatnonzero(genes == 2354)
    if len(fosb_pos) != 1 or np.count_nonzero(probe_ids == 2354) != 1 or physical.shape[1] != N_FEATURES:
        raise ValueError("Historical FOSB sentinel mismatch")
    sentinel = recompute(scores[:, fosb_pos], case, strata)
    if (sentinel[2][0], sentinel[-1]) != (1925, 12096) or abs(sentinel[0][0] - 10.891044576184512) > 1e-10:
        raise ValueError("Failed historical Finnish FOSB control")
    common = np.array(sorted((set(genes) & set(rome_annotation_ids(rome))) - {2354}))
    ix = np.searchsorted(genes, common)
    outcome = recompute(scores[:, ix], case, strata)
    nominated_genes = write_family(output / "genes.tsv", common, outcome)

    sets, _ = load_pathways(pathways, set(common))
    terms = sorted(sets, key=lambda name: int(name.split("-")[-1]))
    donor_means = np.stack([scores[:, np.searchsorted(genes, sets[term][1])].mean(axis=1) for term in terms], axis=1)
    outcome = recompute(donor_means, case, strata)
    nominated_terms = write_family(output / "pathways.tsv", terms, outcome, sets)
    report = {"gene_count": len(common), "pathway_count": len(terms),
              "sentinel_delta": float(sentinel[0][0]), "sentinel_extreme": int(sentinel[2][0]),
              "assignment_count": int(sentinel[-1]), "genes": nominated_genes, "pathways": nominated_terms}
    if compare is not None:
        report["verified_rows"] = {
            family: compare_family(output / f"{family}s.tsv", compare / f"discovery_{family}.tsv")
            for family in ("gene", "pathway")
        }
        lock = json.loads((compare / "nominations.json").read_text())
        independent = []
        for family, nominees in (("gene", nominated_genes), ("pathway", nominated_terms)):
            with (output / f"{family}s.tsv").open(newline="") as fh:
                feature = {row["id"]: float(row["effect"]) for row in csv.DictReader(fh, delimiter="\t")}
            independent.extend((i, "higher" if feature[i] > 0 else "lower") for i in nominees)
        if [(r["id"], r["direction"]) for r in lock["candidates"]] != independent:
            raise ValueError("Frozen candidate order or direction differs between implementations")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--finland", type=Path, required=True)
    parser.add_argument("--rome", type=Path, required=True, help="Platform annotation only")
    parser.add_argument("--pathways", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True, help="Directory for aggregate verification tables")
    parser.add_argument("--compare", type=Path, help="Check all rows and ordered nominations against primary outputs")
    args = parser.parse_args()
    run(args.finland, args.rome, args.pathways, args.output, args.compare)
