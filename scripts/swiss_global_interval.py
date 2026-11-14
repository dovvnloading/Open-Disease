"""Post-outcome inversion of the Swiss global score's exact conditional test.

Run with the same pinned RDS, workbook and HGNC directory as the main study.
Only aggregate inference is written to stdout.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
from itertools import combinations, product
import json
import os
from pathlib import Path
import subprocess
import tempfile

from audit_gse266265_epithelial import verify, write_sample_map
from swiss_epithelial_mapping import load_projection, write_projection
from verify_swiss_glandular import TOLERANCE, load_pathways, reconstruct


def allocation_contrasts(outcomes, labels, batches):
    """For each assignment return its equal-batch contrasts of y and observed d."""
    if len(outcomes) == 0 or not len(outcomes) == len(labels) == len(batches):
        raise ValueError("Donor scores, labels and capture batches must align")
    outcomes = [Fraction(value) for value in outcomes]
    if not all(label in (0, 1) for label in labels):
        raise ValueError("Disease labels must be binary")
    observed = Fraction(0)
    conditions = []
    for batch in sorted(set(batches)):
        rows = [i for i, b in enumerate(batches) if b == batch]
        cases = [i for i in rows if labels[i]]
        n, k = len(rows), len(cases)
        if not 0 < k < n:
            raise ValueError("Each capture batch must include both groups")
        total = sum(outcomes[i] for i in rows)
        weight = Fraction(1, k) + Fraction(1, n-k)
        observed += sum(outcomes[i] for i in cases) * weight - total / (n-k)
        conditions.append([(sum(outcomes[i] for i in selected) * weight - total / (n-k),
                            sum(labels[i] for i in selected) * weight - Fraction(k, n-k))
                           for selected in combinations(rows, k)])
    n_batches = len(conditions)
    contrasts = [(sum(a for a, _ in assignment)/n_batches, sum(b for _, b in assignment)/n_batches)
                 for assignment in product(*conditions)]
    return observed / n_batches, contrasts


def shift_limits(observed, contrasts, alpha=Fraction(1, 20), tolerance=TOLERANCE):
    """Invert |a−bθ|+ε ≥ |observed−θ| exactly; equality belongs to the set."""
    required = int(alpha * len(contrasts)) + 1  # p > alpha
    def radius(center_value, slope):
        # At distance t from the observed contrast, inclusion means
        # |center_value+slope*t|+tolerance >= t. Either sign can qualify.
        width = Fraction(0)
        for sign in (-1, 1):
            numerator = sign * center_value + tolerance
            denominator = 1 - sign * slope
            if denominator == 0:
                if numerator >= 0:
                    return None  # All distances qualify on this side.
            else:
                width = max(width, numerator / denominator)
        return width

    left, right = [], []
    for a, b in contrasts:
        if not -1 <= b <= 1:
            raise ValueError("Label contrast outside [-1,1]")
        c = a - b * observed
        left.append(radius(c, b))
        right.append(radius(c, -b))
    # Each included interval contains observed. The set where ≥required
    # intervals overlap is connected, even when several boundary points tie.
    def boundary(widths, direction):
        finite = sorted((x for x in widths if x is not None), reverse=True)
        index = required - (len(widths) - len(finite)) - 1
        return None if index < 0 else observed + direction * finite[index]
    return boundary(left, -1), boundary(right, 1), required


def direct_tail(outcomes, labels, batches, shift, tolerance=TOLERANCE):
    adjusted = [y - shift * d for y, d in zip(outcomes, labels, strict=True)]
    observed, contrasts = allocation_contrasts(adjusted, labels, batches)
    return sum(abs(value) + tolerance >= abs(observed) for value, _ in contrasts)


def summarize(outcomes, labels, batches):
    observed, contrasts = allocation_contrasts(outcomes, labels, batches)
    lower, upper, required = shift_limits(observed, contrasts)
    if lower is None or upper is None:
        raise ValueError("Unexpected unbounded interval for this source")
    delta = Fraction(1, 10**8)
    checks = {"zero": 0, "lower_minus_1e-8": lower-delta, "lower": lower,
              "lower_plus_1e-8": lower+delta, "upper_minus_1e-8": upper-delta,
              "upper": upper, "upper_plus_1e-8": upper+delta}
    tails = {key: direct_tail(outcomes, labels, batches, shift) for key, shift in checks.items()}
    if (len(contrasts), tails['zero']) != (144, 24):
        raise ValueError("Null test disagrees with the frozen study")
    if not (tails['lower_minus_1e-8'] < required <= min(tails['lower'], tails['lower_plus_1e-8'])
            and tails['upper_plus_1e-8'] < required <= min(tails['upper'], tails['upper_minus_1e-8'])):
        raise ValueError("Failed direct enumeration at the two analytic boundaries")
    return {"endpoint": "GLOBAL", "status": "post-outcome descriptive",
            "assumption": "Constant additive shift on donor global score and null exchangeability within capture batch",
            "units": "percentile points", "donors": len(labels), "batches": len(set(batches)),
            "assignments": len(contrasts), "alpha": 0.05, "minimum_extreme_for_inclusion": required,
            "absolute_tail_tolerance_pp": float(TOLERANCE), "difference": float(observed),
            "closed_interval": [float(lower), float(upper)], "exact_bounds": [str(lower), str(upper)],
            "direct_extreme_counts": tails}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('rds', type=Path)
    parser.add_argument('supplement', type=Path)
    parser.add_argument('annotation_directory', type=Path)
    parser.add_argument('--rscript', default='Rscript')
    args = parser.parse_args()
    verify(args.rds)
    env = dict(os.environ, R_DEFAULT_PACKAGES='NULL', OPENBLAS_NUM_THREADS='1')
    r_dir = Path(__file__).parent / 'r'
    with tempfile.TemporaryDirectory(dir=args.rds.parent) as folder:
        temp = Path(folder)
        mapping, features = temp / 'mapping.tsv', temp / 'rows.txt'
        clinical, stream = temp / 'study.tsv', temp / 'columns.bin'
        write_sample_map(args.supplement, clinical)
        subprocess.run([args.rscript, r_dir/'extract_swiss_feature_labels.R', args.rds, features], env=env, check=True)
        rows, _ = load_projection(args.annotation_directory, labels_path=features)
        write_projection(mapping, rows)
        pathways = load_pathways(mapping)
        subprocess.run([args.rscript, r_dir/'export_swiss_verification_columns.R',
                        args.rds, clinical, mapping, stream], env=env, check=True)
        labels, batches, _, primary, _ = reconstruct(stream, pathways)
    result = summarize([row[0] for row in primary], [int(d) for d in labels], batches.tolist())
    print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
