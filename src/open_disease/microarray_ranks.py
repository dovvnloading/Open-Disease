"""Within-array percentile ranks and equal-weight collapse of probe measurements."""
from __future__ import annotations

import numpy as np
from scipy.stats import rankdata


def gene_percentiles(signals, entrez):
    """Return sorted positive Entrez IDs and a donor-by-gene percentile matrix.

    All supplied feature columns, including those with an Entrez ID of zero,
    form the rank reference. Give tied signals their average rank. An ID of
    zero contributes only to the reference. Average all probes of each gene
    within each donor so there is exactly one measurement per gene.
    """
    values = np.asarray(signals, dtype=np.float64)
    ids = np.asarray(entrez)
    if (values.ndim != 2 or 0 in values.shape or ids.ndim != 1
            or len(ids) != values.shape[1] or ids.dtype.kind not in "iu"):
        raise ValueError("Expected a nonempty donor-by-feature matrix and aligned integer Entrez IDs")
    if not np.isfinite(values).all() or (ids < 0).any():
        raise ValueError("Expected finite feature signals and nonnegative Entrez IDs")
    mapped = np.flatnonzero(ids)
    if not mapped.size:
        raise ValueError("No eligible gene mapping")
    order = mapped[np.argsort(ids[mapped], kind="stable")]
    genes, starts, counts = np.unique(ids[order], return_index=True, return_counts=True)

    # The breakpoints aggregate contiguous probes of each gene, after each array
    # has already been ranked against its full physical feature universe.
    ranks = rankdata(values, axis=1, method="average")
    doubled_midranks = 2 * ranks - 1
    totals = np.add.reduceat(doubled_midranks[:, order], starts, axis=1)
    return genes, totals * (50.0 / (values.shape[1] * counts))
