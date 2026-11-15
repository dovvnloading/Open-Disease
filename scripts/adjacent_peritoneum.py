"""Score the frozen lesion panel across three matched sites in GSE179640."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import re
from tempfile import TemporaryDirectory
from urllib.request import urlopen

import numpy as np
from scipy.stats import rankdata

from spatial_lesion_panel import PANEL_HASH

SOFT = "GSE179640_family.soft.gz"
LIST = "filelist.txt"
SITES = ("Ectopic", "Ectopic Adjacent", "Eutopic")
MAPPING = {
    55026: ("TMEM255A", "ENSG00000125355"),
    9077: ("DIRAS3", "ENSG00000162595"),
    2568: ("GABRP", "ENSG00000094755"),
    6299: ("SALL1", "ENSG00000103449"),
    3512: ("JCHAIN", "ENSG00000132465"),
}
METADATA = {
    SOFT: ("soft", "1eeaa9b3588effa0dff04ef625f4f3d5299bce6562701503630f7f57bb48b4f9"),
    LIST: ("suppl", "f63e02a493e7b11ac8df3e43da4272adc6716ae0476370a2511f6ef620cc8d0f"),
}


def acquire(directory, name, url, limit, fetch):
    path = directory / name
    if not path.is_file() and fetch:
        with urlopen(url, timeout=90) as stream:
            data = stream.read(limit + 1)
        if len(data) > limit:
            raise ValueError("Source object exceeds specified budget")
        directory.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory(dir=directory) as tmp:
            staged = Path(tmp) / name
            staged.write_bytes(data)
            staged.replace(path)
    return path


def design(directory, fetch):
    for name, (section, digest) in METADATA.items():
        url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/{section}/{name}"
        path = acquire(directory, name, url, 1_000_000, fetch)
        if sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError("Unexpected metadata")
    with (directory / LIST).open() as f:
        sizes = {r["Name"]: int(r["Size"]) for r in csv.DictReader(f, delimiter="\t")}
    blocks = gzip.open(directory / SOFT, "rt").read().split("^SAMPLE = ")[1:]
    if len(blocks) != 59:
        raise ValueError("Unexpected sample count")
    donors = defaultdict(list)
    counts = Counter()
    for block in blocks:
        fields = dict(
            x.split(": ", 1)
            for x in re.findall(r"^!Sample_characteristics_ch1 = (.*)$", block, re.M)
        )
        counts[fields["method"], fields["sample location"]] += 1
        site = fields["sample location"]
        if (
            fields["condition"] != "Endometriosis"
            or fields["method"] != "scRNA-seq"
            or site not in SITES
        ):
            continue
        if fields["tissue"] != ("Endometrium" if site == "Eutopic" else "Peritoneum"):
            continue
        title = re.search(r"^!Sample_title = (.*)$", block, re.M)[1]
        patient = re.fullmatch(
            r"Endometriosis Patient ([0-9]+) - (?:Ectopic(?: Adjacent.*)?|Eutopic)", title
        )
        if not patient:
            raise ValueError("Unrecognized patient key")
        (url,) = re.findall(
            r"^!Sample_supplementary_file_1 = (ftp://ftp.ncbi.nlm.nih.gov/geo/samples/.*\.h5)$",
            block,
            re.M,
        )
        name = url.rsplit("/", 1)[1]
        if not name.startswith(block.splitlines()[0] + "_"):
            raise ValueError("GEO accession inconsistent with H5 mapping")
        donors[patient[1]].append((site, url.replace("ftp:", "https:", 1), name, sizes[name]))
    matched = {
        p: records for p, records in donors.items() if {v[0] for v in records} >= set(SITES[:2])
    }
    if counts["RNA-seq", "Ectopic Adjacent"] != 0 or len(matched) != 6:
        raise ValueError("Unexpected pairing")
    observed = Counter(x[0] for records in matched.values() for x in records)
    if observed != dict(zip(SITES, [6, 7, 6], strict=True)):
        raise ValueError("Incomplete three-site design")
    return matched


def counts(path, identity, targets):
    import h5py

    with h5py.File(path, "r") as f:
        g = f["matrix"]
        labels = g["features/name"][:]
        genes = g["features/id"][:]
        if set(g["features/feature_type"][:]) != {b"Gene Expression"} or len(genes) != len(
            set(genes)
        ):
            raise ValueError("Unexpected assay features")
        digest = sha256(b"\n".join(genes) + b"\0" + b"\n".join(labels)).digest()
        if identity is not None and digest != identity:
            raise ValueError("Different feature order across samples")
        idx = []
        for symbol, ensembl in targets:
            ix = np.flatnonzero((labels == symbol.encode()) & (genes == ensembl.encode()))
            if len(ix) != 1:
                raise ValueError("Missing or ambiguous target")
            idx.append(ix[0])
        nrow, ncol = map(int, g["shape"][:])
        ix = g["indices"][:]
        values = g["data"][:]
        ptr = g["indptr"][:]
        if (
            nrow != len(genes)
            or len(ptr) != ncol + 1
            or ptr[0] != 0
            or ptr[-1] != len(values)
            or len(values) != len(ix)
            or np.any(np.diff(ptr.astype(np.int64)) < 0)
            or np.any(values < 0)
            or ix.max() >= nrow
        ):
            raise ValueError("Invalid sparse counts")
        totals = np.bincount(ix, weights=values, minlength=nrow).astype(np.int64)
        if totals.sum() != int(values.sum(dtype=np.uint64)):
            raise ValueError("Pseudobulk sum mismatch")
        return digest, idx, totals, ncol


def ranks(totals, idx):
    return (rankdata(totals, method="average")[idx] - 0.5) * 100 / len(totals)


def paired_shift_interval(scores):
    """Invert the two-sided test; each sign allocation contributes one closed interval."""
    scores = np.asarray(scores, dtype=float)
    n = len(scores)
    mean = scores.mean()
    left, right = [], []
    for signs in product([-1, 1], repeat=n):
        s = np.array(signs)
        if np.all(s == s[0]):
            continue  # Identical or reversed observed assignment always enters both tails.
        subset = s == 1
        # Equality with the absolute observed mean occurs at these two subset means.
        # The permuted tail exceeds the observed tail precisely between them.
        endpoints = [scores[subset].mean(), scores[~subset].mean()]
        left.append(min(endpoints))
        right.append(max(endpoints))
    needed = int(np.floor(0.05 * (2**n))) + 1 - 2
    if needed <= 0:
        return [None, None]
    bounds = [sorted(left)[needed - 1], sorted(right, reverse=True)[needed - 1]]
    for endpoint in bounds:
        centered = scores - endpoint
        perm = np.array(list(product([-1, 1], repeat=n))) @ centered / n
        if np.mean(abs(perm) >= abs(mean - endpoint) - 1e-10) <= 0.05:
            raise ValueError("Inverted endpoints excluded by paired exact test")
    return [float(v) for v in bounds]


def summarize(diffs, signs, exact=False):
    scores = diffs @ signs / len(signs)
    mean = float(np.mean(scores))
    out = dict(
        pairs=len(scores),
        score_pp=mean,
        positive=int(sum(scores > 0)),
        leave_one_range=[
            float(min((scores.sum() - scores) / (len(scores) - 1))),
            float(max((scores.sum() - scores) / (len(scores) - 1))),
        ],
    )
    if exact:
        shifts = np.array(list(product([-1, 1], repeat=len(scores)))) @ scores / len(scores)
        out.update(
            extreme_one=int(sum(shifts >= mean - 1e-12)),
            extreme_two=int(sum(abs(shifts) >= abs(mean) - 1e-12)),
            assignments=len(shifts),
        )
    return out


def run(root, fetch=False):
    panel = root / "evidence/lesion-transfer/panel.json"
    if sha256(panel.read_bytes()).hexdigest() != PANEL_HASH:
        raise ValueError("Changed panel")
    locked = json.loads(panel.read_text())
    targets = [MAPPING[i["entrez"]] for i in locked]
    signs = np.array([i["direction"] for i in locked])
    directory = root / "data/private/gse179640"
    matched = design(directory, fetch)
    file_hashes = []
    groups = defaultdict(list)
    cells = Counter()
    digest = None
    for patient in sorted(matched):
        for site, url, name, size in sorted(matched[patient]):
            path = acquire(directory, name, url, size, fetch)
            if path.stat().st_size != size:
                raise ValueError("Unexpected H5 size")
            file_hashes.append((name, sha256(path.read_bytes()).hexdigest()))
            digest, idx, totals, ncol = counts(path, digest, targets)
            cells[site] += ncol
            groups[patient, site].append(totals)
    pin = sha256("".join(v for _, v in sorted(file_hashes)).encode()).hexdigest()
    if pin != "2a2d62c1531355122be2de5bcfe1f748fef86a6e59e6e9e566f10880b0bc1bcb":
        raise ValueError("H5 source digest differs from the audited 19-file subset")
    people = sorted(matched)
    ranks_by_site = {
        site: np.array([np.mean([ranks(t, idx) for t in groups[p, site]], axis=0) for p in people])
        for site in SITES
    }
    core, adjacent, uterus = [ranks_by_site[k] for k in SITES]
    primary = summarize(core - adjacent, signs, True)
    primary["exact_95pct_constant_shift_pp"] = paired_shift_interval(
        (core - adjacent) @ signs / len(signs)
    )
    diagnostics = {}
    for branch in ["pool_first", "first_adjacent_only", "last_adjacent_only"]:
        if branch == "pool_first":
            values = [ranks(np.sum(groups[p, SITES[1]], axis=0), idx) for p in people]
        else:
            pos = 0 if branch == "first_adjacent_only" else -1
            values = [ranks(groups[p, SITES[1]][pos], idx) for p in people]
        diagnostics[branch] = summarize(core - np.array(values), signs, True)
    return dict(
        source_subset_digest=pin,
        compressed_input_bytes=sum(r[3] for records in matched.values() for r in records),
        features=len(totals),
        cells_per_site=dict(cells),
        counts_of_libraries={s: sum(len(groups[p, s]) for p in people) for s in SITES},
        core_minus_adjacent=primary,
        core_minus_uterus=summarize(core - uterus, signs),
        adjacent_minus_uterus=summarize(adjacent - uterus, signs),
        sensitivities=diagnostics,
        gene_effects_pp={
            name: dict(zip([s for s, _ in targets], map(float, arr.mean(axis=0)), strict=True))
            for name, arr in [
                ("core_minus_adjacent", core - adjacent),
                ("core_minus_uterus", core - uterus),
                ("adjacent_minus_uterus", adjacent - uterus),
            ]
        },
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args.repo, args.fetch), indent=2, allow_nan=False))
