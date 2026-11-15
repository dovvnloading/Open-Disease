"""Challenge the frozen ovarian lesion panel in matched spatial compartments."""

import argparse
from collections import Counter, defaultdict
import gzip
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import re
import tarfile
from urllib.request import urlopen
from tempfile import TemporaryDirectory
import numpy as np
from openpyxl import load_workbook
from scipy.stats import rankdata

from cross_cohort_discovery import holm

SOFT = "GSE263897_family.soft.gz"
PKC = "GSE263897_Hs_R_NGS_WTA_v1.0.pkc.gz"
DCC = "GSE263897_RAW.tar"
PANEL_HASH = "4c634a65d648a83af8418e7ad1bf50e03d1a83982d0265d47d9ccbd06aba0114"
SOURCES = {
    SOFT: "b68340e69de0f28b233d524f87fd4b228d4d8ae2dd22950c3ac9067701f8a778",
    PKC: "a146d7b9ccac568737128213f2492ae4459e005b5820d2d9a4e87327f1651e11",
    DCC: "e8e9b604e004cec73772fad8838bc1cc2bc70c2608de4cefdb5c3354f7f0bca8",
    "mmc2.xlsx": "acf4cfc547f9aaacfa44b907f6e8ef02e2fce70888137eeb0254735cdb2993cf",
    "mmc3.xlsx": "7ab22a2ae7a4cf320930150faf12ed9ef9f70c966f62fb409c3cb83d64f4cf32",
    "mmc4.xlsx": "56067fb4404924e3104f03b51ff667fb1db65f159585b5dd754da3df769ad9d4",
}
TISSUES = ("eutopic endometrium", "endometriotic lesion")
CELL_TYPES = ("Stroma", "Epithelium", "Macrophages")


def verify_sources(directory, fetch):
    geo = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE263nnn/GSE263897"
    for name, digest in SOURCES.items():
        path = directory / name
        url = (
            f"https://pmc-oa-opendata.s3.amazonaws.com/PMC11810701.1/{name}"
            if name.endswith(".xlsx")
            else f"{geo}/{'soft' if name == SOFT else 'suppl'}/{name}"
        )
        if fetch and not path.is_file():
            with urlopen(url, timeout=60) as f:
                data = f.read(4_000_001)
            if len(data) > 4_000_000 or sha256(data).hexdigest() != digest:
                raise ValueError(f"Source changed: {name}")
            path.parent.mkdir(parents=True, exist_ok=True)
            with TemporaryDirectory(dir=path.parent) as staging:
                temporary = Path(staging) / name
                temporary.write_bytes(data)
                temporary.replace(path)
        if sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Unverified input: {name}")


def design(path):
    text = gzip.open(path, "rt").read()
    records = {}
    for block in text.split("^SAMPLE = ")[1:]:
        fields = dict(
            v.split(": ", 1)
            for v in re.findall(r"^!Sample_characteristics_ch1 = (.*)$", block, re.M)
        )
        (filename,) = re.findall(
            r"^!Sample_supplementary_file_1 = .*/([^/]+\.dcc\.gz)$", block, re.M
        )
        if (
            filename in records
            or fields["cell type"] not in CELL_TYPES
            or fields["tissue"] not in TISSUES
        ):
            raise ValueError("Unexpected or duplicate region")
        records[filename] = (fields["sampleID"], fields["tissue"], fields["cell type"])
    cells = Counter(records.values())
    if (
        len(records) != 60
        or len(cells) != 30
        or set(cells.values()) != {2}
        or len({x[0] for x in cells}) != 5
    ):
        raise ValueError("Incomplete five-pair spatial design")
    return records


def configuration(path, panel):
    gene_ids, negatives, symbols, indices = [], [], {}, {}
    for target in json.load(gzip.open(path))["Targets"]:
        probes = target["Probes"]
        if target["CodeClass"] == "Negative01":
            negatives.extend(p["RTS_ID"] for p in probes)
            continue
        if target["CodeClass"] != "Endogenous01" or len(probes) != 1:
            raise ValueError("Unexpected spatial target class")
        (probe,) = probes
        for gene in panel:
            if gene["entrez"] in probe.get("GeneID", []):
                if len(probe["GeneID"]) != 1 or gene["entrez"] in indices:
                    raise ValueError("Ambiguous panel probe")
                symbols[gene["entrez"]] = target["DisplayName"]
                indices[gene["entrez"]] = len(gene_ids)
        gene_ids.append(probe["RTS_ID"])
    all_ids = gene_ids + negatives
    if (
        len(gene_ids) != 18676
        or len(negatives) != 139
        or len(set(all_ids)) != len(all_ids)
        or len(indices) != 5
    ):
        raise ValueError("Inconsistent panel/whole-transcriptome reference")
    return (
        all_ids,
        len(gene_ids),
        [indices[k["entrez"]] for k in panel],
        [symbols[k["entrez"]] for k in panel],
    )


def assay(directory, metadata, ids, genes, indices):
    position = {identifier: i for i, identifier in enumerate(ids)}
    results = []
    with tarfile.open(directory / DCC) as tar:
        members = tar.getmembers()
        if len(members) != len(metadata) or {m.name for m in members} != metadata.keys():
            raise ValueError("Archive and GEO area identities disagree")
        for member in members:
            if not member.isfile() or member.size > 200000:
                raise ValueError("Unexpected DCC archive member")
            data = gzip.decompress(tar.extractfile(member).read()).decode("utf-8")
            block = data.split("<Code_Summary>")[1].split("</Code_Summary>")[0]
            counts = np.zeros(len(ids), dtype=int)
            seen = set()
            for line in block.strip().splitlines():
                name, raw = line.strip().split(",")
                if name in seen or name not in position or not raw.isdigit():
                    raise ValueError("Invalid DCC probe or count")
                seen.add(name)
                counts[position[name]] = int(raw)
            endogenous = counts[:genes]
            negative = np.log(np.maximum(1, counts[genes:]))
            loq = max(2.0, np.exp(negative.mean() + 2 * negative.std(ddof=1)))
            passing = np.mean(endogenous > loq) >= 0.1
            values = 100 * (rankdata(endogenous) - 0.5) / genes
            results.append(
                (metadata[member.name], passing, values[indices], endogenous[indices] > loq)
            )
    return results


def summarize_regions(records, symbols, signs, filtered):
    groups = defaultdict(list)
    detection = {c: [] for c in CELL_TYPES}
    for (patient, tissue, cell), passes, values, detected in records:
        if not filtered or passes:
            groups[patient, tissue, cell].append(values)
            detection[cell].append(detected)
    results = {}
    for cell in CELL_TYPES:
        diffs = []
        for patient in sorted({p for p, _, c in groups if c == cell}):
            lhs = groups.get((patient, TISSUES[0], cell))
            rhs = groups.get((patient, TISSUES[1], cell))
            if lhs and rhs:
                diffs.append(np.mean(rhs, axis=0) - np.mean(lhs, axis=0))
        unit_scores = np.array(diffs) @ np.array(signs) / len(signs)
        n = len(unit_scores)
        obs = float(np.mean(unit_scores))
        flips = np.array(list(product([-1, 1], repeat=n))) @ unit_scores / n
        extreme = int(np.sum(flips >= obs - 1e-12))
        results[cell] = dict(
            pairs=n,
            regions=len(detection[cell]),
            positive_pairs=int(np.sum(unit_scores > 0)),
            effect_pp=obs,
            extreme=extreme,
            assignments=len(flips),
            one_sided_p=extreme / len(flips),
            gene_effects_pp=dict(zip(symbols, np.mean(diffs, axis=0), strict=True)),
            regions_above_loq=dict(
                zip(symbols, map(int, np.sum(detection[cell], axis=0)), strict=True)
            ),
        )
    for row, p in zip(
        results.values(), holm([r["one_sided_p"] for r in results.values()]), strict=True
    ):
        row["holm_p"] = float(p)
    return results


def summary_tables(directory, symbols):
    output = {}
    family = []
    for cell, filename in zip(CELL_TYPES, ["mmc2.xlsx", "mmc3.xlsx", "mmc4.xlsx"], strict=True):
        book = load_workbook(directory / filename, read_only=True, data_only=True)
        try:
            rows = iter(book.active.values)
            header = next(rows)
            if header != ("Gene_name", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"):
                raise ValueError("Source columns changed")
            records = {r[0]: r for r in rows}
        finally:
            book.close()
        if (
            len(records) != 7946
            or "NegProbe-WTX" not in records
            or (cell == "Epithelium" and records["C3"][1] <= 0)
        ):
            raise ValueError("Incorrect gene universe or source contrast")
        output[cell] = {}
        for name in symbols:
            if name in records:
                row = records[name]
                if not np.isfinite(row[1:]).all() or not (0 <= row[4] <= 1 and 0 <= row[5] <= 1):
                    raise ValueError("Invalid summary statistic")
                result = dict(log2FC=row[1], raw_p=row[4], source_bh=row[5])
            else:
                result = dict(log2FC=None, raw_p=None, source_bh=None)
            output[cell][name] = result
            family.append(result)
    for row, adj in zip(
        family, holm([1 if r["raw_p"] is None else r["raw_p"] for r in family]), strict=True
    ):
        row["targeted_holm_p"] = float(adj) if row["raw_p"] is not None else None
    return output


def run(root, fetch=False):
    path = root / "evidence/lesion-transfer/panel.json"
    if sha256(path.read_bytes()).hexdigest() != PANEL_HASH:
        raise ValueError("Unexpected upstream panel")
    panel = json.loads(path.read_text())
    directory = root / "data/private/gse263897"
    verify_sources(directory, fetch)
    metadata = design(directory / SOFT)
    ids, genes, indices, symbols = configuration(directory / PKC, panel)
    records = assay(directory, metadata, ids, genes, indices)
    signs = [g["direction"] for g in panel]
    return dict(
        reference_genes=genes,
        negative_probes=len(ids) - genes,
        panel=dict(zip(symbols, signs, strict=True)),
        source_model=summary_tables(directory, symbols),
        all_regions=summarize_regions(records, symbols, signs, False),
        detection_filtered=summarize_regions(records, symbols, signs, True),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path("."))
    parser.add_argument("--fetch", action="store_true", help="Acquire missing pinned public files")
    args = parser.parse_args()
    print(json.dumps(run(args.repo, args.fetch), indent=2, allow_nan=False))
