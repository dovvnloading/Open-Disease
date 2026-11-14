"""Map the Mumbai transcript FPKMs to gene-level, within-biopsy percentiles.

Only the annotation receipt and aggregate study design are printed by the CLI.
The scorer returns per-donor arrays in memory for the frozen pathway analysis.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import csv
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
from pathlib import Path
import re
import urllib.request
import zlib

import numpy as np
from scipy.stats import rankdata

from open_disease.geo_microarray import sha256_file
from swiss_epithelial_mapping import HGNC_SHA256, HGNC_URL, LOCK_FILE, LOCK_SHA256

MATRIX = "GSE153739_GT_SO_5238_Transcript_Expression_Matrix.txt.gz"
SOFT = "GSE153739_family.soft.gz"
ATTRS = "wgEncodeGencodeAttrsV19.txt.gz"
ENA = "ena-run-metadata.tsv"
DIGESTS = {
    MATRIX: "627634eb5df06434cec03b12f177b7c90e4ea54438eb468b0fc50c0df20a9f4d",
    SOFT: "d26e97f9c15e62a8141874daca4cb33305a452152d1dbb59b845be49052b49f1",
    ATTRS: "d464efa58f6b72f0fcfea646addfe121c9b79a93cc1c02148562a841e492c941",
    ENA: "dce12fbbfb8b35392814e1431c82120e4736a004e8b6db516d4df1200355511c",
}
ENA_URL = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA643835&result=read_run"
           "&fields=run_accession,experiment_accession,sample_accession,secondary_sample_accession,"
           "instrument_model,library_layout,submitted_ftp,fastq_ftp,first_public,library_name&format=tsv")
GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE153nnn/GSE153739"
URLS = {MATRIX: f"{GEO_URL}/suppl/{MATRIX}", SOFT: f"{GEO_URL}/soft/{SOFT}",
        ATTRS: f"https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/{ATTRS}", ENA: ENA_URL}
HEADER_SHA = "3d899ce0d389b6a453270432bb3e1f5206ea016691e763e1e6c6457d8ab88711"
TRANSCRIPT = re.compile(r"ENST[0-9]{11}")
FPKM_SCALE = 10**10


def inputs(directory, hgnc):
    directory = Path(directory)
    for filename, digest in DIGESTS.items():
        if sha256_file(directory / filename) != digest:
            raise ValueError(f"Unexpected {filename} source digest")
    if sha256_file(hgnc) != HGNC_SHA256:
        raise ValueError("HGNC snapshot differs from pinned release")


def sample_groups(directory):
    """Return exact source-description -> diagnosis after GEO and ENA reconciliation."""
    with gzip.open(directory / SOFT, "rt", encoding="utf-8") as source:
        blocks = source.read().split("^SAMPLE = ")[1:]
    with (directory / ENA).open(encoding="utf-8", newline="") as source:
        ena = list(csv.DictReader(source, delimiter="\t"))
    if len(ena) != 7 or len(blocks) != 7:
        raise ValueError("Expected seven libraries in both metadata sources")
    for field in ("run_accession", "experiment_accession", "secondary_sample_accession"):
        if len({r[field] for r in ena if r[field]}) != 7:
            raise ValueError("Run, experiment and BioSample keys must each be unique")
    pairs = {(r["experiment_accession"], r["sample_accession"]) for r in ena}
    if any(r["instrument_model"] != "NextSeq 500" or r["library_layout"] != "PAIRED" for r in ena):
        raise ValueError("Unexpected library design")
    groups = {}
    for block in blocks:
        def one(pattern):
            found = re.findall(pattern, block, re.MULTILINE)
            if len(found) != 1:
                raise ValueError("Missing or repeated source field")
            return found[0]

        specimen = one(r"^!Sample_description = ([0-9]+)$")
        group = one(r"^!Sample_characteristics_ch1 = group: (Control|Endometriosis)$")
        one(r"^!Sample_characteristics_ch1 = phase of menstrual cycle: Mid- proliferative$")
        one(r"^!Sample_source_name_ch1 = Eutopic endometrium$")
        exp = one(r"^!Sample_relation = SRA: https://www\.ncbi\.nlm\.nih\.gov/sra\?term=(SRX[0-9]+)$")
        sample = one(r"^!Sample_relation = BioSample: https://www\.ncbi\.nlm\.nih\.gov/biosample/(SAMN[0-9]+)$")
        if specimen in groups or (exp, sample) not in pairs:
            raise ValueError("Clinical and archive sample keys do not agree")
        pairs.remove((exp, sample))
        groups[specimen] = group == "Endometriosis"
    if pairs or sum(groups.values()) != 4:
        raise ValueError("Unexpected case or library count")
    return groups


def mapping(directory, hgnc, features):
    """Exclude ambiguous links at transcript, stable gene, and HGNC entry levels."""
    tx_genes = defaultdict(set)
    with gzip.open(directory / ATTRS, "rt", encoding="utf-8", newline="") as source:
        for row in csv.reader(source, delimiter="\t"):
            if len(row) != 14:
                raise ValueError("Unexpected GENCODE v19 field count")
            transcript = row[4].split(".")[0]
            gene = row[0].split(".")[0]
            if TRANSCRIPT.fullmatch(transcript) and re.fullmatch(r"ENSG[0-9]{11}", gene):
                tx_genes[transcript].add(gene)
    gene_entrez = defaultdict(set)
    identities = defaultdict(set)
    with Path(hgnc).open(encoding="utf-8", newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            entrez = row["entrez_id"]
            if row["status"] == "Approved" and re.fullmatch(r"[1-9][0-9]*", entrez):
                identities[int(entrez)].add(row["hgnc_id"])
                gene_entrez[row["ensembl_gene_id"]].add(int(entrez))
    mapping = {}
    observed_gene_ids = defaultdict(set)
    unique_gencode = 0
    for tx in features:
        genes = tx_genes.get(tx, set())
        if len(genes) != 1:
            continue
        gene, = genes
        unique_gencode += 1
        ids = gene_entrez[gene]
        if len(ids) == 1:
            entrez, = ids
            if len(identities[entrez]) == 1:
                mapping[tx] = entrez
                observed_gene_ids[entrez].add(gene)
    unique_ids = {entrez for entrez, stable_ids in observed_gene_ids.items() if len(stable_ids) == 1}
    usable = {tx: entrez for tx, entrez in mapping.items() if entrez in unique_ids}
    projected = np.array([usable.get(tx, 0) for tx in features], dtype=np.int64)
    return projected, unique_gencode


def fixed_value(text):
    # The source's smallest displayed decimal unit is 10^-10 FPKM, including
    # numbers in scientific notation. Integer sums preserve exact printed ties.
    try:
        scaled = Decimal(text) * FPKM_SCALE
        if not scaled.is_finite() or scaled < 0 or scaled != scaled.to_integral_value():
            raise ValueError("FPKM is invalid or has unsupported precision")
        return int(scaled)
    except InvalidOperation as exc:
        raise ValueError("Invalid deposited FPKM") from exc


def matrix(directory, groups, *, read_values):
    with gzip.open(directory / MATRIX, "rt", encoding="utf-8", newline="") as source:
        rows = csv.reader(source, delimiter="\t")
        header = next(rows)
        if header[0] != "Transcript_ID" or len(header) != 8 or set(header[1:]) != groups.keys():
            raise ValueError("Unknown or duplicated matrix specimen labels")
        features, values = [], []
        for row in rows:
            if not row:
                continue
            if len(row) != len(header) or not TRANSCRIPT.fullmatch(row[0]):
                raise ValueError("Unrecognized transcript or matrix shape")
            features.append(row[0])
            if read_values:
                values.append([fixed_value(text) for text in row[1:]])
    if len(features) != 145630 or len(set(features)) != len(features):
        raise ValueError("Transcript rows must be complete and unique")
    return features, values, np.array([groups[key] for key in header[1:]])


def gene_scores(ids, expression):
    """Rank summed transcript abundance and, separately, average transcript-first ranks."""
    selected = ids > 0
    genes, index = np.unique(ids[selected], return_inverse=True)
    # Absolute source totals are an upper bound for every positive subset sum.
    # Summing in Python integers checks int64 arithmetic bounds before add.at.
    if any(sum(map(int, column)) > np.iinfo(np.int64).max for column in expression.T):
        raise ValueError("Fixed-point FPKM exceeds the supported integer range")
    totals = np.zeros((len(genes), expression.shape[1]), dtype=np.int64)
    np.add.at(totals, index, expression[selected])
    primary = 100 * (rankdata(totals, axis=0, method="average") - 0.5) / len(genes)
    ranks = rankdata(expression, axis=0, method="average") - 0.5
    ranks_per_gene = np.zeros(totals.shape, dtype=np.float64)
    np.add.at(ranks_per_gene, index, ranks[selected])
    sensitivity = 100 * ranks_per_gene / (len(ids) * np.bincount(index)[:, None])
    return genes, primary.T, sensitivity.T


def load_mumbai(directory, hgnc):
    directory = Path(directory)
    inputs(directory, hgnc)
    groups = sample_groups(directory)
    features, values, cases = matrix(directory, groups, read_values=True)
    ids, _ = mapping(directory, hgnc, features)
    genes, primary, sensitivity = gene_scores(ids, np.array(values, dtype=np.int64))
    return genes, primary, cases, sensitivity


def read_header_audit(directory, fetch=False):
    """Keep only the first archive read-name line from a bounded gzip prefix per run."""
    if sha256_file(directory / ENA) != DIGESTS[ENA]:
        raise ValueError("Unexpected ENA metadata")
    with (directory / ENA).open(encoding="utf-8", newline="") as source:
        rows = sorted(csv.DictReader(source, delimiter="\t"), key=lambda r: r["run_accession"])
    digest = hashlib.sha256()
    for record in rows:
        run = record["run_accession"]
        path = directory / f"{run}.read-header"
        if not path.exists() and fetch:
            locator = "https://" + record["fastq_ftp"].split(";")[0]
            request = urllib.request.Request(locator, headers={"Range": "bytes=0-4095"})
            with urllib.request.urlopen(request, timeout=30) as response:
                if response.status != 206 or not re.fullmatch(r"bytes 0-4095/[0-9]+", response.headers.get("Content-Range", "")):
                    raise ValueError("ENA did not return the requested header range")
                prefix = response.read(4097)
            if len(prefix) != 4096:
                raise ValueError("Unexpected compressed prefix size")
            start = zlib.decompressobj(31).decompress(prefix, max_length=1024)
            name, newline, _ = start.partition(b"\n")
            if not newline or not name.startswith(b"@"):
                raise ValueError("Expected complete first FASTQ name line")
            path.write_bytes(name + newline)
        if not path.is_file():
            raise ValueError("Missing private header receipt. Run with --fetch-headers.")
        header = path.read_bytes()
        if not re.fullmatch((r"@" + re.escape(run) + r"\.1 1/1\n").encode(), header):
            raise ValueError("RNA read-header format differs from the audited archive snapshot")
        digest.update(header)
    if digest.hexdigest() != HEADER_SHA:
        raise ValueError("Aggregate archive header digest differs")
    return {"distinct_runs": len(rows), "bounded_bytes_per_run": 4096,
            "checked_first_mate_1_headers": len(rows), "native_flowcell_identifiers": 0,
            "submitted_fastq_links": sum(bool(r['submitted_ftp']) for r in rows),
            "ena_library_names": sum(bool(r['library_name']) for r in rows),
            "header_set_sha256": digest.hexdigest(),
            "interpretation": "Each ENA-generated first read name contains only an archive identifier and mate number. "
                              "The common NextSeq 500 platform and seven distinct archive runs do not identify shared "
                              "physical flowcells, lanes or library-preparation batches across diagnostic groups."}


def annotation_receipt(directory, hgnc):
    directory = Path(directory)
    inputs(directory, hgnc)
    groups = sample_groups(directory)
    features, _, _ = matrix(directory, groups, read_values=False)
    ids, unique_gencode = mapping(directory, hgnc, features)
    genes = set(ids) - {0}
    if sha256_file(LOCK_FILE) != LOCK_SHA256:
        raise ValueError("Frozen pathway membership changed")
    rows = json.loads(LOCK_FILE.read_text())["candidates"]
    return {
        "source_rows": len(features), "unique_gencode_transcripts": unique_gencode,
        "entrez_projected_transcripts": int(np.count_nonzero(ids)), "entrez_genes": len(genes),
        "cases": sum(groups.values()), "controls": len(groups)-sum(groups.values()),
        "pathway_coverage": {r["id"]: len(genes.intersection(r["entrez_ids"])) for r in rows if r["family"] == "pathway"},
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("hgnc", type=Path)
    parser.add_argument("--fetch-headers", action="store_true", help="Request the first 4,096 compressed bytes per missing ENA header")
    args = parser.parse_args()
    receipt = annotation_receipt(args.directory, args.hgnc)
    receipt["sources"] = {name: {"url": URLS[name], "sha256": digest, "bytes": (args.directory / name).stat().st_size}
                          for name, digest in DIGESTS.items()}
    receipt["hgnc"] = {"url": HGNC_URL, "sha256": HGNC_SHA256}
    receipt["technical_metadata"] = read_header_audit(args.directory, args.fetch_headers)
    receipt["technical_metadata"]["sequencing_provider"] = "Genotypic Technology Pvt. Ltd., as reported for cases and controls in the original article."
    receipt["mapping"] = {
        "method": "Unversioned ENST to a single GENCODE v19 stable ENSG, then a single approved HGNC Entrez identifier. "
                  "Require one approved HGNC identity per Entrez identifier and one observed mapped ENSG per Entrez identifier.",
        "primary": "Sum mapped transcript FPKMs per gene using exact fixed-point units, then midrank all mapped gene totals within biopsy.",
        "sensitivity": "Midrank all 145630 original transcripts within biopsy, including unmapped rows, then average ranks within mapped genes.",
    }
    receipt["clinical_metadata"] = {
        "article": "https://doi.org/10.1093/humrep/deaa255",
        "tissue": "Whole eutopic endometrium, histologically dated mid-proliferative phase.",
        "controls": "Study controls underwent laparoscopy for interval tubal ligation or male-factor infertility and had no visible endometriosis.",
        "medication": "Original study excluded participants reporting any medication or oral contraceptive in the previous three months. Lifetime contraceptive exposure was incompletely recorded.",
        "staging": "The original study used revised ASRM surgical staging. GEO records disease site but supplies no stage field for the sequenced subgroup.",
    }
    receipt["transfer_bytes"] = {
        "required_study_and_annotation_sources": sum(item["bytes"] for item in receipt["sources"].values()),
        "archive_read_name_ranges": receipt["technical_metadata"]["checked_first_mate_1_headers"] * 4096,
        "optional_series_metadata_and_annotation_schema": 2690 + 2098,
        "reused_hgnc_snapshot": 0,
    }
    print(json.dumps(receipt, indent=2))
