"""Audit physical array reuse in the published matrisome classifier split."""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import csv
from hashlib import sha256
import json
from pathlib import Path
import urllib.request

from audit_us_bulk_sources import linked_files, native_signature
from bulk_gse6364 import fetch_prefix, scanner_fields
from open_disease.geo_microarray import sha256_file
from rome_scan_dates import HEADERS_HASH, METADATA_HASH

AUTHOR_REV = "df538583652d3359d71959437837ea3e0d1195e6"
AUTHOR = f"https://raw.githubusercontent.com/fogg-lab/characterizing-ecm-transcriptome-of-endometriosis/{AUTHOR_REV}/"
FILES = {
    "coldata.tsv": (AUTHOR + "data_prep/metadata/coldata.tsv",
                    "a2696216f28f6f80a8e7db716f280c45e4d2e4264730723fb3775b94751ebc4a"),
    "test.tsv": (AUTHOR + "condition_stratification_results/tabulated_results/all_phases_all_matrisome_test_results.tsv",
                 "35cbd873dabf6f5567f9a01a2c6fad29850a2920597328edeac762cca1f4f273"),
}
# study: raw record count, SOFT view, SOFT digest, ordered concatenated prefix digest
SERIES = {
    "GSE4888": (27, "brief", "48f89117347a351ef19c412df717af437be6bf305ffc596c04e11985f9f8f50b",
                "86555a082e3322982b17f3f595e416b56b3cea796225e9bf7bd4b9b080c534b4"),
    "GSE6364": (37, "brief", "663f7c276d364c252b36b3a2dc8d89eee37a9cac250ce662197506746c44e413",
                "73928e6928400b48b192db0ccf058374f77f36c123f0c4383002cde062da7167"),
    "GSE51981": (148, "quick", "767499b710e33598cffcc7cdaf414925630ec241d90fcc04bcfec93b903bd5d8",
                 "73d9dc597c10563744e06d02a7800bce7336e302b1f8948d870e8837c4f02909"),
    "GSE29981": (20, "brief", "228c904d5a4d984709ea1d658aba24a8e6bef1852606187917c87d3e76a4a75c",
                 "f002d1fb7980b735a23a9c166c5faf684304bcc34b088f0137b1279a53ab8933"),
}
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"


def pinned_text(path, url, digest, fetch):
    if not path.is_file() and fetch:
        with urllib.request.urlopen(url, timeout=45) as source:
            raw = source.read(1_000_001)
        if len(raw) > 1_000_000 or sha256(raw).hexdigest() != digest:
            raise ValueError("Unexpected source document")
        path.write_bytes(raw)
    if sha256_file(path) != digest:
        raise ValueError("Source document changed")


def study_prefixes(directory, study, specification, fetch):
    count, view, soft_hash, header_hash = specification
    soft = directory / (study + ".soft")
    pinned_text(soft, f"{GEO}?acc={study}&targ=gsm&view={view}&form=text", soft_hash, fetch)
    links = linked_files(soft, study, count)
    target = directory / study
    target.mkdir(exist_ok=True)

    def load(record):
        key, url = record
        path = target / (key + ".part")
        if path.is_file():
            return path.read_bytes()
        if not fetch:
            raise ValueError("Missing scan headers. Use --fetch once.")
        return fetch_prefix(url)

    with ThreadPoolExecutor(max_workers=4) as pool:
        parts = list(pool.map(load, links))
    if any(len(part) != 2048 for part in parts) or sha256(b"".join(parts)).hexdigest() != header_hash:
        raise ValueError("Header ranges or source order changed")
    output = {}
    for (key, _), part in zip(links, parts, strict=True):
        path = target / (key + ".part")
        if fetch and not path.exists():
            path.write_bytes(part)
        output[key] = part
    return output


def read_keyed(path):
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
    keyed = {r["sample_name"]: r for r in rows}
    if len(keyed) != len(rows):
        raise ValueError("Repeated classifier array identifier")
    return keyed


def audit(directory, fetch=False):
    directory.mkdir(parents=True, exist_ok=True)
    for name, (url, digest) in FILES.items():
        pinned_text(directory / name, url, digest, fetch)
    membership = read_keyed(directory / "coldata.tsv")
    predictions = read_keyed(directory / "test.tsv")
    holdout = {key for key, row in membership.items() if row["batch"] == "2"}
    if len(membership) != 216 or len(holdout) != 37 or holdout != predictions.keys():
        raise ValueError("Author preprocessing split and prediction rows disagree")

    training_images = defaultdict(list)
    test_images, test_years, calvin_years = {}, set(), set()
    cohort_counts = {}
    for study, spec in SERIES.items():
        rows = study_prefixes(directory, study, spec, fetch)
        accession_set = set(rows) & membership.keys()
        batch = {"GSE4888": "1", "GSE6364": "2", "GSE51981": "3", "GSE29981": "4"}[study]
        if accession_set != {k for k, r in membership.items() if r["batch"] == batch}:
            raise ValueError("Metadata study and classifier batch disagree")
        cohort_counts[study] = len(accession_set)
        for key in accession_set:
            _, signature, year = native_signature(rows[key])
            if key in holdout:
                if signature is None or signature in test_images:
                    raise ValueError("Unresolved or repeated scan in test set")
                test_images[signature] = key
                test_years.add(scanner_fields(rows[key])[0].year)
            elif signature is None:
                calvin_years.add(year)
            else:
                training_images[signature].append((study, key))
    if test_years & calvin_years:
        raise ValueError("Calvin training scan years overlap test archive; matching needs native pixels")
    by_phase = defaultdict(Counter)
    by_group = defaultdict(Counter)
    by_source = Counter()
    label_disagreements = phase_disagreements = 0
    correct = Counter()
    for signature, key in test_images.items():
        matches = training_images.get(signature, [])
        source_group = "+".join(sorted({study for study, _ in matches})) or "unmatched"
        by_source[source_group] += 1
        category = "matched" if matches else "unmatched"
        case = membership[key]["condition"] == "endometriosis"
        label = "case" if case else "comparator"
        by_phase[membership[key]["phase"]][category] += 1
        by_group[label][category] += 1
        label_disagreements += sum(membership[other]["condition"] != membership[key]["condition"] for _, other in matches)
        phase_disagreements += sum(membership[other]["phase"] != membership[key]["phase"] for _, other in matches)
        pred = predictions[key]
        if pred["true_label"] != str(int(case)) or pred["predicted_label"] not in {"0", "1"}:
            raise ValueError("Prediction and specimen labels disagree")
        correct[category] += pred["predicted_label"] == pred["true_label"]
    return {"author_revision": AUTHOR_REV, "included_by_accession": cohort_counts,
            "training_records": len(membership) - len(holdout), "test_records": len(test_images),
            "matched_test_scans": len(test_images.keys() & training_images.keys()),
            "matched_training_rows": sum(len(training_images[sig]) for sig in test_images if sig in training_images),
            "study_combinations": dict(by_source), "phase": dict(by_phase), "condition": dict(by_group),
            "cross_split_label_disagreements_among_matching_row_pairs": label_disagreements,
            "cross_split_phase_disagreements_among_matching_row_pairs": phase_disagreements,
            "reported_predictions_correct": dict(correct), "test_scan_years": sorted(test_years),
            "calvin_train_scan_years": sorted(calvin_years)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--fetch", action="store_true", help="Opt in to bounded public data transfers")
    parser.add_argument("--include-rome", action="store_true", help="Also cache the 22 native Rome headers for scan-order analyses")
    args = parser.parse_args()
    result = audit(args.directory, args.fetch)
    if args.include_rome:
        study_prefixes(args.directory, "GSE25628", (22, "brief", METADATA_HASH, HEADERS_HASH), args.fetch)
    print(json.dumps(result, sort_keys=True, indent=2))
