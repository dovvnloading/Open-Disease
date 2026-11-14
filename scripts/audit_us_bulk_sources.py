"""Replay the GSE6364 scan-day and cross-archive physical array audits."""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import struct
import zlib

from bulk_gse6364 import (BRIEF_URL, MATRIX_URL, checked_header, checked_platform,
                          fetch_prefix, scan_audit, scanner_fields)
from open_disease.geo_microarray import sha256_file
from reconcile_gse51981_sources import (CHAR, DIAGNOSIS, GROUP_CODES, PATHOLOGY, SEVERITY,
                                       TISSUE, read_quick_view)


GPL_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
DATASETS = (
    ("GSE4888", "brief", "48f89117347a351ef19c412df717af437be6bf305ffc596c04e11985f9f8f50b",
     27, "86555a082e3322982b17f3f595e416b56b3cea796225e9bf7bd4b9b080c534b4"),
    ("GSE51981", "quick", "767499b710e33598cffcc7cdaf414925630ec241d90fcc04bcfec93b903bd5d8",
     148, "73d9dc597c10563744e06d02a7800bce7336e302b1f8948d870e8837c4f02909"),
)


def file_record(path, url):
    return {"url": url, "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def linked_files(metadata_path, study, count):
    records = []
    current = None
    for line in metadata_path.read_text(encoding="utf-8").splitlines():
        key, sep, value = line.partition(" = ")
        if not sep:
            continue
        if key == "^SAMPLE":
            current = value
        if key == "!Sample_supplementary_file":
            if not current or not re.fullmatch(r"GSM[0-9]+", current):
                raise ValueError(f"Invalid {study} sample key")
            if not value.startswith("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/") or not value.endswith(".CEL.gz"):
                raise ValueError(f"Invalid {study} native source URL")
            records.append((current, value.replace("ftp://", "https://", 1)))
    if (len(records) != count or len({k for k, _ in records}) != count
            or len({url for _, url in records}) != count):
        raise ValueError(f"Unexpected {study} raw CEL links")
    return records


def native_signature(prefix):
    """Full native DAT image header and four grid corners or Calvin's actual scan year."""
    if len(prefix) != 2048:
        raise ValueError("Invalid compressed CEL range")
    raw = zlib.decompressobj(31).decompress(prefix, 32768)
    if raw.startswith(b"\x3b\x01"):
        name = "affymetrix-scan-date".encode("utf-16-be")
        if raw.count(name) != 1:
            raise ValueError("Unrecognized Calvin scan date")
        after = raw.split(name)[1]
        length = struct.unpack_from(">I", after)[0]
        date = datetime.fromisoformat(after[4:4 + length].decode("utf-16-be"))
        return "Calvin", None, date.year
    if raw.startswith(b"[CEL]"):
        fmt = "text"
        header, boundary, _ = raw.partition(b"[INTENSITY]")
        if not boundary:
            raise ValueError("Truncated text CEL header")
    elif raw.startswith(b"@\x00\x00\x00"):
        fmt = "XDA_v4"
        magic, ver, cols, rows, entries, size = struct.unpack_from("<6i", raw)
        if (magic, ver, cols, rows, entries) != (64, 4, 1164, 1164, 1164**2) or not 0 < size <= len(raw) - 24:
            raise ValueError("Invalid XDA header size or format")
        header = raw[24:24 + size]
    else:
        raise ValueError("Unsupported native CEL version")
    fields = [b"DatHeader", b"GridCornerUL", b"GridCornerUR", b"GridCornerLR", b"GridCornerLL"]
    values = {}
    for line in header.splitlines():
        key, sep, value = line.partition(b"=")
        if sep and key in fields:
            if key in values or not value.strip():
                raise ValueError("Ambiguous native CEL identity field")
            values[key] = value.strip()
    if set(fields) != values.keys():
        raise ValueError("Missing native CEL identity fields")
    return fmt, tuple(values[key] for key in fields), None


def matched_labels(study, records, matched_keys, clinical):
    counts = Counter()
    volunteers = set()
    for identity, accession in matched_keys.items():
        chars = records[accession][CHAR]
        if study == "GSE4888":
            counts["proliferative"] += "Proliferative phase PE" in chars
            diagnostic = [entry.lower() for entry in chars if entry.startswith("Age: ")]
            if len(diagnostic) != 1:
                raise ValueError("Ambiguous GSE4888 clinical metadata")
            if "nml volunteer" in diagnostic[0]:
                counts["normal_volunteer"] += 1
                volunteers.add(identity)
            elif "leiomyomata" in diagnostic[0]:
                counts["leiomyomata"] += 1
            else:
                counts["other_or_unresolved"] += 1
        else:
            fields = dict(entry.split(": ", 1) for entry in chars)
            case = fields[DIAGNOSIS] == "Endometriosis"
            counts["proliferative"] += fields[TISSUE] == "Proliferative Endometrial tissue"
            counts["case_control_agreement"] += case == clinical[identity]
            if case:
                counts["moderate_severe"] += fields.get(SEVERITY) == "Moderate/Severe"
            else:
                group = GROUP_CODES[fields[DIAGNOSIS], fields[PATHOLOGY]]
                counts["control_other_pelvic_pathology"] += group == "other"
                counts["control_no_pelvic_pathology"] += group == "control"
    return dict(counts), volunteers


def source_audit(directory, annotation, talbi, larger, fetch=False):
    _, selected = checked_header(directory)
    scan = scan_audit(directory, selected, fetch)
    probes = checked_platform(annotation)
    sources = {"series.txt.gz": file_record(directory / "series.txt.gz", MATRIX_URL),
               "brief.soft": file_record(directory / "brief.soft", BRIEF_URL),
               "GPL570.annot.gz": file_record(annotation, GPL_URL)}
    clinical = {}
    us_prefixes = set()
    us_years = set()
    for _, accession, case, _ in selected:
        part = (directory / (accession + ".header.gzpart")).read_bytes()
        us_prefixes.add(part)
        us_years.add(scanner_fields(part)[0].year)
        _, identity, _ = native_signature(part)
        if identity in clinical:
            raise ValueError("Repeated physical scan within the selected US cohort")
        clinical[identity] = case
    studies = {}
    hits = []
    volunteers = set()
    for (study, view, metadata_sha, count, prefix_sha), metadata_path in zip(DATASETS, (talbi, larger), strict=True):
        if sha256_file(metadata_path) != metadata_sha:
            raise ValueError(f"Unexpected {study} GEO metadata")
        sources[study + "_" + view + ".soft"] = file_record(metadata_path, f"{GEO}?acc={study}&targ=gsm&view={view}&form=text")
        paths = linked_files(metadata_path, study, count)
        work = directory / study.lower()
        if fetch:
            work.mkdir(exist_ok=True)
        formats, keys, calvin_years = Counter(), Counter(), Counter()
        digest = hashlib.sha256()
        same_prefixes = 0
        matched_keys = {}
        for accession, url in paths:
            dest = work / (accession + ".header.gzpart")
            if not dest.exists() and fetch:
                dest.write_bytes(fetch_prefix(url))
            part = dest.read_bytes()
            fmt, identity, year = native_signature(part)
            formats[fmt] += 1
            digest.update(part)
            same_prefixes += part in us_prefixes
            if identity is None:
                calvin_years[year] += 1
            else:
                keys[identity] += 1
                if identity in clinical:
                    matched_keys[identity] = accession
        if digest.hexdigest() != prefix_sha or us_years & calvin_years.keys():
            raise ValueError(f"Unverified {study} source ranges or ambiguous Calvin scan dates")
        matches = clinical.keys() & keys.keys()
        if any(keys[ident] != 1 for ident in matches):
            raise ValueError(f"Nonunique matched physical scan in {study}")
        labels, study_volunteers = matched_labels(study, read_quick_view(metadata_path), matched_keys, clinical)
        volunteers.update(study_volunteers)
        hits.append(matches)
        studies[study] = {
            "cross_archive_clinical_labels": labels, "array_headers": count, "formats": dict(formats), "calvin_native_scan_years": dict(calvin_years),
            "compressed_header_bytes": count * 2048, "ordered_compressed_header_sha256": digest.hexdigest(),
            "byte_identical_compressed_prefixes": same_prefixes, "native_image_matches": len(matches),
            "matched_gse6364_cases": sum(clinical[ident] for ident in matches),
            "matched_gse6364_controls": sum(not clinical[ident] for ident in matches),
        }
    return {
        "sources": sources, "primary_paper": "https://doi.org/10.1210/en.2006-1692",
        "clinical_context": "Burney et al. report surgically ascertained disease and historical comparator specimens from the earlier Talbi et al. study. GSE6364 group and cycle-phase fields were verified against its individual GEO records.",
        "source_assay": "Individual uterine whole-tissue samples. GCOS v1.2 with within-array trimmed mean intensity scaled to 500.",
        "reference_probes": len(probes), "uniquely_mapped_probes": sum(v > 0 for v in probes.values()),
        "entrez_genes": len(set(probes.values()) - {0}), "native_scan_audit": scan,
        "comparison_method": "Use bytes 0–2047 of each original .CEL.gz linked from the sampled GEO records, ordered as in source metadata. Compare complete native DAT header text and all four grid corners after decoding CEL and trimming surrounding whitespace. Calvin files in GSE51981 record nonoverlapping scan years. A match identifies a shared physical array scan. Unmatched scans do not establish distinct donors.",
        "normal_volunteers_among_scans_matching_both_archives": len(volunteers & set.intersection(*hits)),
        "gse6364_scan_years": sorted(us_years), "cross_archive_native_scans": studies,
        "gse6364_matching_either_archive": len(set.union(*hits)), "gse6364_matching_both_archives": len(set.intersection(*hits)),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("annotation", type=Path)
    parser.add_argument("talbi_brief", type=Path)
    parser.add_argument("gse51981_quick", type=Path)
    parser.add_argument("--fetch-headers", action="store_true", help="Fetch missing bounded native header ranges from GEO")
    args = parser.parse_args()
    print(json.dumps(source_audit(args.directory, args.annotation, args.talbi_brief, args.gse51981_quick, args.fetch_headers),
                     indent=2, ensure_ascii=False))
