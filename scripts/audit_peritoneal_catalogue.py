"""Cross-check the published Gunther study catalogue against the prior tissue inventory."""

import argparse
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import re

from openpyxl import load_workbook


WORKBOOK_URL = "https://cdn.elifesciences.org/articles/103825/elife-103825-supp1-v1.xlsx"
WORKBOOK_SHA256 = "03f6d4130afce2184c52d45dc43d5b894048006cda351802c2fd7b331b40cde4"
PRIOR_INVENTORY = Path(__file__).resolve().parents[1] / (
    "docs/research-notes/PERITONEAL_CONTROL_AUDIT_2026-11-16.md"
)
HEADERS = (
    "Series", "Year", "Platform", "ENDOMETRIUM", "NON-LESION ENDO",
    "PRIMARY CELL CULTURE(n, phenotype)", "CELL LINE", "ENDO TISSUE",
    "BIOLOGICAL CONTROL", "CITATION", None,
)
ACCESSION = re.compile(r"\b(?:GSE\d+|E-(?:MTAB|MEXP)-\d+)\b")

# Manual adjudication from study records, distinct from the catalogue's control labels.
ADJUDICATION = {
    "GSE141549": ("Finnish EndometDB", "https://doi.org/10.1038/s41597-020-00623-x",
                  "Shared recruitment with the Finnish audit. Biopsy-level negative pathology and matching pelvic subsites remain unverified."),
    "GSE15309": ("CGAP SAGE library collection", "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15309",
                 "One archival peritoneal library also appears in GSE14. No matched superficial-lesion cohort."),
    "GSE179640": ("Tan", "https://doi.org/10.1038/s41556-022-00961-5",
                  "Adjacent case material. Noncase specimens are uterine."),
    "GSE213216": ("Fonseca", "https://doi.org/10.1038/s41588-022-01254-1",
                  "Four pathology-negative profiles from two case donors. No noncase peritoneum."),
    "GSE56414": ("Saare", "https://doi.org/10.1371/journal.pone.0112630",
                 "Small RNA sequencing of two cases and adjacent material cannot measure the frozen mRNA score."),
    "E-MTAB-694": ("Sohler", "https://doi.org/10.1177/1933719112451147",
                   "Case-only Erlangen collection. Histology and anatomical subsite of negative peritoneal biopsies remain unverified."),
    "GSE197928": ("Mortlock", "https://pmc.ncbi.nlm.nih.gov/articles/PMC9040176/",
                  "Endometrioma and adjacent stromal chromatin profiling. Wrong organ and molecular assay."),
    "GSE86534": ("Zhang", "https://pubmed.ncbi.nlm.nih.gov/30496345/",
                 "Four case donors with ovarian lesions, adjacent ovary and uterus on an mRNA/lncRNA array. No peritoneal biopsy."),
}


def screen(rows):
    """Return accession-level discovery flags from the catalogue's control column only."""
    rows = iter(rows)
    if tuple(next(rows, ())) != HEADERS:
        raise ValueError("Unexpected Included Series headers")
    seen = set()
    lesions = 0
    hits = {}
    for row in rows:
        if len(row) != len(HEADERS):
            raise ValueError("Unexpected catalogue row length")
        accession = str(row[0]).strip().removeprefix("Series ")
        if not ACCESSION.fullmatch(accession) or accession in seen:
            raise ValueError("Invalid or duplicated series accession")
        seen.add(accession)
        lesions += bool(row[7])
        control = str(row[8] or "").lower()
        flags = [term for term in ("adjacen", "periton") if term in control]
        if flags:
            hits[accession] = flags
    return len(seen), lesions, hits


def audit(path, inventory):
    with path.open("rb") as handle:
        raw = handle.read(100_001)
    if len(raw) > 100_000 or sha256(raw).hexdigest() != WORKBOOK_SHA256:
        raise ValueError("Unexpected catalogue workbook checksum or size")
    book = load_workbook(BytesIO(raw), read_only=True, data_only=True)
    try:
        n_series, n_lesion_series, hits = screen(book["Included Series"].values)
    finally:
        book.close()
    if hits.keys() != ADJUDICATION.keys():
        raise ValueError("Candidate accessions do not match the source-reviewed set")
    prior = set(ACCESSION.findall(inventory.read_text()))
    candidates = {}
    for accession, flags in sorted(hits.items()):
        cohort, paper, reason = ADJUDICATION[accession]
        record = (f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"
                  if accession.startswith("GSE")
                  else f"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/{accession}")
        candidates[accession] = {
            "catalogue_control_matches": flags, "named_in_prior_report": accession in prior,
            "collection": cohort, "source": paper, "official_record": record,
            "primary_source_disposition": reason,
        }
    peritoneal_hits = {accession for accession, flags in hits.items() if "periton" in flags}
    return {
        "catalogue": "https://doi.org/10.7554/eLife.103825",
        "catalogue_search_date": "2025-01-27", "workbook_url": WORKBOOK_URL,
        "workbook_sha256": WORKBOOK_SHA256, "workbook_bytes": len(raw),
        "sheet": "Included Series", "included_series": n_series,
        "series_with_nonempty_endo_tissue_field": n_lesion_series,
        "discovery_rule": "Case-insensitive substring adjacen or periton in BIOLOGICAL CONTROL across every series.",
        "matched_series": len(hits), "peritoneal_control_label_series": len(peritoneal_hits),
        "peritoneal_control_label_series_in_prior_report": len(peritoneal_hits & prior),
        "candidates": candidates,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path, help="Local unmodified public Gunther supplement 1")
    args = parser.parse_args()
    print(json.dumps(audit(args.workbook, PRIOR_INVENTORY), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
