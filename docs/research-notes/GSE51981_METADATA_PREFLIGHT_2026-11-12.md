# GSE51981 metadata preflight: research status

**2026-11-12 · Task #15 · Evidence: public-source audit and synthetic-tested engineering only.** No real GEO series matrix, patient record, or expression value was downloaded or analyzed for this change.

## Cohort facts that affect the next analysis

The original 2014 study reports **77 endometriosis cases, 37 other-pathology samples and 34 pathology-free controls** (148 total). In its Table 1, four participants have indeterminate phase: **two cases and two other-pathology**. The phase-assigned cohort therefore contains 75 cases, 35 other-pathology, and 34 controls (144 total). These are published reference counts, not verified per-GSM classifications by this project.

The original GEO overall-design field describes samples from **148 women**, which supports subject independence at the study-design level. It does not establish the sample-to-person mapping for every GSM. Numeric sample titles and unique GSM accessions must not be assumed to be donor identifiers. The source study says hormonal treatment in the prior three months was excluded, but that cohort-level statement is not a reviewed treatment observation for every GSM.

Primary sources: https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/ and https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981.

The later published 2022 reanalysis uses **71 cases plus 34 pathology-free controls** after additional exclusions and different processing. It is from the **same cohort**, not independent validation: https://pmc.ncbi.nlm.nih.gov/articles/PMC8766492/.

## What the new utility verifies

The standard-library script at scripts/audit_gse51981_metadata.py reads GEO Series Matrix **header metadata up to the expression-table begin marker**. It checks exact GSE51981 accession, GPL570 platforms, 148 distinct GSM records and per-sample metadata alignment. It reports only the aggregate number of source records with characteristic **field keys** for tissue, endometriosis, pathology, explicit phase, subject identity, hormone/treatment, and processing. It does not infer clinical labels from values, titles, or sample IDs. SHA-256 calculation necessarily streams the original input bytes, but **does not parse or analyze expression intensities**.

If the user independently created a local reviewed manifest, the optional check verifies the same 148 accessions, explicit review/inclusion decisions, exclusion reasons, independent donor keys *as declared*, supported group/phase/tissue/treatment strings, and a matching original GSM URL for each included record. These structural checks **cannot authenticate donor identity, labels, informed consent or legal rights**; approval always remains NO_GO until separately documented source verification.

The input may be the source's local, unmodified processed file (never put individual records in Git). Example after explicit public data acquisition and permission review:

    python scripts/audit_gse51981_metadata.py \
      --matrix data/raw/GSE51981_series_matrix.txt.gz \
      --out results/gse51981/header_preflight.json

To preflight an independently reviewed, local-only manifest:

    python scripts/audit_gse51981_metadata.py \
      --matrix data/raw/GSE51981_series_matrix.txt.gz \
      --manifest data/private/GSE51981_review.csv \
      --out results/gse51981/review_preflight.json

Outputs are aggregate JSON under ignored results/, with compressed-file SHA-256, no overwrite, no individual clinical/source text or expression values. Run synthetic tests without network:

    python -m pytest -q tests/test_audit_gse51981_metadata.py

## Scientific gate

**NO-GO for #18 until** the actual 148 GSM records have source-backed group, tissue, phase and treatment review; the sample-to-donor identity assumption is defensible; public data use and hashes are checked; and the comparison and gene-mapping policy are frozen *before* inspecting outcomes. Ambiguous records are excluded with explicit reasons rather than guessed. This utility **does not close #15**. No grants, data-custodian contact, clinician recruitment or separate issue is requested.