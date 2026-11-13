# GSE51981: auditable public microarray path (local and opt-in)

**Status:** software pathway with **synthetic-only tests**. The maintainers have not yet downloaded, manually phenotype-reviewed, or analyzed GSE51981 with these scripts. There are **no biological findings, diagnostic claims, or validated targets** from this feature. This is engineering infrastructure for [the first real-data milestone (#18)](https://github.com/dovvnloading/Open-Disease/issues/18), **not** evidence that its scientific acceptance criteria are met.

GSE51981 is a public Affymetrix GPL570 endometrial **microarray** study. The cited original paper describes 77 endometriosis cases, 34 pathology-free controls, and 37 *other pelvic-pathology* samples among 148 arrays, but these are **published cohort totals, not machine-verified per-GSM labels**. A GSM is a sample accession, not necessarily an independent person. GEO values are processed arrays, **not RNA-seq/edgeR counts**.

Sources:

- [GEO GSE51981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981)
- [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/)
- [NCBI's programmatic access convention](https://www.ncbi.nlm.nih.gov/geo/info/geo_paccess.html)
- [GEO submitter metadata specification](https://www.ncbi.nlm.nih.gov/geo/info/soft.html)

## 1. Verify rights, then optionally acquire the pinned processed matrix

Do not use the full CEL archive for this process. Check the primary GEO record, attached platform and publication for source processing, permitted reuse, technology and metadata quality *before* interpreting results. GEO public availability is not a substitute for reviewing data rights. No data are downloaded by installation, the default CLI, or CI; the following explicit action requires an Internet connection:

```bash
python -m pip install -e '.[dev]'
open-disease geo-download --out data/raw/GSE51981_series_matrix.txt.gz
```

The downloader uses the HTTPS URL:

`https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51981/matrix/GSE51981_series_matrix.txt.gz`

It rejects non-gzip responses, enforces a 150 MB compressed transfer limit, streams an SHA-256 digest, and refuses overwrites. A same-directory `*.provenance.json` sidecar records the URL, UTC retrieval time, compressed byte count and digest. If GEO changes its delivery layout, confirm the new official source independently rather than silently substituting files. Save the sidecar as part of the **local** reproducibility bundle. All raw and derived human sample data must remain in ignored, controlled directories.

Alternatively, obtain the official matrix manually after verifying the source URL and run the completely offline step below. In that case, the prep audit still computes its exact file SHA-256 (but no independent retrieval timestamp is available).

## 2. Extract without inventing sample diagnoses

```bash
open-disease geo-prepare \
  --matrix data/raw/GSE51981_series_matrix.txt.gz \
  --manifest data/private/GSE51981_review.csv \
  --expression data/processed/GSE51981_intensities.tsv \
  --audit data/processed/GSE51981_audit.json
```

This parser requires the actual `!Series_geo_accession`, the exact `!Sample_geo_accession` order, 148 unique GSM columns, unique probe identifiers, a complete matrix table, and finite numeric intensities or explicit missing values. It preserves submitter titles/characteristics as **source text only**, creates a probe-by-GSM table, and records sample/probe/missing-value totals plus the compressed input SHA-256. It makes **no** case/control, phase, tissue, or donor inference; every row starts `include=false`, `reviewed=false`, and `exclusion_reason=pending_manual_review`.

Missing intensities remain `NA`, never zero. All three outputs are local and are refused if the destination already exists. They are *not* committed or published in CI. The 148-sample check is an integrity gate for this chosen study; any mismatch requires human investigation of GEO file/platform versions, not an unquestioned `--expected-samples` override.

## 3. Create a defensible reviewed manifest *before* inspecting effects

Source-audit sample attributes and the original paper. An agent can prepare a citation-backed map, but included values should be marked reviewed only when sources genuinely verify them; unresolved fields require exclusion. An independent second check is recommended. In the ignored `data/private/GSE51981_review.csv`, record every GSM and complete the following fields only when verified:

| Field | Required handling |
| --- | --- |
| `geo_accession` | Preserve the original GSM ID; do not add/remove samples silently |
| `donor_id` | Independently verified, pseudonymous **person-level** ID; repeated donors are refused, and a GSM must not be assumed to equal a donor |
| `include`, `reviewed` | Set to `true` only after evidence review; excluded rows must have `exclusion_reason` |
| `disease_group` | Exactly `case`, `control`, or `other`; the 37 other-pathology samples are **never** healthy controls by default |
| `phase` | Verified, harmonized menstrual-phase category; no guessed or imputed cycle phase |
| `tissue` | Exactly `eutopic_endometrium` for this protocol; other tissues need another model |
| `pathology` | Document the verified specific/other pathology, including pathology-free controls |
| `hormone_treatment` | Exactly `none` for included samples; unknown/exposed samples require a separate protocol |
| `review_source_url` | HTTPS URL to GEO GSM / paper / eligible clinical metadata supporting the decision |
| `review_notes` | Source ambiguity, harmonization decisions and explicit exclusion rationale |

The `source_*` columns are **untrusted submitter text** for human review, not fitted covariates. The script cannot establish whether each submitted sample is independent or each published group label maps correctly to its GSM: an evidence audit must substantiate those assertions; an independent cross-check is preferred. If metadata are insufficient, leave samples excluded, document the blocker and **do not run** an underpowered or mislabeled analysis. Do not place identifying clinical records in the manifest or source repository.

Before the analysis, freeze the selected phase definitions, inclusion/exclusion rules, sample independence audit, primary case–control question, expected known expression theme, excluded probes, outcome interpretation and multiplicity family. A dated, frozen pre-analysis note and source-backed cross-check are needed for this exploratory milestone; this is not formal scientific preregistration.

## 4. Fit the separate limma pathway on a reviewed log2 array scale

Install [Bioconductor limma](https://bioconductor.org/packages/limma/) in a suitable R environment, confirm **from GEO processing documentation** that the prepared intensities are normalized, comparable **log2 microarray intensities**, and only then explicitly acknowledge this with `--log2-confirmed yes`:

```bash
Rscript scripts/r/analyze_geo_microarray.R \
  --expression data/processed/GSE51981_intensities.tsv \
  --manifest data/private/GSE51981_review.csv \
  --out-dir results/GSE51981-reviewed-v1 \
  --log2-confirmed yes \
  --min-group 3 \
  --seed 20261112
```

`--log2-confirmed yes` is an **operator attestation**, not automatic normalization verification. The model fails closed on missing/ambiguous diagnoses, unreviewed included samples, duplicate donor IDs, incompatible tissue/treatment, absent case/control overlap within phase, non-estimable design, insufficient independent donors and mismatched GSM IDs.

**Pre-specified exploratory model:** log2 probe intensity ~ group + menstrual-phase fixed effect (phase term omitted only when all reviewed included samples have one shared phase). The first contrast is `case - control`. If at least three independent `other` donors with shared-phase control comparisons are reviewed, the separate `other - control` contrast estimates *nonspecific pathology effects*, **not** a new control arm. The additive model assumes a constant group effect across phases; it does not test group-by-phase interactions, and it is not a reimplementation of the paper's exact phase-specific analysis.

- Each probe with **any** missing value among included samples is marked `excluded_missing_intensity` in the full output. No missing-value imputation, silent feature omission or outcome-based filtering.
- `case_vs_control.csv` and optional `other_vs_control.csv` report probe-level limma log2 fold differences, moderated tests, 95% confidence intervals, BH q values **within each contrast**, and a separate BH q across *both* planned contrast families. All input probes appear, including those excluded for missingness.
- `phase_permuted_control.csv` uses one fixed-seed, **within-phase shuffled diagnosis** assignment as a negative-control stress check, *not* external replication or a calibrated permutation FDR.
- `phase_group_counts.csv`, `sample_review_audit.csv`, `analysis_report.txt`, and `session_info.txt` capture cohort counts, exclusions, design, exact software version and the negative-control diagnostic. Keep all local; they are sample-derived.
- Probe IDs are **not gene IDs**. Gene annotation must use a versioned GPL570 map with one-to-many ambiguity retained before any gene-level or pathway claim. No target ranking is executed here.

## 5. Honest release / remaining science gates

All reported results must disclose the exact reviewed donor/sample counts, possible source selection, original publication overlap, hormone/cycle confounding, individual sample exclusions, nulls and missing probe rates, scale assumption, the complete false-discovery family and failures of the phase-permuted check. The other-pathology contrast is only run if its required reviewed records exist; otherwise the report **explicitly says it was not tested**.

A public **exploratory computational report** requires a real controlled execution, source-backed phenotype and reuse-rights audit, provenance, predeclared comparison, careful interpretation and a genuinely independent compatible comparison **or an explicit validation gap**. Independent specialist review is needed before promoting a result as validated biological or clinically consequential evidence, not before honestly reporting preliminary computation. The tests here use invented sample and probe IDs exclusively. No cure, diagnosis, causal target or clinical use is supported.

For offline verification:

```bash
python -m pytest -q tests/test_geo_microarray.py
# In an R environment with limma installed:
Rscript tests/r/test_geo_microarray.R
```

## Public metadata audit status (2026-11-12)

Source-only first-pass reports are available for [GSE51981 sample labels and donor independence](research-notes/GSE51981_SOURCE_AUDIT_2026-11-12.md) and [GPL570 probe mapping](research-notes/GPL570_ANNOTATION_AUDIT_2026-11-12.md). The original 2014 paper reports **four indeterminate cycle-phase samples**, while a later 2022 reanalysis used **71 cases + 34 no-pathology controls** after exclusions; these are **not** proof of a finished per-GSM manifest and the latter is not an independent cohort. Original GEO per-sample `VALUE` describes log2 GC-RMA, while the 2022 reanalysis independently normalized/adjusted arrays. Do not silently equate the preprocessing.

The annotation audit command takes a locally sourced GPL570 table and classifies ambiguous or missing mappings without publishing probe-level or human-data results:

```bash
python scripts/audit_gpl570.py --input data/raw/GPL570_platform_table.tsv \
  --out results/gpl570/annotation_audit.json
python -m pytest -q tests/test_audit_gpl570.py
```

**Source checks and 148-record audit remain incomplete.** The first actual GSE51981 limma analysis remains gated and should not be described as completed.


## Strict metadata-only preflight (before expression inspection)

An optional **offline, standard-library** preflight can inspect the official local GSE51981 Series Matrix **header**, without parsing expression rows. It rejects an incorrect GEO accession, platform mismatch, non-148/duplicate GSM set or malformed metadata and only emits **aggregate** source-field coverage. The input file is still fully byte-hashed (SHA-256); this is not a gene-level analysis and cannot establish clinical truth, source rights or independent donors.

    python scripts/audit_gse51981_metadata.py \
      --matrix data/raw/GSE51981_series_matrix.txt.gz \
      --out results/gse51981/header_preflight.json

After an individually **source-reviewed** manifest is generated, the same command can add "--manifest data/private/GSE51981_review.csv" and report structural problems without making guesses. The tool **always returns NO_GO for clinical/scientific approval**; the operator's substantive donor, pathology, phase, reuse and original-result review must remain separate. Its output is aggregate-only under ignored results/, with no source titles, GSM IDs or individual medical information. It does not automatically download data.

Published 2014 phase-eligible group counts are **75 cases / 35 other pathology / 34 pathology-free controls**, plus two indeterminate cases and two indeterminate other-pathology samples. This is a literature reference, **not** a verified 144-person manifest produced by the tool. The [2022 same-cohort reanalysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8766492/) used 71 cases and 34 controls after different exclusions and processing and does not constitute external validation.

See the [source audit](research-notes/GSE51981_METADATA_PREFLIGHT_2026-11-12.md), [pre-analysis plan](research-notes/GSE51981_PREANALYSIS_PLAN_2026-11-12.md) and [19 synthetic tests](../tests/test_audit_gse51981_metadata.py). These additions advance #15/#18 but **do not close either issue**.
