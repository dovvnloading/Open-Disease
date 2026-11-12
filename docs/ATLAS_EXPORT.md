# Offline Seurat atlas export — research ingestion v0.1

**Status:** Independent, locally runnable adapter; as of 12 November 2026 **NOT YET TESTED on the GSE266265 human atlas**. The base Matrix aggregation is tested against completely synthetic data in R CI. This is a software engineering component, not completion of the source-rights, phase-stratification or model-validation gates.

## Purpose

The 2026 60-donor endometrium atlas is distributed as an approximately 39 GB Seurat RDS. We must turn its **full raw RNA feature × cell count layer** into one complete, pseudonymous **donor × cell-type × gene** count table without mistaking thousands of cells for thousands of patients.

Use a controlled high-memory R environment with a verified, locally obtained data object. The adapter deliberately does not fetch any human files, copy upstream analysis scripts, label donor phenotypes automatically or print genomic counts to CI logs.

## Explicit local sample map — REQUIRED

Prepare a private, authorized CSV file outside the repository (for example data/private/sample_mapping.csv) containing exactly one mapping row per original sample:

| Column | Example value, fictional | Explanation |
| --- | --- | --- |
| sample_id | SYN_SAMPLE_01 | Must match the Seurat cell metadata sample field |
| study_id | GSE266265 | Source study identifier |
| donor_id | SYN_DONOR_01 | Pseudonymous donor identifier; must reflect genuine biological participant, not a cell or sample |
| diagnosis | case | Exactly case or control based on reviewed clinical metadata |
| phase | proliferative | Exactly proliferative or secretory for this first baseline |
| hormone_treatment | none | Only untreated participants in first baseline |

Do not invent patient labels, infer donor independence from sample name patterns or publish this mapping file to GitHub. The same donor may have multiple sample records but may not change diagnosis, phase or treatment within a study. Confirm the exact assay metadata columns and field meanings with the source paper before running.

## Example (only after independent rights and metadata review)

    Rscript scripts/r/export_seurat_pseudobulk.R \
      --rds=/approved/EndoAtlas.rds \
      --mapping=data/private/sample_mapping.csv \
      --out=data/processed/donor_pseudobulk.csv \
      --audit=data/processed/donor_export_audit.csv \
      --assay=RNA \
      --sample-col=sample \
      --celltype-col=AnnotationMain \
      --min-cells=20

Requires R, Matrix and SeuratObject installed in your approved environment. It must run with the correct Seurat object version. For multi-layer Seurat v5 assays, resolve the true raw counts layers and sample ordering **under scientific review** before export, rather than silently replacing counts with normalized data. Expect high RAM/storage needs; benchmark memory with approved subsets.

Audit files and output CSV both include donor-level pseudonymous metadata and must remain in approved controlled storage. This command will refuse to overwrite an existing file. It emits one explicit zero row for **every truly measured feature in every exported donor/cell-type group** and sets library_size to the full source-feature count total for that group. Groups with fewer than min-cells are written to local audit with excluded_low_cell_count; they are not included in the output.

Verify the audit report and run:

    python scripts/audit_cohort.py \
      --pseudobulk data/processed/donor_pseudobulk.csv \
      --out results/atlas_cohort_audit.csv

A passing structural audit is NOT evidence of independent donors, adequate power, consent, batch balance, true phenotype labels, or permission to distribute the files. A very small cell-type stratum may still fail the donor-level analysis. Counts from this adapter do not yet support causal interpretation.

## CI smoke test without human data

    Rscript tests/r/test_export_pseudobulk.R

The R test validates sparse Matrix aggregation, integer counts, total library-size accounting, actual zero rows, donor mapping, refusal to overwrite outputs, and failure on mismatched metadata. No Seurat installation or real dataset is required.

## Next scientific step

Audit GSE266265 object version/Seurat layer semantics, map donor and cell labels with an independent research specialist, verify sample-by-phase distribution, and compare output to an authors-reported published result. The current Python logCPM Welch test is only a computational baseline: a pre-specified covariate-adjusted count model and genuinely independent cohort are needed before making gene-level research claims. See issues #2–#6.
