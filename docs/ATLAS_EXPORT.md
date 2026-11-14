# Rights-reviewed Seurat raw-count export to donor pseudobulk

**Research-only adapter.** Synthetic tests verify the sparse export. The complete 38.9 GB atlas still requires a suitable high-memory R environment. The [public summary feasibility audit](research-notes/GSE266265_SUMMARY_FEASIBILITY_2026-11-14.md) and [epithelial derivative audit](research-notes/GSE266265_EPITHELIAL_FEASIBILITY_2026-11-14.md) document which smaller releases retain donor information, and the source of their count layers.

## Why

The 2026 atlas has 466,371 cells from only **60 donors**, stored in a processed ~38.9 GB Seurat RDS. Our scientific units are donors. No cell should be counted as a separate participant.

This adapter reads an explicitly named **RNA counts layer**, groups by verified donor and cell type, and writes every assayed gene's raw summed count, including zeros. It does not calculate significant genes or perform clinical analysis.

## Before running: rights and phenotype audit

- Confirm the dataset's access/reuse terms, study version, biological sample identities, case definition, cycle phase, and hormone exposure. Document consent/ethics and any controlled-access restrictions.
- Use a **secured local/high-memory computer**. RDS on disk is about 39 GB; peak RAM can be much higher. No automatic downloads or cloud upload are performed here.
- Review Seurat metadata **within the approved environment** and construct a separate CSV mapping with *exactly one row per cell*, with columns:

```
cell_id,study_id,donor_id,diagnosis,phase,hormone_treatment,cell_type
```

- `cell_id` must exactly match a column name from the **raw RNA count matrix**. A mismatch or duplicate is a hard error.
- `diagnosis` must be **case** or **control**, explicitly harmonized from trusted clinical metadata; do not infer disease label from gene expression or use filenames as an undocumented classifier.
- `phase` is an audited menstrual phase label; values must be consistent within each donor. A missing phase is a hard error.
- `hormone_treatment` must be `none` for baseline v0, backed by an audited study inclusion rule. The program will not infer absence from a missing field.
- `donor_id` must group all technical libraries/cells from the same participant; verify cross-sample aliases. This is sensitive metadata and must never be committed.
- `cell_type` must have an explicit biological annotation/ontology. The script does not invent cell types or collapse them automatically.

The source lab's published analysis repository has **no formal source-code license**, so do not copy its implementation into this repository without separate permission.

## Run (on approved research infrastructure)

```bash
Rscript scripts/seurat_to_pseudobulk.R \
  --rds /approved/source/GSE266265_EndoAtlas.rds \
  --metadata /approved/derived/reviewed_cell_mapping.csv \
  --out /approved/derived/pseudobulk_raw_counts.csv

python -m open_disease sha256 /approved/derived/pseudobulk_raw_counts.csv

python scripts/audit_cohort.py \
  --pseudobulk /approved/derived/pseudobulk_raw_counts.csv \
  --out /approved/derived/cohort_audit.csv
```

Outputs must remain under approved access controls. Do not upload donor-level counts, donor labels or CSVs into GitHub or CI.

## How to review metadata without assumptions

Inside a rights-cleared R session:

```r
obj <- readRDS("/approved/source/GSE266265_EndoAtlas.rds")
colnames(obj[[]])    # inspect available metadata fields only
head(obj[[]], 3)      # sensitive; view locally, do NOT paste into public issues
```

Choose mapping fields only after reading their definitions. If metadata have several technical libraries per donor, harmonize them **before** export. A local mapping CSV can be generated using base R, but we deliberately do not supply a hard-coded upstream field translation because it has not been verified.

## Output and failure behavior

For each study/donor/cell-type group, the exporter emits one row for **every gene** in the source RNA count matrix:

- `counts` = sparse column-group sum of raw UMI counts.
- `library_size` = sum of **all source genes**, not just selected targets.
- `n_cells` = number of distinct matrix cell columns in the group.
- `diagnosis`, `phase`, `hormone_treatment` = reviewed, consistent donor annotations.

It refuses fractional counts, duplicate identifiers, missing mapping cells, inconsistent donor labels/phases, non-untreated cohorts, zero libraries, unsupported non-sparse/normalized assay layers and overwrite of existing output.

There is no automatic `JoinLayers` or assay conversion. If a Seurat v5 object has split/ambiguous counts layers, resolve the layer semantics manually with the data custodian. An "integrated" or log-normalized assay must not be described as raw counts.

## Synthetic test

```bash
Rscript scripts/test_seurat_to_pseudobulk.R
```

This tests the sparse aggregation with small **fake** data; it does not install SeuratObject or read any human data. The actual RDS/SeuratObject version and memory behavior remain open validation tasks. Next: independently rerun the original published analysis with an appropriate counts-aware method and covariate plan, then compare to this export.

---

## Alternative: sample-to-donor mapping interface and exclusion audit

The original export approach above expects a privately reviewed mapping row for every cell. A complementary offline adapter accepts a **one-row-per-sample** donor mapping joined with the chosen Seurat sample and cell-type fields. This can avoid creating a second large per-cell mapping CSV when the data custodians can certify the sample metadata. It does not infer donor or diagnosis labels. Both paths must produce the same gene-by-donor counts on equivalent inputs; tests/r/test_cross_exporters.R confirms this using artificial count matrices.

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

Audit GSE266265 object version/Seurat layer semantics, verify the donor-to-cell metadata with an independently implemented check, verify sample-by-phase distribution, and compare output to an authors-reported published result. The current Python logCPM Welch test is only a computational baseline: a pre-specified covariate-adjusted count model and genuinely independent cohort are needed before making gene-level research claims. See issues #2–#6.
