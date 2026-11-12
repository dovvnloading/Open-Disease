# Rights-reviewed Seurat raw-count export to donor pseudobulk

**Research-only adapter; not tested on the real GSE266265 atlas.** The export algorithm has synthetic sparse-matrix tests; the large Seurat RDS reading path requires confirmation in a compatible, high-memory R environment by a qualified researcher.

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