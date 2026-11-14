# What the Swiss epithelial trajectory file supports

**Source and assay assessment, November 14, 2026.** The small GSE266265 epithelial derivative preserves donor-to-sample clinical and capture-batch information through an exact join to Supplementary Data 1. Its matrices follow a depth-corrected SCT workflow. We could not establish original RNA UMI counts in this release. Restricting the authors' proliferative differential-expression eligibility to retained patients in capture batches with both study groups leaves **six cases and seven controls**. This audit produced no disease-associated expression estimates. [Source study](https://doi.org/10.1038/s41467-026-73020-4) · [GEO release](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265) · [Issue #35](https://github.com/dovvnloading/Open-Disease/issues/35)

## Physical measurement and selection

| Property of the deposited object | Verified finding |
| --- | --- |
| Monocle object | One `cell_data_set`, 29,034 distinct feature labels and 20,549 distinct cell labels. Both matrices have the same order as the native row and cell records. |
| Assays | Two compressed-sparse-column matrices, `counts` and `logcounts`. All stored values are finite and nonnegative. Every stored count is integer-valued and the log layer has fractional values. There are zero alternative experiments. |
| Feature naming | 6,638 labels match the human Ensembl gene ID pattern. The other 22,396 have symbol-compatible syntax. The rows contain exact symbol labels for ABO, KIF3C, ZNF835 and GLTP. F11 lacks an exact symbol label. Resolving alternate identifiers requires a versioned gene crosswalk. |
| Cell selection | All 20,549 cells have one of seven refined epithelial labels. One lacks the broad annotation. Of the 58 retained samples, 53 have at least 100 remaining cells. |

The source workflow selected epithelial cells, attempted to retain samples with more than 99 epithelial cells, then ran SCT v2 with mitochondrial percentage and sample specified for regression and split the data by 10x capture batch. The authors' integration script uses SCT and the trajectory script sets the active assay to SCT immediately before converting to a Monocle data set. Seurat documents the default SCT count layer as depth-corrected UMI counts; the source call does not override that default. Numeric integrality therefore cannot establish original library depth. The published Muscat differential-expression analysis used the RNA assay with capture batch as a covariate. The small object supplies no independently identifiable raw-RNA assay for rebuilding that count model. [Authors' normalization](https://github.com/DuempelmannLea/endometriosis_endometrium_scRNA_atlas/blob/e00da6a525c7300022f0e1d66fdb6199e6910efc/01_Atlas_generation/Integration_epithelial_only/01b_run_normalization_split.R) · [Integration](https://github.com/DuempelmannLea/endometriosis_endometrium_scRNA_atlas/blob/e00da6a525c7300022f0e1d66fdb6199e6910efc/01_Atlas_generation/Integration_epithelial_only/02_run_seuratv3.R) · [Monocle conversion](https://github.com/DuempelmannLea/endometriosis_endometrium_scRNA_atlas/blob/e00da6a525c7300022f0e1d66fdb6199e6910efc/04_Time_trajectory_analysis/trajectoryAnalysis_monocle3.r) · [Seurat field definitions](https://satijalab.org/seurat/reference/sctransform)

The paper's trajectory Methods specify removal of ciliated, mesothelial and MUC5B-positive cells. The generic source script's cluster-subset block is guarded by a null check and does not provide a populated exclusion list. The preceding >99-cells criterion refers to preintegration selection. Exact final per-sample inclusion probabilities and a per-subtype sampling cap are undocumented in these scripts. The object can describe its recorded epithelial states, but lacks the full epithelial and whole-tissue denominators required to estimate depletion or enrichment of all such cell types. [Study, trajectory Methods](https://doi.org/10.1038/s41467-026-73020-4) · [Pinned workflow](https://github.com/DuempelmannLea/endometriosis_endometrium_scRNA_atlas/tree/e00da6a525c7300022f0e1d66fdb6199e6910efc)

## Clinical design

The native cell metadata contain sample, disease status and grade, cycle and phase, epithelial annotations and cell keys. There is no separate native donor or batch field. An exact match from the native sample label to the authors’ `Patient sample ID` recovered **58 of 58** samples among the 67 unique IDs in Supplementary Data 1. Disease and fine menstrual phase agreed for all 58; the 58 sample labels individually have consistent clinical group, grade and phase. The linked design has ten capture batches and two sequencing batches. Only the source-defined patient sample keys were joined. [Supplementary Data 1](https://pmc-oa-opendata.s3.amazonaws.com/PMC13385355.1/41467_2026_73020_MOESM3_ESM.xlsx)

| Native phase in the derivative | Cases | Controls |
| --- | ---: | ---: |
| Proliferative | 13 | 13 |
| Periovulatory | 13 | 5 |
| Early secretory | 2 | 2 |
| Mid-secretory | 5 | 2 |
| Late secretory | 0 | 3 |
| Total | 33 | 25 |

| Proliferative design or restriction | Cases | Controls | Shared capture batches |
| --- | ---: | ---: | ---: |
| All retained proliferative samples | 13 | 13 | — |
| Within capture batches containing both groups | 9 | 10 | 5 |
| Source DEG inclusion flag in the full clinical supplement | 12 | 11 | — |
| Source DEG flag among the retained samples | 10 | 11 | — |
| Source DEG flag and within-batch case–control overlap in the derivative | **6** | **7** | **4** |

The donor is the independent unit for these design counts. The 20,549 cells do not enlarge the number of independent patients. Any future analysis of corrected expression would need to specify a new endpoint, the eligible donors and cell subtypes, the capture-batch comparison, and treatment of genes missing under literal symbol mapping. A proportion computed from retained cells would describe the trajectory subset.

## Reproduce

The local audit uses base R serialization and direct attribute inspection with default R packages disabled. It checks feature/cell alignment and scans sparse numeric arrays in bounded chunks solely for finiteness, nonnegativity and integrality. It checks clinical consistency, performs an exact sample join, and emits only aggregate design and storage information. The optional Python dependency for the clinical workbook is openpyxl. Inputs and a temporary key projection stay in an ignored private directory; the projection is removed after each run. The separate study summary workbook is documented in the [public summary assessment](GSE266265_SUMMARY_FEASIBILITY_2026-11-14.md).

| Input | Bytes | SHA-256 |
| --- | ---: | --- |
| [Epithelial Monocle RDS](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE266nnn/GSE266265/suppl/GSE266265_epithelial_cells_monocle3.rds) | 199,996,542 | `50715d0e2f5a09b8f4075668569a75baff9ce24bf8e62503f259858bb1c7d683` |
| [Supplementary Data 1 XLSX](https://pmc-oa-opendata.s3.amazonaws.com/PMC13385355.1/41467_2026_73020_MOESM3_ESM.xlsx) | 21,019 | `21d7996ad209662a550fceb9694b3e9b70eb1a7e9cb44cce26275a60a10891dd` |

Download the two inputs above into `data/private/gse266265/` using their published filenames. From the repository root, with the project and openpyxl installed, run:

```sh
d=data/private/gse266265
python scripts/audit_gse266265_epithelial.py \
  "$d/GSE266265_epithelial_cells_monocle3.rds" \
  "$d/41467_2026_73020_MOESM3_ESM.xlsx"
```

The audit used Python 3.12.14, openpyxl 3.1.5 and base R 4.5.0. The RDS was retrieved once in 22.44 seconds. The complete pinned RDS and workbook check ran in 9.49 seconds with peak RSS 1,570,352 KiB under a 2,000,000 KiB address-space limit. The plan allowed 230 MB of downloads and 2 GB of RAM; no random procedure was used. [Python entry point](../../scripts/audit_gse266265_epithelial.py) · [R inspector](../../scripts/r/audit_gse266265_epithelial.R) · [Synthetic integrity tests](../../tests/test_audit_gse266265_epithelial.py)
