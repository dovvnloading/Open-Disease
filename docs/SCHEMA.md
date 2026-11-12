# Input and output contracts (v0.1)

## Cell-level CSV (small datasets only)

One row per measured cell/gene count, with columns:

| Field | Meaning |
| --- | --- |
| study_id | Independent study identifier |
| donor_id | Pseudonymous within-study donor identifier; never a name or patient identifier |
| diagnosis | Exactly case or control |
| phase | Documented menstrual phase, consistent within donor |
| hormone_treatment | Exactly none in current untreated baseline |
| cell_type | Harmonized cell-type annotation |
| cell_id | Cell barcode/ID unique within donor |
| gene_id | Gene identifier, from one recorded genome/annotation release |
| count | Nonnegative integer count |

The aggregate command computes a donor-by-gene output with *explicit zeros*, along with library_size (sum of all included gene counts) and n_cells per donor/cell type. Only appropriate if the cell input contains the full count space used for normalization. For large atlases, aggregate in a suitable scientific tool first, verifying the full gene universe and library-size semantics.

## Donor-level pseudobulk CSV

Columns: study_id, donor_id, diagnosis, phase, hormone_treatment, cell_type, gene_id, counts, library_size, n_cells.

Each row is one gene in one donor/cell type. Counts is nonnegative; library_size is the **total unfiltered transcriptome count** for that donor/cell-type sample, not just the chosen candidate genes. n_cells is the number of included cells. The gene matrix must include one row for every assayed gene in every donor in the same study/phase/cell-type comparison. Never encode unmeasured as zero.

The code rejects duplicate donor/gene rows, inconsistent library sizes, incomplete within-stratum gene matrices, changed diagnosis within donor and insufficient groups. The schema alone is not evidence that data are consented for this use.

## Differential expression output

study_id, phase, cell_type, gene_id, n_case, n_control, effect_log2_cpm, p_value, status, method, q_value.

A negative effect means lower average log2(1+CPM) in cases. q is BH-adjusted **within** study/phase/cell type. Status other than ok is untestable. There is no global testing-family correction or covariate adjustment in v0.1.

## Genetics evidence CSV — manually curated, NOT automatically calculated

gene_id, evidence_type, support_score, source_id, safety_review.

Allowed evidence_type: fine_mapping, colocalization, integrated_multiomics, locus_only. A locus_only record cannot enter a ranked target table. support_score is a declared **heuristic** in [0,1]; do not interpret it as the probability a gene causes disease. source_id must trace to actual material and evidence method, including genome build and variant-to-gene mapping. safety_review is not_assessed, needs_review, or known_concern. A not_assessed flag does not mean safe.

The supplied demo genetics table is entirely fabricated and clearly marked SYNTHETIC_NOT_REAL_GENETICS.

## Validation output

Rank, frozen discovery information, validation_study, validation effect/q/donor counts, replicated Boolean and validation_status (replicated_exploratory, not_replicated or not_measured). Holdout results are never incorporated into rank sorting.

**Privacy:** even donor-level outputs may expose protected human subject information. All actual human-derived datasets and output folders must remain in approved local/controlled storage, not GitHub.
