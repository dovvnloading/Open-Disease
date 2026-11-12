# edgeR v4 donor-level quasi-likelihood model — experimental implementation

**NOT FOR CLINICAL USE.** This implementation has only been exercised with simulated SYN### genes and fake donor data. The existence of an edgeR model is not evidence of an endometriosis mechanism, reliable biomarker, or cure. It has not been approved by a domain biostatistician for the actual 2026 atlas.

## Scientific rationale

The initial Welch test on log2-counts-per-million is a **software baseline only**. For future adequately powered donor-level research, the project now implements a count-aware negative-binomial quasi-likelihood (QL) workflow using the established edgeR v4 package. It is not a new mathematical method created by Open-Disease.

Public method references:

- Official Bioconductor edgeR: https://bioconductor.org/packages/release/bioc/html/edgeR.html
- 2025 edgeR v4 methods paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11754124/
- Current package manual and v4 QL workflow: use R command browseVignettes("edgeR") after installation.

## Required data

Supply the **complete raw integer count** donor-by-cell-type, donor-by-gene CSV produced by a rights-cleared atlas export. All genes in the assayed universe must have explicit rows even when observed count is zero. Every donor-cell-type raw library_size must exactly equal the sum over the exported genes, avoiding the misleading normalization of candidate-only gene subsets.

Minimum structural software gates are 4 independent cases and 4 controls, at least 20 raw cell observations per donor/cell-type, at least 20 genes and at least four residual degrees of freedom; these deliberately low gates permit computational tests, **not** biologically powered inference.

The human dataset must first pass rights, independent donor ID, anatomical-site, phase, exogenous hormone, assay and cell-type provenance review. Do not automatically infer hormone non-exposure from a missing field; the code requires documented untreated status and valid phase metadata.

### Allowed adjustment designs

The default is not to run at all without an explicit adjustment choice. Available choices:

| --adjust | Model | Required additional covariate CSV columns |
| --- | --- | --- |
| none | diagnosis | none; exploratory, unadjusted only |
| batch | diagnosis + batch | study_id, donor_id, batch |
| age | diagnosis + age_years | study_id, donor_id, age_years |
| batch+age | diagnosis + batch + age_years | study_id, donor_id, batch, age_years |

The model has control as the diagnosis reference, so a positive logFC means higher expression in cases. All covariates must be complete and uniquely matched to the donor ID. Non-estimable models, such as all cases processed in batch A and all controls in batch B, fail **rather than silently produce disease associations**. If covariates are missing or the contrast is not estimable, pre-specify and document the alternative—never secretly drop a confounder after inspecting results.

These limited covariate choices do not imply that stage, sequencing site, ancestry, fertility, surgery history, cell composition, medication or other bias are controlled. Phase is stratified; some phases may be underpowered. Do not interpret an unadjusted result as robust discovery.

## Run

Install an R version compatible with a stable Bioconductor release and edgeR v4 or newer, then run locally in an authorized high-memory data environment:

    Rscript scripts/edger_pseudobulk_ql.R \
      --pseudobulk /approved/derived/pseudobulk.csv \
      --study GSE266265 --phase proliferative --cell-type epithelial \
      --adjust batch+age --covariates /approved/derived/covariates.csv \
      --out /approved/derived/differential_expression.csv

A synthetic-only software demo runs with:

    bash scripts/run_demo.sh
    Rscript scripts/test_edger_pseudobulk_ql.R

The pipeline: donor-level raw count matrix -> filterByExpr using the model design -> TMM normalization -> edgeR glmQLFit with modern dispersion estimation -> glmQLFTest testing the case-control coefficient -> Benjamini-Hochberg correction within the **preselected single study × phase × cell-type stratum**.

It produces every assayed gene with its tested or filtered status, the edgeR QL effect/p-values/FDR where appropriate, donor sample sizes, adjustment design, residual degrees of freedom and package version. It refuses duplicate donor/gene pairs, incomplete matrices, fractional counts, inconsistent donor metadata, library-size mismatch, missing covariates, rank-deficient designs and insufficient software sample gates. It never overwrites an existing result file.

## Why this does NOT yet mean biological validation

1. There is **no actual atlas expression analysis** in this project. Current CI only uses simulated data.
2. Statistical tests across dozens of cell types and phases require multiplicity control across families, not just within one stratum.
3. A minimum of four donors/group is a technical stop gate, not an adequate power calculation; real phase-specific groups may be much smaller.
4. The same human cohort cannot be both novel discovery and independent external confirmation. Identifying a comparable, independently sampled holdout remains an open blocker.
5. A log-fold change is a statistical association—not causal gene action, a treatable target or evidence of clinical benefit.
6. The edgeR dependency is GPL-licensed by its maintainers and is installed separately. This repository contains only an independently written script calling its APIs; our MIT license does **not** replace edgeR's license.

### Scientific promotion gate

A qualified biostatistician and endometriosis domain researcher must review sample eligibility, covariates, power, residual/model diagnostics, data rights, multiplicity and realistic biological endpoints. Only then can this script be used as a component in a frozen prospective hypothesis evaluation. Record every exclusion and failed result.
