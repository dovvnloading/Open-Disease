# Endometriosis Protocol v0.1 — pre-human-data working specification

Created 12 November 2026. This is a **working protocol**, not formal preregistration. Freeze cohort selections, contrasts, code revision and endpoints **before viewing independent real-world validation outcomes**. Changes must be versioned and disclosed. Today all executable experiments use fully synthetic data.

## Question and hypotheses

Question: Can disease-associated eutopic-endometrial genes show convergence of meaningful gene-level genetic evidence and disease-specific cell-state expression that replicates in an independent donor cohort?

H1 (unproven): Some cell-type-by-cycle-phase changes will have concordant direction and survive independent donor-cohort validation, once covariates are adequately addressed.

H0: Apparent signals are explained by donor selection, phase imbalance, technical batch, hormonal exposure, altered cell composition, or chance.

## Endpoint hierarchy

**Primary research endpoint** (not a clinical endpoint): freeze a ranked gene table using discovery data and independently curated variant-to-gene evidence; test direction and the prechosen FDR threshold in a genuinely independent donor study. Count assessed, replicated, unmeasured and conflicting candidates; publish all.

Secondary: assess reproducibility across cell types and cycle phases; track known normal-reproductive-tissue risks, fertility considerations, and eventual transfer to menstrual fluid if separately approved.

A clinical diagnostic AUROC is explicitly **not** a primary endpoint. Ten menstrual-fluid pilot donors cannot support such a claim.

## Units, cohorts and splits

- Each participant/donor is one independent biological unit per study; repeated samples and aliases must be identified, and no donor may appear in both discovery and validation.
- The currently verified atlas has 35 cases and 25 controls and includes different menstrual phases. The released code must not split its individual cells across train/test and call that donor-independent validation.
- An external cohort needs independently sampled donors and compatible assay annotations. A held-out subset of the same atlas is internal validation only.
- Hormone-treatment exclusion is enforced for the first baseline (hormone_treatment=none). Cross-treatment effects need their own protocol.
- Report phase definitions, age, batch, infertility history if permitted, tissue site, disease confirmation/stage, ancestry and missingness. Confounder-adjusted models are required prior to scientific interpretation.

## Current v0.1 baseline

1. Validate study/donor IDs and metadata consistency.
2. Aggregate cell-level gene counts to donor × cell type; explicit zero counts required.
3. Use full-transcriptome library size for counts-per-million normalization, never the library size of a preselected gene list unless the source truly includes the whole library.
4. Separate by study, menstrual phase and cell type.
5. Exploratory Welch test on log2(1+CPM), donor as unit; BH correction within each declared study/phase/cell-type family. Both groups require at least three donors for a computable test; this technical minimum does **not** imply adequate scientific power.
6. A target can be provisionally ranked only when discovery q <= 0.10 and supported by **curated gene-level** fine-mapping, colocalization or integrated evidence. An isolated GWAS locus is not adequate.
7. The transparent score combines curated support, absolute effect and within-family q-value. Scores are **not** probabilities and should not be compared across different input assessments without calibration.
8. Lock top candidates and test them in a separately labeled study. Replicated-exploratory means same effect direction and holdout q <= 0.10; it does not mean confirmed therapeutic action.
9. Report complete source trail, sample counts, effect direction, non-replications and unmeasured genes.

### Required improvements before running a real scientific study

- Predefine study/phase/cell-type families and global multiple-testing strategy, not only within-family BH.
- Benchmark against edgeR/DESeq2 or another validated donor-pseudobulk count model; adjust for available pre-specified confounders and assess dispersion.
- Independently validate phenotype labels and phase metadata, library sizes, batch handling, cell-type annotations, and possible patient overlap.
- Perform power analysis; avoid confident conclusions from small phase strata.
- Ensure genetics inputs are attributable to actual fine-mapping/colocalization analyses, not invented ratings, proximity to a GWAS lead SNP or synthetic demo scores.
- Test stable/negative-control genes, phenotype permutations at the donor level, potential expression-composition artifacts and multiplicity of attempted analyses.
- Have an independent computational biologist and clinician review interpretation and any public candidate claims.

## Failure / stop gates

**Stop before any human-data run** if study-level permissions, donor metadata, phase definitions, sample overlap, or feature completeness cannot be established. Stop candidate promotion for a result driven by technical batch, ambiguous variant-to-gene assignments, unsafe target biology or failure on the untouched holdout.

Publish negative results. A software test passing is not biological verification. A low p-value is not a cure.

## Excluded actions

No uncontrolled patient enrollment, specimen collection, patient genetic uploads, treatment suggestions, unsupervised lab work, human gene-editing proposals, dosing experiments or clinical deployment within this software program.
