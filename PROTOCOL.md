# Endometriosis research protocol — draft v0.1

**Date:** 2026-11-12  
**Status:** PUBLIC DRAFT; **not formally preregistered**, not a published protocol and not validated on human data. Freeze a versioned protocol before accessing outcome labels for the first full human analysis. Changes after that are exploratory and must be explicitly labeled.

## Question and scope

Which endometriosis-associated genes are supported **simultaneously** by credible human genetic mapping, interpretable donor-level cell-specific expression changes, and comparable independent human cohorts, while presenting tractable experimental and reproductive-safety questions?

We are looking for **hypotheses for wet-lab testing**. Distinct endpoints: (a) mechanistic target discovery, (b) biomarker development, (c) clinically effective therapy. A good (a) cannot be represented as achievement of (b) or (c).

## Hypotheses, including null

H1 (exploratory): a nonempty group of genes with fine-mapping/colocalization support will show cell-type-specific case/control expression differences in untreated endometrial tissue after accounting for menstrual cycle phase.

H2 (prospective replication): at least one frozen candidate may show a consistent direction and statistically supported effect in an **independent, tissue-matched** cohort. A candidate that fails is a valid negative outcome.

H0: observed signals arise from chance, patient selection, menstrual-cycle effects, platform/batch effects, tissue composition, prior treatment, ancestry, or other confounders, and do not replicate when these are addressed.

No hypothesis should be evaluated with a pipeline tuned on validation outcomes.

## Evidence unit, cohort boundaries, and data provenance

- Independent biological unit = donor/participant, not cell, nucleus, biopsy, library or aliquot. Multiple libraries from one donor must share the donor key.
- The primary atlas comprises **35 cases / 25 controls**; 466,371 cells are technical/biological observations nested inside those 60 donors.
- Primary discovery tissue = **eutopic endometrium**, untreated with exogenous hormones under the atlas's published exclusion definition.
- Preserve study ID, stable donor key, diagnosis, cell type, sample and anatomical site, menstrual-cycle phase, sequencing batch, known exogenous hormone exposure, tissue preparation, inclusion/exclusion notes and count matrix provenance.
- Phase is mandatory for the v0 baseline. Missing or conflicting phase metadata prevents use until an explicit alternative protocol is written.
- Any follow-up study must be audited for **anatomical tissue**, phase, hormones, donor identity overlap, platform and case definition. Comparing lesions against healthy uterine tissue measures multiple differences and is NOT straightforward external replication of case/control endometrial tissue.
- Independent does not merely mean a different GEO accession. Confirm study/cohort overlap and exact donor identity where possible; otherwise classify independence as unverified.
- Record file SHA-256, original URL, retrieval date, sample exclusions and version of parsing/export code locally (do not upload participant-level artifacts).

## Stages and frozen decision rules

### Stage A: structural preflight

Input: a long table with columns `study_id, donor_id, diagnosis, phase, hormone_treatment, cell_type, gene_id, counts, library_size, n_cells`. The schema rejects duplicates, missing keys, negative/noninteger counts, donor metadata inconsistencies and invalid library sizes. An explicit zero row is required when a gene is measured but unexpressed.

The default software baseline requires at least **three donors per class in each study × phase × cell type stratum**, >=20 cells per donor/cell type and >=20 genes for technical smoke tests. **These are minimal software guards, NOT sample-size/power sufficiency thresholds for biological inference.** A statistician must define meaningful group sizes and power before any scientific result is interpreted.

The preflight tool does not automatically establish covariate balance, consent/data rights, biological independence, gene universe completeness, batch correction or cell annotation quality. These require documented human review.

### Stage B: reference reproduction (before hypothesis generation)

Replicate at least one published cohort result using its original documented gene universe/filters. Perform stringent identity and QC checks and preserve negative results. Compare counts-aware pseudobulk methods (such as edgeR/DESeq2) to the shipping Welch/log-CPM *software baseline*. Do not promote baseline q-values to confirmatory inference.

Check phase, case severity, diagnosis modality, hormone treatment, age, ancestry where available, QC batch, donor-level sampling biases, cell type abundance, cell-type-specific gene detection, and tissue source. Evaluate confounder imbalance rather than assume it absent.

### Stage C: gene mapping and candidate nomination

An endometriosis GWAS risk locus is not automatically a disease gene. A gene receives genetic support only with documented variant-to-gene evidence (e.g. robust fine mapping, colocalization, independent gene-expression/transcript evidence). Nearest-gene-only or unannotated locus evidence is listed but **not eligible** for ranking by the current software.

For each gene: genetic source and mapping method, tissue, transcript version, effect size/sign, source provenance, uncertainty, GWAS sample population and overlap assumptions, and potential adverse effects on normal reproductive tissue. Numerical support scores in the proof-of-concept software are transparently **subjective heuristics** with no probabilistic calibration. Freeze any ranking formula before independent testing; sensitivity-test its weights rather than choose weights using holdout outcomes.

### Stage D: independent validation

Lock the discovery study, genes, cell types, phase strata, inclusion/exclusion, signal direction and multiple-testing family before reading the independent study's labels and outcomes.

The shipping function checks direction and adjusted q-value <=0.10 within each holdout stratum, with sufficient data/status. This is only an **exploratory software convention**, not proof of biological validation. For inferential claims require formal power assessment, effect size confidence intervals, a well-defined multiplicity correction accounting for tested strata, and independent biological/statistical review.

Primary failure outcome: no candidate remains robust and testable. Report null results, target exclusions, potential confounds and limitations.

### Stage E: laboratory handoff

**Not a clinical trial.** If a target survives independent evaluation, qualified experimental collaborators must confirm relevant expression and functional effects in appropriate tissues/models, test effects on lesions, fertility and normal tissue, and examine toxicity/off-target and immune consequences. Positive cell experiments alone do not demonstrate clinical cure.

## Analysis pitfalls and required negative controls

1. **Pseudoreplication:** never use cell counts as independent n; donor-level sampling and grouped cross-validation.
2. **Donor overlap:** patient/library technical replicates cannot appear on both sides of a split; leave-one-donor-out is not external validation.
3. **Phase/hormone imbalance:** calculate within-phase effects and show exclusions; avoid extrapolating across menstrual stages.
4. **Tissue/assay confounding:** eutopic tissue ≠ ectopic lesion ≠ menstrual fluid. A different platform may support directional triangulation, not identical-effect replication.
5. **Gene selection leakage:** do not use holdout to select genes, filters, scores, normalization, thresholds or report examples.
6. **False-positive burden:** pre-specify all tests and adjustment family; show tested-gene counts, model limitations, complete effect tables and negative comparisons.
7. **Sample size:** pilot cohorts with 5 cases/5 controls cannot establish diagnosis or reproducible biomarkers; multiple RNA libraries per donor do not increase independent sample size.
8. **Molecular causality:** association with expression or genetic risk does not show that inhibiting/activating a gene reverses disease. Require directional perturbation and safety evidence with lab partners.
9. **Reproducibility:** seed stochastic steps, version packages, hash input files and archive methodology (where distribution permitted).
10. **Safety and equity:** prioritize analyses considerate of symptom burden, fertility, excluded groups, ancestry, endometriosis phenotype and accessibility; record important representation gaps.

## Deliverables and gates

- **Gate 0 (current):** code/tests/CI operate on simulated data only; science status = not started.
- **Gate 1:** data rights, provenance, donor cohort reconciliation and analysis plan independently checked; first extraction produces validated *local* pseudobulk with checksums.
- **Gate 2:** reproduce a published result using appropriate count-based methods; publish full failures and QA report.
- **Gate 3:** independently tested, frozen candidate shortlist with case/control donor counts and uncertainty, versus baseline/null models; avoid causal claims.
- **Gate 4:** documented interest from domain scientist, patient representative and qualified laboratory; experimental plan includes negative controls and normal reproductive tissue safety.

**Stop or change scope** if permission is unclear, confounding dominates, independent cohorts are not comparable, overlap cannot be ruled out, sample size fails pre-specified power needs, candidates fail replication, or proposed experiments cannot address reproductive safety.

## Protocol amendment policy

Amend with a dated Git commit documenting what changed, why, whether new outcomes were already seen, and which results become exploratory. Do not rewrite history to make post-hoc choices appear prospective. The first public commit of this draft is **not** a registered trial/protocol.
