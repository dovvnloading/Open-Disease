# Research roadmap and decision gates

Version 0.1, November 2026. This is a plan, not a claim that scientific milestones have been achieved.

## Milestone 0 — build an honest research foundation

- [x] Preserve Open-Disease's MIT license and multi-disease mission.
- [x] Add a working donor-level pseudobulk prototype and baseline analysis.
- [x] Add curated-evidence ranking separated from holdout validation.
- [x] Add deterministic synthetic data and automated regression tests.
- [x] Document the actual studies, access restrictions, patent interests and science limitations.
- [ ] Verify CI is green and have an independent scientific reviewer assess protocol.

## Milestone 1 — make genuine data use permissible and reproducible

- [ ] Review GSE266265 processed object and its license, metadata dictionary, cell counts, sample IDs and menstrual phase categories.
- [ ] Verify Zenodo GWAS file layout and permitted uses; design reproducible gene-mapping schema with genomic build and provenance.
- [ ] Freeze the first real discovery/validation comparison and define independent donor cohorts and statistical power.
- [x] Implement a synthetic-tested sparse raw-count aggregator and explicit-metadata Seurat adapter.
- [ ] Validate the adapter's SeuratObject input behavior, high-memory performance and source field mappings on the actual rights-cleared atlas.
- [ ] Produce checksums, source inventory and an independently reproduced published baseline result.
- **Stop gate:** no actual human-data analysis if consent/terms or donor independence cannot be established.

## Milestone 2 — statistically defensible discovery

- [ ] Predefine multiple testing families, exact phenotypes, covariates and tissue contrasts.
- [ ] Fit validated count-based pseudobulk model (edgeR/DESeq2) and benchmark against v0.1 baseline.
- [ ] Perform phase-, treatment-, age-, batch- and cell-composition sensitivity analyses where permitted.
- [ ] Integrate actual fine-mapping/colocalization evidence rather than heuristic invented scores.
- [ ] Apply donor-level permutations, negative controls and leakage audit.
- **Stop gate:** no ranked therapeutic candidates promoted without source-level genetic evidence and credible sensitivity results.

## Milestone 3 — untouched external replication and experimental handoff

- [ ] Secure an actually independent eligible cohort; verify no donor overlap.
- [ ] Freeze the ranked list and evaluate all candidates without reselecting on holdout.
- [ ] Publish complete outcomes, nonreplications and uncertainty.
- [ ] Review normal reproductive-tissue and fertility safety with relevant experts.
- [ ] Partner with a gynecologist/scientist, biostatistician and organoid or tissue laboratory for appropriate functional experiments.
- **Stop gate:** no claims of disease reversal or safe therapy without qualified biological and eventual clinical evidence.

## Later: expand architecture to the remaining four diseases

ALS, type 1 diabetes, pulmonary fibrosis and ME/CFS remain research ideas in CHARTER.md. Do not equate their distinct data, biomarkers or causal biology with endometriosis. New programs require their own disease-specific protocols, controlled datasets, and domain review.

## What would count as progress

A real external result, a reliable negative finding, an improved computational benchmark or a laboratory-accepted testable hypothesis. Number of model parameters, GitHub stars or synthetic classifier performance are not scientific endpoints.
