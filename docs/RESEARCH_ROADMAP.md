# Research roadmap: compute-first, evidence-first

**November 2026.** The first genome-wide external screen is complete. Among 12,044 shared genes and 1,049 human Reactome pathways, none cleared discovery multiplicity control and none of ten locked candidates validated in Rome. All five nominated pathways reversed direction. [Full screen, validation and phase diagnostic](research-notes/CROSS_COHORT_DISCOVERY_RESULTS_2026-11-14.md)

Four earlier [FOSB comparisons](../README.md) remain documented. The original UCSF comparison requires individual eligibility evidence.

## Now — make agents productive

- [x] Baseline synthetic Python/R pipelines, source catalog, tests and CI.
- [x] Agent instructions, machine-readable task queue, resource limits and local run-provenance scaffold.
- [x] Audit the full GPL570 platform and reconcile the GSE51981 group and phase fields across all 148 array records. [Annotation results](research-notes/GPL570_ANNOTATION_RESULTS_2026-11-13.md) · [Cohort evidence and unresolved donor eligibility](research-notes/GSE51981_COHORT_AUDIT_2026-11-13.md)
- [x] Run the synthetic Python/R CI and test the implemented rejection gates. Keep further regression testing in routine CI; it is not a separate research milestone.

Exit: an agent can select a bounded task, reproduce or report a blocker, and leave a source-cited handoff with its commands and limitations.

## Next — first real computational research result

- [ ] Resolve the study's source rights, independent sample identities, diagnosis categories and menstrual phases, **before** inspecting an outcome.
- [x] Pin the complete GSE51981 processed array and specify the FOSB probe, contrasts, phase rules, negative control and multiplicity. [Frozen comparison](research-notes/GSE51981_FOSB_ENDPOINT_2026-11-13.md) · [Full matrix provenance and blocked outcome](research-notes/GSE51981_REPRODUCTION_STATUS_2026-11-13.md)
- [ ] Run a real public-data reproduction using the existing opt-in GEO + limma path and publish a permissible **aggregate** account of positive, null or failed results.
- [ ] Compare findings with source-backed published evidence, record failures and distinguish reanalysis of an existing cohort from genuinely independent replication.

Exit: a credible, reproducible computation or a documented blocker, **not** a validated disease target. See `docs/INDEPENDENT_RESEARCH_PLAN.md` and the current [reproduction issue #18](https://github.com/dovvnloading/Open-Disease/issues/18). The ~39 GB Seurat atlas, new ML models and novel drugs are **not** prerequisites.

## External cohort comparisons

- [x] Execute the single-probe GSE25628 proliferative-phase case–control comparison. The estimate was positive and its interval included zero. Controls had benign surgical indications and individual recent treatment was unknown. [Report](research-notes/GSE25628_FOSB_RESULTS_2026-11-13.md)

- [x] Resolve Finnish native assay identity for all 38 eligible donors and run the prospectively defined physical-probe test in seven shared chip/run/phase groups. The percentile difference was positive with an inconclusive exact test. [Chip census](research-notes/GSE141549_CHIP_CENSUS_2026-11-14.md) · [Result](research-notes/GSE141549_WG6_FOSB_RESULTS_2026-11-14.md)
- [x] Test the secretory-phase contrast in GSE120103. Verify the original scanner headers and retain the complete separation of clinical groups across slides as an attribution limit. [Result](research-notes/GSE120103_FOSB_RESULTS_2026-11-14.md)
- [x] Test a separately specified RNA-seq transcript contrast in eight proliferative biopsies. Retain the change of sign on participant omission. [Result](research-notes/GSE315857_FOSB_RESULTS_2026-11-14.md)
- [ ] Resolve the authors' cross-chip derivation of the pooled normalized Finnish FOSB row. This question is separate from the completed native WG-6 analysis. [Issue #26](https://github.com/dovvnloading/Open-Disease/issues/26)
- [ ] Audit GWAS fine mapping and public target evidence without treating nearest genes or uncalibrated scores as causal confidence.
- [x] Freeze each external comparison before inspecting its target outcomes and retain inconclusive results alongside positive estimates.
- [ ] If public data cannot answer a question, disclose the limitation and change research direction rather than automatically recruit a collaborator.

## Next computational question

A [Melbourne expression-summary check](research-notes/MELBOURNE_FIVE_GENE_RESULTS_2026-11-14.md) tested the one nominated gene represented in its high-detection panel. The [Swiss atlas supplement audit](research-notes/GSE266265_SUMMARY_FEASIBILITY_2026-11-14.md) found significance-filtered cell-type tables. The Swiss [epithelial object audit](research-notes/GSE266265_EPITHELIAL_FEASIBILITY_2026-11-14.md) verified sample-to-batch links and four capture batches containing six affected and seven comparator samples under the authors’ strict clinical rules. Its corrected count layer leaves raw-UMI provenance unresolved. A donor-level raw-count study now needs a verified original RNA count source and a frozen cell-type contrast. Reusing a discovery donor, a second paper on the same arrays or a phase subset of an evaluated cohort cannot supply new replication. [Prior cohort reuse and literature](research-notes/CROSS_COHORT_PRIOR_EVIDENCE_2026-11-14.md)

The [Swiss glandular study](research-notes/SWISS_GLANDULAR_RESULTS_2026-11-14.md) now reports a prospectively defined corrected-expression comparison in those 13 donors. The positive combined score had exact p = 0.167. All five earlier whole-tissue reversals persisted after matching the genes.

The [completed Human Endometrial Cell Atlas audit](research-notes/HECA_GLANDULAR_FEASIBILITY_2026-11-14.md) located 13 donors across two contributing studies in the proliferative functionalis gland states. Verified cell and library annotations left only two independent Oxford case–control pairs sharing any qualifying gland cells within a capture. Huang had no mixed capture. Even the most permissive Oxford design has a two-sided p-value floor of 0.50. A higher-resolution computational test needs a different cohort or a prospectively justified cross-library design. [Cell and library census](../evidence/heca-glandular-2026-11-14/cell-library-support.json)

The [two additional proliferative whole-tissue comparisons](research-notes/EXTERNAL_BULK_PATHWAYS_RESULTS_2026-11-14.md) gave opposing combined estimates. A physical scan audit found extensive reuse across the U.S. archives. Matching 111 pathway genes across all five analyzed groups preserved the direction differences. A next computational question is whether published endometriosis GWAS–eQTL colocalization implicates any member of the frozen pathways, with tissue, locus and study overlap documented. Such evidence addresses inherited susceptibility, a different question from a general expression shift.

## Possible future disease tracks

ALS (TDP-43 splicing), type 1 diabetes (islet/immune evidence), idiopathic pulmonary fibrosis (fibroblast states) and ME/CFS (cohort-reproducible phenotypes) remain **unstarted computational questions**, each needing its own protocol. We are not building a universal cure engine.
