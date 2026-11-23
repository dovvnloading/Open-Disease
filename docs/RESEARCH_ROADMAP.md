# Research roadmap: compute-first, evidence-first

## Current priorities

1. Test whether the evidence for higher MMP2 in endometriosis also supports higher **enzyme activity in peritoneal fluid**. A protein measurement alone cannot answer that question. Record how the conclusion changes when restricted to assays that can measure MMP2 activity in the presence of endogenous inhibitors. [MMP2 task](../research/tasks/peritoneal-mmp2-activity.json)
2. The [completed anatomical comparison](research-notes/OUTFLOW_OBSTRUCTION_ENDOMETRIOSIS_2026-11-23.md) quantified how many missed diagnoses would explain the two observed differences. The remaining question requires evidence on how consistently women in each anatomical group were examined for endometriosis. The published counts cannot establish that.

For the next proposal, identify what the primary sources already established, which observation the project can add, and which finding would change the decision. Check the study history before selecting evidence. Reuse needs a specific scientific justification, including for a different analysis of the same donors. [Task selection and review](AGENT_WORKFLOW.md#what-would-this-task-add)

## Completed expression screen

The first genome-wide external screen is complete. Among 12,044 shared genes and 1,049 human Reactome pathways, none cleared discovery multiplicity control and none of ten locked candidates validated in Rome. All five nominated pathways reversed direction. [Full screen, validation and phase diagnostic](research-notes/CROSS_COHORT_DISCOVERY_RESULTS_2026-11-14.md)

Four earlier [FOSB comparisons](../README.md) remain documented. The original UCSF comparison requires individual eligibility evidence.

## Workbench milestones

- [x] Baseline synthetic Python/R pipelines, source catalog, tests and CI.
- [x] Agent instructions, machine-readable task queue, resource limits and local run-provenance scaffold.
- [x] Audit the full GPL570 platform and reconcile the GSE51981 group and phase fields across all 148 array records. [Annotation results](research-notes/GPL570_ANNOTATION_RESULTS_2026-11-13.md) · [Cohort evidence and unresolved donor eligibility](research-notes/GSE51981_COHORT_AUDIT_2026-11-13.md)
- [x] Run the synthetic Python/R CI and test the implemented rejection gates. Keep further regression testing in routine CI; it is not a separate research milestone.

Exit: an agent can select a bounded task, reproduce or report a blocker, and leave a source-cited handoff with its commands and limitations.

## Original UCSF reproduction (blocked on donor eligibility)

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
- [x] Assess whether the public records establish a common measurement for the pooled Finnish FOSB row. The audit found insufficient provenance for a cross-chip comparison. [Decision and reopening criteria](research-notes/GSE141549_FOSB_CROSS_PLATFORM_DECISION_2026-11-16.md)
- [ ] Audit GWAS fine mapping and public target evidence without treating nearest genes or uncalibrated scores as causal confidence.
- [x] Freeze each external comparison before inspecting its target outcomes and retain inconclusive results alongside positive estimates.
- [ ] If public data cannot answer a question, disclose the limitation and change research direction rather than automatically recruit a collaborator.

## Follow-up evidence from expression and genetics

A [Melbourne expression-summary check](research-notes/MELBOURNE_FIVE_GENE_RESULTS_2026-11-14.md) tested the one nominated gene represented in its high-detection panel. The [Swiss atlas supplement audit](research-notes/GSE266265_SUMMARY_FEASIBILITY_2026-11-14.md) found significance-filtered cell-type tables. The Swiss [epithelial object audit](research-notes/GSE266265_EPITHELIAL_FEASIBILITY_2026-11-14.md) verified sample-to-batch links and four capture batches containing six affected and seven comparator samples under the authors’ strict clinical rules. Its corrected count layer leaves raw-UMI provenance unresolved. A donor-level raw-count study now needs a verified original RNA count source and a frozen cell-type contrast. Reusing a discovery donor, a second paper on the same arrays or a phase subset of an evaluated cohort cannot supply new replication. [Prior cohort reuse and literature](research-notes/CROSS_COHORT_PRIOR_EVIDENCE_2026-11-14.md)

The [Swiss glandular study](research-notes/SWISS_GLANDULAR_RESULTS_2026-11-14.md) now reports a prospectively defined corrected-expression comparison in those 13 donors. The positive combined score had exact p = 0.167. All five earlier whole-tissue reversals persisted after matching the genes.

The [completed Human Endometrial Cell Atlas audit](research-notes/HECA_GLANDULAR_FEASIBILITY_2026-11-14.md) located 13 donors across two contributing studies in the proliferative functionalis gland states. Verified cell and library annotations left only two independent Oxford case–control pairs sharing any qualifying gland cells within a capture. Huang had no mixed capture. Even the most permissive Oxford design has a two-sided p-value floor of 0.50. A higher-resolution computational test needs a different cohort or a prospectively justified cross-library design. [Cell and library census](../evidence/heca-glandular-2026-11-14/cell-library-support.json)

The [two additional proliferative whole-tissue comparisons](research-notes/EXTERNAL_BULK_PATHWAYS_RESULTS_2026-11-14.md) gave opposing combined estimates. A physical scan audit found extensive reuse across the U.S. archives. Matching 111 pathway genes across all five analyzed groups preserved the direction differences. The [subsequent locus audit](research-notes/GWAS_3P21_COLOCALIZATION_AUDIT_2026-11-16.md) examined what published colocalization at 3p21 could establish about inherited susceptibility.

## Possible future disease tracks

ALS (TDP-43 splicing), type 1 diabetes (islet/immune evidence), idiopathic pulmonary fibrosis (fibroblast states) and ME/CFS (cohort-reproducible phenotypes) remain **unstarted computational questions**, each needing its own protocol. We are not building a universal cure engine.
