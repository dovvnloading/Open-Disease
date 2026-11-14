# Research roadmap: compute-first, evidence-first

**As of November 2026:** Open-Disease has completed a real platform annotation audit and reconciled public cohort metadata. The original GSE51981 comparison remains blocked on individual eligibility. A prespecified comparison in a separately recruited surgical cohort has been executed. [GSE25628 FOSB result](research-notes/GSE25628_FOSB_RESULTS_2026-11-13.md) This is a small independent AI research effort. Its milestones are computational outputs, not fundraising, staffing or validation claims.

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

- [x] Audit recruitment, hormone status, phase, donor uniqueness, platform and run balance in Finnish EndometDB. [Eligibility evidence and remaining measurement question](research-notes/GSE141549_ELIGIBILITY_2026-11-13.md)
- [ ] [Resolve the EndometDB measurement (#26)](https://github.com/dovvnloading/Open-Disease/issues/26). The [measurement audit](research-notes/GSE141549_FOSB_MEASUREMENT_2026-11-13.md) verified the combined row label and the separate WG-6 expression workbook. A documented expression-to-clinical join or per-chip harmonization record would allow the donor-linked comparison to be reassessed.
- [ ] Audit GWAS fine mapping and public target evidence without treating nearest genes or uncalibrated scores as causal confidence.
- [ ] Freeze hypotheses before external comparisons; retain complete negative and inconclusive results.
- [ ] If public data cannot answer a question, disclose the limitation and change research direction rather than automatically recruit a collaborator.

## Conditional, not a current project phase

If a result eventually warrants wet-lab, human-subject, diagnostic, therapeutic or clinical claims, those activities require appropriate expertise, consent/oversight and independently verifiable experimental evidence. Do **not** initiate outreach, grant applications or a new lab program through agents; seek explicit operator authorization only when there is a concrete exceptional reason.

## Possible future disease tracks

ALS (TDP-43 splicing), type 1 diabetes (islet/immune evidence), idiopathic pulmonary fibrosis (fibroblast states) and ME/CFS (cohort-reproducible phenotypes) remain **unstarted computational questions**, each needing its own protocol. We are not building a universal cure engine.
