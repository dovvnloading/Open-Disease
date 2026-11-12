# Research roadmap: five diseases, one disciplined first program

**Status, 2026-11-12:** the endometriosis research software foundation is under development; test data are simulated. The other disease tracks are not yet under active experimental analysis. No claims of cure or validated targets.

| Disease | First credible computational milestone | Essential experimental partner | Key scientific failure mode |
|---|---|---|---|
| Endometriosis | Donor-level, phase-stratified replication of genetically supported targets | Endometrial tissue/organoid + gynecology lab | Confounding and sample/tissue mismatch |
| ALS | Independent TDP-43 cryptic splicing map and functional-target evidence review | iPSC motor-neuron/ASO specialist | Splice correction may not rescue full disease |
| Type 1 diabetes | Evidence atlas of immune–islet cell vulnerabilities and immune protection | Islet/transplant immunology lab | Immune evasion creates unacceptable risks |
| Idiopathic pulmonary fibrosis | Cross-cohort fibroblast-state and antifibrotic pathway replication | Human lung fibroblast/ex vivo lung lab | Stopping scarring ≠ reversing it |
| ME/CFS | Cross-study molecular reproducibility and rigorous phenotype alignment | Cohort and immunology/metabolomics collaborators | Cohort heterogeneity and no robust subtypes |

## Phase 0: infrastructure (now)

- Project charter, ethics and data catalog.
- Strict donor-level schema, audit, artificial test data and CI.
- A documented, falsifiable endometriosis protocol.
- No proprietary or restricted source code copied.

## Phase 1: verified human-data reproduction (next)

1. Obtain data-rights and clinical metadata review (GEO + Zenodo); evaluate original analysis artifacts.
2. Extract local donor-level raw count tables. Preserve batch, phase, cell type and tissue provenance; never commit human data.
3. Reproduce one published result with a fit-for-purpose model and negative controls. Produce a source-versioned QC audit.
4. Resolve case/control and cohort overlap in external candidates **before** approving a holdout.
5. Freeze a target-selection model and keep holdout outcomes unobserved until that freeze.

**Exit gate:** an auditable report explains where reference and reproduction agree/disagree. If rights or phenotype metadata prevent this, publish blockers and pivot to a feasible cohort rather than fabricate results.

## Phase 2: target evidence and external evaluation

- Curate genetics-to-gene mappings from 2026 study with strong mapping support, not nearest-locus shortcut.
- Create a source-backed evidence table and sensitivity-analyze heuristic weights.
- Test candidates in a cohort with matching anatomy, hormone context, menstrual phase and adequate participant n.
- Assess effects on non-endometriosis reproductive tissue, off-target/broad-function risks and competing explanations.
- Publish nulls and conditionally shortlisted hypotheses with uncertainty, not a diagnostic product.

## Phase 3: translation by qualified partners

- Domain scientist, biostatistician, patient advocate and ethics/data-governance reviewers.
- Independently designed and documented experiments, including appropriate controls.
- Screen for assay reproducibility and safety before further laboratory development.
- Reassess whether other disease tracks present stronger reproducible opportunities.

**Operating rule:** fund and launch one next bottleneck only after the preceding evidence gate is satisfied.
