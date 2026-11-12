# Candidate external endometriosis cohorts: why a holdout is not yet ready

Research check dated 12 November 2026. Public accession metadata were reviewed; data use, donor overlap and biological comparability have **not** been verified. None of these cohorts has yet been used in this repository.

| Dataset | Potential use | Major limit |
| --- | --- | --- |
| GSE266265 | Primary untreated eutopic endometrium discovery (60 donors, 35 case/25 control) | Donor/phase strata are smaller; same cohort cannot provide genuine independent replication |
| GSE213216 | Possible independent single-cell eutopic tissue direction check | 10 eutopic samples are reported among many mixed anatomical sites; small control count; hormone, donor and phase alignment need auditing |
| GSE179640 | Possible cross-study single-cell sensitivity analysis | Study includes eutopic endometrium, ectopic lesions and organoids across ~14 participants; control counts are very small and tissue matching is difficult |
| GSE51981 | External sample-level direction check on microarray | 148 array samples, not RNA-seq counts; distinct technology and phenotype definitions, so current pseudobulk pipeline cannot process it |

Useful source URLs:

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213216
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179640
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981
- https://www.nature.com/articles/s41556-022-00961-5

## External validation rule

A source accession is **not** automatically an independent replication set. Verify individual-level independence, tissue compartment, disease definition, untreated exposure, menstrual-cycle phase, sequencing assay, effect direction, batch and cell-type ontology before the holdout is ever accessed for target selection.

The next scientific gate is not a new AI model. It is finding a sufficiently powered, rights-cleared, tissue-comparable cohort for a prospective frozen-target test. If no such cohort exists, report that limitation explicitly and seek a laboratory collaboration to generate the necessary samples rather than claiming success from mismatched public data.
