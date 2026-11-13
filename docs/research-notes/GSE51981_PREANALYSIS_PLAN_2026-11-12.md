# GSE51981 first computational reproduction: analysis-plan freeze

**Draft 2026-11-12, before any GSE51981 expression analysis by Open-Disease.** Issue #18 remains gated by #15 (sample audit) and #16 (annotation audit). This dated plan is *not* a registered protocol, evidence of a completed human study, or external replication. Any modification after inspecting an outcome must be versioned and marked exploratory.

## Hypothesis and comparisons

Question: can an association in the original processed GSE51981 log2 GCRMA **probe-level** data be reproduced for endometriosis against pathology-free endometrial controls, with menstrual-phase adjustment? Secondary question: do comparable expression differences also occur in the *other pelvic pathology* group?

Primary contrast: **case minus pathology-free control**.
Secondary contrast (if approved data and estimable design permit): **other pelvic pathology minus pathology-free control**.
Never call the 37 other-pathology records healthy controls. The 2014 paper's immune-activation and steroid/thyroid-signaling interpretations are **prior findings, not new Open-Disease discoveries**. No confirmatory gene panel or pathway-testing family has yet been defensibly frozen, so any later gene/theme comparison without a prespecified mapping and testing plan must be labeled exploratory.

Sources: https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981, https://pmc.ncbi.nlm.nih.gov/articles/PMC8766492/.

## Locked cohort/model choices before observing expression effects

1. **Biological unit:** a source-reviewed, non-overlapping participant, not a GSM accession alone. One reviewed included array per unique donor key for the existing analysis design.
2. **Eligibility:** individually source-supported eutopic endometrium, disease/pathology category and no recent hormonal treatment; verified proliferative, early-secretory or mid-secretory phase. Exclude the four original indeterminate-phase samples and any other unresolved or incompatible entries. Record all exclusions before reading outcomes.
3. **Input scale:** confirm the official submitted *processed log2 GCRMA* matrix and exact source SHA-256. Do not substitute RNA count edgeR, raw CEL or the differently normalized 2022 reanalysis.
4. **Main model:** limma moderated case-control difference with phase fixed effects, subject to per-phase group overlap and estimability. This additive model is **not** an exact reproduction of the original paper's separate within-phase analyses.
5. **Secondary pathology specificity:** retain other-pathology as an independent class, not a pooled control. Include it only when data support the declared contrast.
6. **Multiplicity and reporting:** return all assayed probes, exclusions, signed log2 effect sizes, uncertainty, BH q-values within contrast and across performed contrast families; preserve nulls and failures.
7. **Negative control:** one fixed-seed within-phase label permutation as a stress check, not calibrated permutation significance or external validation.
8. **Gene mapping:** pin GPL570 release and decide a gene/probe aggregation policy *before* reading signed effects. Do not select the lowest-p-value probe after results, or treat a nearest gene or a multi-gene probe as causal.

## Unfinished before the first run

The source-verified sample/donor manifest, rights review, file bytes/checksums, eligible cohort count table, annotation release and **specific original published outcome or gene/pathway test family** remain unresolved. These must be independently established and frozen by a dated commit before a result can be presented as a planned reproduction. If not available, the correct result is a documented **blocked study**, not fabricated discovery.

No actual human-data analysis or new disease mechanism has been generated under this note. The project is independent and computational; no recruitment, grants or unsolicited outreach are requested.