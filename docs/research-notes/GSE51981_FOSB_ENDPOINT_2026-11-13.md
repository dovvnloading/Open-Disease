# GSE51981 prespecified FOSB comparison

**November 13, 2026.** This amendment to the [analysis plan](GSE51981_PREANALYSIS_PLAN_2026-11-12.md) fixes the published finding, probe and decision rule for [issue #18](https://github.com/dovvnloading/Open-Disease/issues/18) before inspection of expression outcomes. It adds a single primary endpoint. Eligibility still depends on the sample audit.

## Published reference and probe

Tamaresis et al. list FOSB first in the upregulated section of Table 4. For endometriosis versus pathology-free controls they report fold changes of 34.1, 10.2 and 17.0 in proliferative, early-secretory and mid-secretory endometrium, respectively. The other-pathology group also had higher values than controls in each phase. These are the authors' fold changes. The planned analysis estimates a common adjusted log2 difference across phases, so it tests directional agreement and does not seek to match those three magnitudes. [Original study, Table 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/#T4)

The NCBI-generated GPL570 annotation dated August 9, 2016 maps the symbol FOSB and Entrez Gene 2354 to exactly one probe set, `202768_at`. The archive SHA-256 is `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394`. Table 4 identifies the gene by symbol. The probe assignment here comes from the pinned platform annotation. [Official GPL570 annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz)

## Decision fixed before the fit

- **Population and model.** Use the original processed log2 GC-RMA intensities and reviewed, independent, eligible donors. Retain only proliferative, early-secretory and mid-secretory phases. Apply the existing additive limma model, complete-case probe rule and seed 20261112. Keep other pelvic pathology separate from pathology-free controls.
- **Primary result.** Read the `202768_at` row of `case_vs_control.csv`. Report its log2 difference, moderated 95% confidence interval, unadjusted p-value, within-contrast Benjamini–Hochberg q-value and q-value across performed contrasts. Within-contrast adjustment includes every tested complete-case probe. There is one designated probe and one primary contrast.
- **Interpretation.** A positive estimate with within-contrast q ≤ 0.05 supports directional agreement under this model. A negative estimate meeting that threshold contradicts the expected direction under this model. Any other estimable result is inconclusive. Missing or excluded probe, unresolved cohort eligibility or an inestimable contrast is a blocked result. Never substitute another probe after seeing the outcome.
- **Specificity and controls.** Report the same probe in the other-pathology versus control contrast if estimable, and the fixed-seed within-phase permutation diagnostic. A shared direction across the clinical groups would limit disease specificity. The secondary result cannot rescue the primary endpoint.

This is a same-cohort reanalysis. A successful fit would still leave independent replication, causality and clinical performance untested.
