# GSE51981 public-metadata audit — first pass

**Date:** 2026-11-12  
**Task:** `gse51981-sample-audit` / [#15](https://github.com/dovvnloading/Open-Disease/issues/15)  
**Evidence level:** source-audit only. **No matrix downloaded, sample-specific manifest verified, donor-level analysis run, or biological result generated.**

## What the sources actually establish

- [GEO GSE51981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981) lists 148 publicly submitted endometrial **microarray samples** on GPL570, and describes samples from women across menstrual phases. The number of GEO samples is not by itself verification of one independent donor per column.
- [Original 2014 paper, Materials and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/) gives **77 endometriosis, 37 without endometriosis but with other uterine/pelvic pathology, and 34 without endometriosis or uterine/pelvic pathology**. The 37 cannot be reassigned to healthy controls.
- That paper states an exclusion for **hormonal treatment within the previous three months**, menstrual-phase assignment based on histology/hormones, and **four indeterminate phase samples** that were only used in phase-unrestricted classifiers. These are *cohort-level publication statements*, not proof that each GSM record meets each analysis field.
- Individual [example control](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1256727), [other-pathology](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1256765) and [endometriosis](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1256685) GEO records expose separate `tissue`, `endometriosis/no endometriosis`, `presence or absence of uterine/pelvic pathology`, and/or `endometriosis severity` characteristics. These *examples* prove such metadata exist; **they do not audit all 148 records**. Original GSM `VALUE` is labeled **base-2 logarithm GC-RMA**, not raw RNA counts.
- The [2022 Bunis et al. reanalysis, Methods 4.1.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC8766492/) excluded other-pathology records and ambiguous severity/lab-source/phase assignments, ultimately analyzing **71 endometriosis cases plus 34 no-pathology controls** (105 total). This is a documented **prior reanalysis of the same GEO cohort**, not independent replication. Its separate processing uses justRMA/batch correction, not necessarily the GEO-provided GC-RMA values.

## Preliminary field-by-field decision

| Intended field | Candidate primary evidence | Decision before analyzing effects |
| --- | --- | --- |
| GEO sample accession | Unique `GSM...` in Series Matrix header, cross-check against `!Sample_geo_accession` | Parser already enforces match; **not yet run on actual matrix** |
| Disease / other pathology | Individually sourced `!Sample_characteristics_ch1` and GEO GSM pages; original group definitions | Review all 148; case/control/other **must not** be inferred from a sample's numeric title |
| Menstrual phase | GSM `tissue` characteristic and original study's histologic phase protocol | Confirm allowable PE/ESE/MSE individually; indeterminate or inconsistent phase excluded from phase-aware analysis |
| Tissue | GSM characteristics, original paper | Require verified eutopic endometrium; do not mix lesion samples |
| Hormonal exposure | Original paper's cohort exclusion | Do not assume every GSM is verified `none` without documented dataset-to-publication linkage and an explicit evidence rule |
| Person-level independence | Original paper's description of samples from women, plus any published linked metadata | Aggregate study description supports but does not prove **per-GSM** uniqueness. Require reproducible donor reconciliation or a documented limitation, not GSM=donor by fiat |
| Processing / expression scale | Individual GEO sample `VALUE base-2 logarithm GC-RMA` and Series Matrix | Verify matrix-provided scale before `--log2-confirmed yes`; never use edgeR on array intensities |
| Reuse/redistribution | GEO records, original article, relevant source terms | Public access verified. **File-specific redistribution terms and final sample-level handling not yet audited**; do not commit matrix/manifest |

## Stop/go decision

**NO-GO for a reported real-data limma result today.** A per-GSM clinical/donor map, rights record, matrix integrity/checksum, fixed exclusions and frozen reference comparison are not yet completed. The current `geo-prepare` tool correctly leaves all samples unreviewed. The experiment can proceed *locally* once a fully cited audit supports inclusion, phase, person-level independence and dataset reuse; a missing field means exclude or disclose/stop, not invent.

## Next bounded agent work

1. In an explicitly authorized public-data environment, obtain the official GSE51981 Series Matrix via the existing opt-in downloader and record compressed SHA-256/version. Processed matrix and source characteristics remain in ignored private data directories.
2. Audit all 148 public GSM entries against individual GEO source characteristics and original paper, documenting discrepancies **locally**, then publish only aggregate field-coverage/exclusion counts and source URLs.
3. Resolve person-level sample mapping and source processing/rights; if unresolved, document blocker and refrain from fitting the model. Avoid using the 2022 reanalysis as independent validation.
4. Freeze first comparison/phase rules and the published outcome to reproduce **before viewing resulting gene effects**. Record an explicit go/no-go decision with reasons.

No outreach, grant, patient contact or private human-data acquisition is called for. Source examples and published totals are not a completed cohort audit.
