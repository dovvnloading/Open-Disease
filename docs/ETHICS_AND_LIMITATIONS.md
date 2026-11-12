# Ethics, patient protection and scientific limitations

Open-Disease is nonclinical, exploratory software. It does not diagnose endometriosis, predict personal disease status, identify a cure, recommend a treatment, or replace a doctor or research ethics committee.

- **No private health data in GitHub.** Do not upload raw reads, person-level pseudobulk counts, raw sequence, donor-level metadata, patient histories, medical images, login credentials, or consent forms. Real data, even pseudonymized, may be sensitive. Gitignore is only an accidental-commit guard, not a permission system.
- **Human subjects:** Any donor enrollment, collection, access to controlled genomic data or prospective diagnostic research must be supervised by qualified institutions with relevant ethics/consent approvals and access agreements.
- **Disability and health equity:** Include patients and advocates in choosing outcomes. Validate across populations, disease subtypes and care settings; do not equate a narrow laboratory signal with general clinical utility.
- **No false causal claims:** Genetic association and differential expression are clues, not demonstrations that gene perturbation changes lesion biology. A target may be necessary for normal fertility or immunity.
- **No clinical efficacy claims:** In vitro, in silico or animal findings are not proof of safe treatment in people. Real interventional work must be designed and supervised by relevant biological and clinical experts.
- **Baseline statistical limits:** No covariate-adjusted count model, external real-world replication, multi-family correction, sample-size design or patient benefit measure is yet implemented. Every output must retain these limitations.
- **Publication:** Predefine outcomes and record failures. Report small sample sizes, missing data, confounding, ancestry, menstruation phase, batch and code/data versions. Do not optimize candidate lists using a nominally untouched holdout.
- **Third-party rights:** The software's MIT license does not grant rights to restricted data, other analysis code, patents, or patient-derived samples.

If you discover an exposure of patient data, stop dissemination and notify the dataset custodian and repository maintainers through a private channel. Do not describe sensitive contents in a public issue.
