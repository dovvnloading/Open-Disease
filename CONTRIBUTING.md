# Contributing to Open-Disease

We welcome expertise in endometriosis biology, human genetics, statistical genomics, computational biology, reproducibility engineering, reproductive health and patient-centered research.

## Good first contributions

- Independently verify source access, licenses, donor counts and menstrual-cycle labels without committing human data.
- Suggest **published**, ethically reusable comparator cohorts and explain anatomical/clinical comparability.
- Contribute tests for scientific failure modes (pseudoreplication, leakage, duplicate donors, partial gene matrices).
- Reproduce a reference finding under a documented protocol with counts-aware methods.
- Review target mapping logic, statistical methods, reproductive safety and accessible research communication.

## Pull-request rules

1. Open an issue first for scientific hypothesis changes and describe the falsifiable question and data rights.
2. Run `python -m pytest -q`; include regression tests for code changes.
3. List source URLs, accession IDs, code/data licenses and versions. Distinguish observed data from synthetic or simulated examples.
4. Report negative, null and nonreplicated outcomes; do not selectively omit them.
5. Do not include records containing participant identifiers or individual medical history, any restricted dataset, secrets, or upstream code without compatible license.
6. No medical advice or clinical promises in README, releases or submissions.

Scientific proposals must state one primary endpoint, inference unit, planned covariate handling, power limitations, negative controls and what would falsify the hypothesis.
