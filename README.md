# Open-Disease

**Open-source research infrastructure for discovering reproducible disease mechanisms.**

**Status: preclinical research software, v0.1.0.** This repository has not discovered a cure, analyzed real endometriosis donors, or validated a diagnostic test. The runnable example uses explicitly simulated data.

## Mission

Build transparent, testable computational methods to help scientific teams investigate serious diseases without broadly available cures. Prioritize independent replication, patient benefit, and biological safety over impressive-looking model scores.

Our five initial research priorities:

1. **Endometriosis — ACTIVE FIRST PROGRAM:** genetics-informed, donor-level cell-state and target discovery.
2. **ALS:** TDP-43 and cryptic RNA-splicing defects (exploratory program).
3. **Type 1 diabetes:** beta-cell replacement and immune protection (exploratory program).
4. **Idiopathic pulmonary fibrosis:** pathological fibroblast states and scarring (exploratory program).
5. **ME/CFS:** reproducibility of genetic, immune and metabolic signatures (exploratory program).

These are research directions, not five parallel clinical cure claims. See docs/CHARTER.md and docs/ROADMAP.md.

## What works today

- Validation of donor identity and metadata across gene rows.
- Cell-to-donor pseudobulk aggregation with explicit zero-count handling for small datasets.
- Exploratory per-study, per-cell-type, per-cycle-phase donor-level differential expression (Welch test on log2 CPM with Benjamini–Hochberg correction).
- Hypothesis ranking from **externally curated gene-level mapping evidence**; a GWAS locus alone is not accepted as gene-level evidence.
- Frozen candidate list evaluated on a **different study**, reporting both success and failure.
- Fully synthetic two-cohort demo, deterministic tests and GitHub Actions CI.

**Not implemented:** automatic download and parsing of the ~39 GB tissue atlas, automatic GWAS fine-mapping, covariate-adjusted count-model analysis, independent human cohort validation, clinical tests or wet-lab assays. Results from the baseline statistical model are exploratory only.

## Quick start

Python 3.11 or later:

    git clone https://github.com/dovvnloading/Open-Disease.git
    cd Open-Disease
    python -m venv .venv
    source .venv/bin/activate
    python -m pip install -e ".[dev]"
    python -m pytest -q
    bash scripts/run_demo.sh

Windows: activate the virtual environment using the corresponding PowerShell command and run each demo command from scripts/run_demo.sh manually.

Read the example outputs under results/synthetic-demo/. Every simulated gene is named SYNxxx and is biologically fictitious. The two simulated studies share planted effects; this exercises software controls and **does not constitute external scientific replication**. Generated results are gitignored.

### Commands

    python -m open_disease simulate --out-dir results/synthetic-demo
    python -m open_disease aggregate --cells data/raw/cell_counts.csv --out data/processed/pseudobulk.csv
    python -m open_disease analyze --pseudobulk data/processed/pseudobulk.csv --out results/de.csv
    python -m open_disease rank --de results/de.csv --evidence data/processed/curated_genetics.csv --discovery-study STUDY_A --phase proliferative --cell-type epithelial --out results/ranked.csv
    python -m open_disease validate --ranked results/ranked.csv --de results/de.csv --validation-study STUDY_B --out results/validation.csv
    python -m open_disease sha256 data/processed/pseudobulk.csv

The aggregate command is intended for *small, permitted datasets*, not to expand 39 GB of single-cell objects in memory. For genuine data, perform pseudobulk aggregation within the original analysis environment, with transcriptome-wide library-size normalization, and export a complete donor-by-gene count matrix in the documented schema.

The rank command requires a **curated** genetics evidence CSV. It does not infer causality from GWAS p-values and must not be run with fabricated scores for real genes. Rank weights are a deliberately simple heuristic for testing a research workflow, not a probability of target success.

## Read before analyzing human data

- docs/PROTOCOL.md — hypothesis, validation rules, stop/go decisions and remaining statistical requirements.
- docs/DATA_SOURCES.md — real source URLs, access controls, citation and licensing review.
- docs/SCHEMA.md — exact data formats and scientifically important assumptions.
- docs/ETHICS_AND_LIMITATIONS.md — patient privacy, research integrity and diagnostic boundaries.
- docs/ROADMAP.md — milestones and collaboration needs.
- CONTRIBUTING.md — review and contribution standards.

**Never commit patient data, genetic sequences, patient-level pseudobulk, consent forms, credentials, raw images or protected annotations.** A pseudonymous donor ID is not an assurance of anonymity. Protected research may require ethics approval, access agreements, controlled computing and institutional collaborators.

## License and scientific use

The repository's existing MIT license remains intact and governs original software, not third-party datasets or patents. Every external dataset and tool requires its own license and consent review. This project does **not** provide medical advice or establish treatment efficacy.

Research questions and corrections are welcome through GitHub issues. Do not post personal medical histories or genetic information.
