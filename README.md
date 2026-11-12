# Open-Disease

**An open-source, evidence-first research effort to accelerate understanding of diseases without reliable, generalizable cures.**

**Status: research infrastructure under active development. First active program: endometriosis. No treatment, diagnostic test, or validated therapeutic target has been produced.**

Our long-term disease shortlist: (1) endometriosis, (2) ALS, (3) type 1 diabetes, (4) idiopathic pulmonary fibrosis, and (5) ME/CFS. We do not claim that software alone can cure any of them. Our immediate deliverable is reproducible biomedical evidence and independent validation of hypotheses that experimental scientists could test.

## First question: endometriosis

*Can genetic evidence and donor-level cell-type expression analysis identify candidate endometriosis mechanisms that persist across comparable, independent studies?*

The primary human single-cell atlas contains 466,371 cells from **60 donors** (35 cases, 25 controls). Its individual cells are not 466,371 independent patients. The 2026 multi-ancestry GWAS identifies 80 risk regions, but a risk region alone does **not** implicate its nearest gene as a therapeutic target.

Sources: [GEO GSE266265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265), [2026 atlas article](https://doi.org/10.1038/s41467-026-73020-4), [2026 GWAS](https://doi.org/10.1038/s41588-026-02582-2), and [GWAS summary statistics](https://doi.org/10.5281/zenodo.18983492).

## What works today

- A **strict donor-level** cell-count aggregation adapter for small long-format input tables.
- An exploratory, phase- and cell-type-stratified donor pseudobulk baseline (Welch test on log2(1+CPM), with Benjamini–Hochberg correction). **Not** a validated disease model or a publication-grade replacement for edgeR/DESeq2.
- Structural cohort preflight auditing: missing genes, minimum donors, cell and library counts.
- A separately tested R/Matrix sparse aggregator and **local-only** Seurat raw-RNA-count adapter ([instructions](docs/ATLAS_EXPORT.md)); real RDS ingestion has not been verified.
- Auditable, manually sourced *candidate* ranking, with a strict separate-study validation API.
- Explicitly synthetic data, seeded tests, and GitHub Actions.
- No background collection, API calls, data downloads, automatic diagnosis, clinical advice, or human data in the repository.

## Run the code locally

Requires Python 3.11 or newer.

```bash
git clone https://github.com/dovvnloading/Open-Disease.git
cd Open-Disease
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
python -m pytest -q
bash scripts/run_demo.sh
python scripts/audit_cohort.py \
  --pseudobulk results/synthetic-demo/pseudobulk.csv \
  --out results/synthetic-demo/audit.csv
```

The demo writes **SIMULATED ONLY** data into an ignored directory. Any discovery, target ranking or "replication" in synthetic data is a software test, **not** biological evidence.

Use `python -m open_disease --help` to see commands. The core input contract is documented in [DATA_SOURCES.md](DATA_SOURCES.md) and [PROTOCOL.md](PROTOCOL.md).

## The first non-synthetic milestone

1. Independently check licenses, metadata, ancestry/cohort overlap, tissue matching and permissions for the human datasets.
2. Export **raw count** donor-by-cell-type pseudobulk from the 2026 atlas locally, preserving the original source file and analysis script. The full Seurat RDS is about 38.9 GB; do not download it in GitHub Actions.
3. Preflight each phase and cell type; mark underpowered strata as ineligible.
4. Replicate an established result with a count-aware statistical model and pre-specified confounders before claiming novelty.
5. Lock a target selection plan, then evaluate a **truly independent and comparable cohort**. Publish negative findings, uncertainty, and exclusions.

Existing code handles only steps of this workflow; it **has not yet accessed or analyzed human datasets**.

## Research materials

- [PROTOCOL.md](PROTOCOL.md): hypothesis, statistical guardrails, success/failure criteria.
- [DATA_SOURCES.md](DATA_SOURCES.md): source identifiers, access and validation caveats.
- [data/sources.json](data/sources.json): machine-readable study inventory.
- [ETHICS_AND_LIMITATIONS.md](ETHICS_AND_LIMITATIONS.md): research ethics, governance and limitations.
- [docs/RESEARCH_ROADMAP.md](docs/RESEARCH_ROADMAP.md): milestones for the five-disease effort.
- [CONTRIBUTING.md](CONTRIBUTING.md): how scientists and software contributors can help.
- [LICENSE](LICENSE): MIT **for our original code only**; third-party data and code have separate terms.

## Boundaries

**NOT FOR CLINICAL USE.** No output from this repository should be used to diagnose, treat, triage, or predict outcomes for individuals. Do not post personal health records, participants' raw data or protected participant metadata to GitHub. No source code is copied from unlicensed upstream research repositories. A gene expression difference is not evidence of causal disease mechanism or drug efficacy.

See [ETHICS_AND_LIMITATIONS.md](ETHICS_AND_LIMITATIONS.md).
