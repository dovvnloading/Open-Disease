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
- A separately tested **experimental edgeR v4 donor-count quasi-likelihood model** with estimability and explicit covariate gates ([methods and limitations](docs/EDGE_R_MODEL.md)); no real human atlas data have been analyzed.
- A separately tested R/Matrix sparse aggregator and **local-only** Seurat raw-RNA-count adapter ([instructions](docs/ATLAS_EXPORT.md)); real RDS ingestion has not been verified.
- Auditable, manually sourced *candidate* ranking, with a strict separate-study validation API.
- Explicitly synthetic data, seeded tests, and GitHub Actions.
- No background collection, automatic data downloads, diagnosis, clinical advice, or human data in the repository. An explicitly invoked, local-only GEO downloader is available.

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

## Immediate milestone: independent, public-data research (no lab prerequisite)

**Our next scientific goal is a real, reproducible analysis of existing publicly released biomedical data—not another synthetic demonstration or a request to recruit a laboratory.** See [Independent research plan](docs/INDEPENDENT_RESEARCH_PLAN.md) and [Issue #11](https://github.com/dovvnloading/Open-Disease/issues/11).

**New opt-in tooling (synthetic-tested, no real analysis yet):** `geo-download` retrieves the pinned GSE51981 Series Matrix with SHA-256 provenance, `geo-prepare` extracts a strictly unreviewed manifest, and a separate R/limma script fits reviewed, phase-aware microarray contrasts. All data and result artifacts stay local. See [GSE51981 microarray pipeline, gates and commands](docs/GEO_MICROARRAY.md).

- Start with **GSE51981**, a manageable 148-sample endometrial *microarray* dataset from GEO. Original paper: 77 endometriosis cases, 34 pathology-free controls, **37 other pelvic-pathology samples**. The last 37 must not be silently treated as healthy controls. Respect clinical and cycle-phase metadata.
- Obtain the released processed **log2 microarray intensity** Series Matrix under applicable source terms, log source SHA256 and sample QC, and reproduce a specific published result using a microarray-appropriate analysis (e.g. limma). **Do not feed array intensities into edgeR or raw RNA count tools.**
- Document exact cohort exclusions, disease-vs-control contrast, menstrual phase, effect directions, FDR, negative controls, uncertainty and mismatches with the original publication. No invented gene discoveries.
- Cross-reference **public Open Targets target/disease evidence**, but report prior findings, mapping ambiguity and association-versus-causation limits.
- Include a truly independent source where scientifically comparable, or explicitly report a validation gap. Honest null results are valuable.
- Reserve the huge GSE266265 Seurat import for a later, sufficiently resourced phase. Independent computational work does **not** require a clinician or laboratory coauthor to start. Laboratory/clinical collaboration remains necessary to claim experimentally validated therapeutic activity or a cure.

**Science status: no actual human-data analysis completed yet.** The existing programs remain synthetic-only and do not support clinical use.

## Research materials

- [PROTOCOL.md](PROTOCOL.md): hypothesis, statistical guardrails, success/failure criteria.
- [docs/EDGE_R_MODEL.md](docs/EDGE_R_MODEL.md): count-model design, assumptions, test gates, and limitations.
- [DATA_SOURCES.md](DATA_SOURCES.md): source identifiers, access and validation caveats.
- [data/sources.json](data/sources.json): machine-readable study inventory.
- [ETHICS_AND_LIMITATIONS.md](ETHICS_AND_LIMITATIONS.md): research ethics, governance and limitations.
- [docs/RESEARCH_ROADMAP.md](docs/RESEARCH_ROADMAP.md): milestones for the five-disease effort.
- [CONTRIBUTING.md](CONTRIBUTING.md): how scientists and software contributors can help.
- [LICENSE](LICENSE): MIT **for our original code only**; third-party data and code have separate terms.

## Boundaries

**NOT FOR CLINICAL USE.** No output from this repository should be used to diagnose, treat, triage, or predict outcomes for individuals. Do not post personal health records, participants' raw data or protected participant metadata to GitHub. No source code is copied from unlicensed upstream research repositories. A gene expression difference is not evidence of causal disease mechanism or drug efficacy.

See [ETHICS_AND_LIMITATIONS.md](ETHICS_AND_LIMITATIONS.md).

### Audited sample-level map (alternative offline R interface)

Alongside the existing cell-mapped Seurat exporter, scripts/r/export_seurat_pseudobulk.R accepts a separately reviewed sample-to-donor map and emits a private exclusion audit. Both exporters must agree on the same artificial raw-count data in R CI. No actual human atlas has been ingested. See [docs/ATLAS_EXPORT.md](docs/ATLAS_EXPORT.md).
