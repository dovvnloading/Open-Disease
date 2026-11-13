# Open-Disease

**An independent, AI-first computational disease research workbench.**

Open-Disease is a place to aim research agents and bounded compute at difficult disease questions. Agents can audit public papers/data, write and test analyses, reproduce known findings, explore falsifiable hypotheses, and publish the *methods, failures, uncertainties and source-backed evidence*. The goal is useful work that future agents or researchers can build on—not a grant-funded lab, research consortium, clinical service or claim to have cured anything.

**Current status (November 2026): software and synthetic tests exist; no genuine human-data analysis, confirmed new molecular target, clinical test or cure has been produced by this repository.** Endometriosis is the only active disease research program. ALS, type 1 diabetes, pulmonary fibrosis and ME/CFS are possible later tracks, not staffed or validated projects.

## Give an agent a concrete task

Read [AGENTS.md](AGENTS.md) for operating rules and [the agent workflow](docs/AGENT_WORKFLOW.md) for evidence/handoff standards. Select a bounded task from [research/tasks](research/tasks) instead of asking an agent to "solve a disease":

    python scripts/research_queue.py validate
    python scripts/research_queue.py ready

Each JSON task defines a falsifiable question, sources, dependencies, acceptance criteria and a rough resource budget. An agent can create a **local, ignored** run record with `python scripts/research_queue.py new-run --task TASK_ID --agent AGENT_ID --out results/TASK_ID/run.json`. The queue is a lightweight coordination format, **not** an automatic scheduler or a compute sandbox. A task does not automatically authorize network access, external API costs, high-memory workloads or contact with people.

Starting tasks:

- Audit public **GSE51981** sample/phenotype/menstrual-phase labels and person-level independence before analysis.
- Audit **GPL570** probe annotations without forcing ambiguous probes to causal genes.
- Run existing synthetic test suites and document failures.
- After source audits pass, attempt a genuinely sourced **GSE51981 microarray reproduction** (not a novel drug discovery). A negative or blocked result is worth reporting.

See [Issue #11](https://github.com/dovvnloading/Open-Disease/issues/11) and [the independent research plan](docs/INDEPENDENT_RESEARCH_PLAN.md). We use existing public computational data **without requiring grants, recruiting people, contacting experts or running a lab**. Validation gaps remain explicit; formal scientific/clinical claims require correspondingly stronger evidence and oversight.

## What works today

- A seeded synthetic end-to-end demo, structural cohort QC, source catalog and CI tests.
- Strict donor-level RNA pseudobulk utilities and exploratory differential-expression baseline, with an experimental edgeR v4 path ([model limitations](docs/EDGE_R_MODEL.md)). These have **not** analyzed the real large Seurat atlas.
- Sparse R/Seurat export utilities tested on artificial counts ([high-memory export guide](docs/ATLAS_EXPORT.md)).
- An opt-in official GEO **GSE51981** processed-array downloader, provenance/manifest prep, and phase-aware R/limma implementation ([local-only workflow](docs/GEO_MICROARRAY.md)); **synthetic-tested only**, no audited real-GEO outcome yet.
- Heuristic evidence ranking and strict study-separation checks as software prototypes, **not** established target biology.

### Run existing software

Requires Python 3.11+; for R workflows, install the relevant Bioconductor dependencies separately.

    git clone https://github.com/dovvnloading/Open-Disease.git
    cd Open-Disease
    python -m venv .venv
    source .venv/bin/activate
    python -m pip install -e '.[dev]'
    python -m pytest -q
    bash scripts/run_demo.sh
    python scripts/verify_catalog.py

The demo generates **SIMULATED ONLY** data in ignored `results/`. No automated real-data collection or biological claim results from installation or tests.

## Research directions and science boundary

The immediate question is whether source-verified, phase-aware endometrial microarray results can reproduce a specific published signal, with other pelvic pathology **distinct from** healthy controls. Study: [GEO GSE51981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981); [published paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/). Never treat arrays as RNA-seq counts. The large 60-donor [GSE266265](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265) atlas and genetic mapping are later optional paths, not blockers to the first experiment.

Computational evidence can suggest or rule against hypotheses. It does not by itself prove causal mechanisms, treatment benefit, clinical safety or a cure. Human/sample-level data and restricted files are never committed. A public source must still be checked for permission and data quality. See [protocol](PROTOCOL.md), [limitations](ETHICS_AND_LIMITATIONS.md), [data sources](DATA_SOURCES.md), [charter](docs/CHARTER.md), and [roadmap](docs/RESEARCH_ROADMAP.md).

## Contribute

Agents and people can take a small task, submit a PR with reproducibility notes and tests, or record a well-supported blocker. See [CONTRIBUTING.md](CONTRIBUTING.md). No outreach, lab partnership or formal affiliation is expected for exploratory computational contributions. **Not for clinical use.** Original code is MIT licensed; upstream datasets and code retain their own rights.
