# The no-lab research plan: a useful result before a larger infrastructure project

**Adopted:** November 12, 2026  
**Active reproduction issue:** https://github.com/dovvnloading/Open-Disease/issues/18 (source audits #15 and #16 must finish first)

## Reset the project constraints

Our available resources are source code, AI-assisted literature review, open public records and ordinary computer hardware. We do **not** assume a hospital, institutional IRB, graduate students, sequencing facility, wet lab, drug funding, or formal expert partnership. A cure for complex endometriosis cannot be credibly promised by this setup. What we **can** contribute is independent computational evidence, replication/falsification work, accessible tools and clearly defined candidate hypotheses for others to validate.

The nearest-term failure mode to avoid is building ever more synthetic-data infrastructure while doing zero actual research on public disease data. A credible, limited real-data result—even a negative one—is more meaningful than a larger codebase without findings.

## Experiment 1: reproduce a real endometriosis molecular result on a normal laptop

### Public starting point

- [NCBI GEO GSE51981](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981): 148 Affymetrix GPL570 **endometrial microarrays**, processed expression values and public sample metadata. A GEO Series Matrix is available without needing 39 GB of single-cell files.
- [Original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/): 77 endometriosis samples, 34 pathology-free controls, and 37 samples with *other* pelvic conditions. These groups are not interchangeable, and the analysis includes menstrual phase.
- [GEO's download instructions](https://www.ncbi.nlm.nih.gov/geo/info/download.html) explain public Series Matrix downloads.

**Crucial:** GSE51981 is **normalized log2 microarray expression**, not RNA-seq raw counts. Do not run it through our edgeR, Seurat or pseudobulk modules. Build an independent, correctly specified microarray/limma pathway.

### The smallest credible scientific question

Can we reproduce one reported endometrial expression signature in women with endometriosis compared with pathology-free controls *within or accounting for menstrual phase*, and test whether it is also present in participants with **other pelvic pathologies**? If the signal is nonspecific, that is useful information.

### Pre-analysis checklist

1. Audit sample metadata and exact pathology coding before viewing gene results. Record case/control/other status, menstrual phase, source specimen, published criteria, unknowns, exclusions and unique sample identifier.
2. Version the processed GEO data URL, download date and SHA256 locally. No participant/sample-level export or clinical notes go into the GitHub repository or Actions artifacts.
3. Select one published biological theme to reproduce (e.g. inflammation/steroid hormone signaling), with the original paper as prior art—not a new hypothesis.
4. Fit a suitable limma model to the processed microarray intensities with independently reviewed, predeclared contrasts and menstrual-phase covariate/stratification treatment if estimable. Do not impute or invent metadata.
5. Report per-group sample counts, expression effects, uncertainty, FDR family, all filtered probes/genes, at least one negative control and the 37 other-pelvic-pathology cases as a separate specificity check.
6. Probe-to-gene identifiers can be one-to-many or uninformative. Use a documented annotation release; never assign ambiguity as unique causal gene.
7. Reproducibility report should name published signals recovered, published signals not recovered, possible confounding and sample-selection sensitivity. Release scripts and aggregate tables where permitted.

### Second-stage computational contribution

- Use [Open Targets](https://platform-docs.opentargets.org/) disease/target evidence (public API and releases), carefully distinguish source-backed disease association from variant-to-gene causality.
- Investigate both **agreement and contradiction** between gene-expression effects, prior GWAS interpretation, menstrual phase, tissue site and known pharmacology.
- Freeze a modest shortlist *before* looking at independent study outcomes. A cross-platform directional check or tissue mismatch is weaker than a true participant-independent replication and must be labeled accordingly.

### What counts as success

- A real public-data experiment that reproduces known signal(s), with exact exclusions/uncertainty and a repeatable script; **or**
- A demonstrable nonreplication or critical data-label/confounding mistake that changes a published interpretation; **or**
- A comprehensive, source-linked negative result narrowing which genes are **not yet** credible candidate targets.

No drug, biomarker, therapeutic causality or cure is established by such an experiment. Even if the output is negative, it is an honest piece of science and can save downstream effort.

## Things not on our critical path

- Downloading GSE266265's approximately 39 GB processed object.
- Finding a surgeon, university affiliation, clinical collaborator, sponsor, wet-lab partner, or restricted-data approval before we can run **public** computational experiments.
- Training a novel deep learning model or building a diagnostic web app.
- Predicting clinical benefit from AI-generated drug candidates.

External expertise becomes materially necessary **later** when claiming functional biological causal effects, experimenting with human tissue or proposing a therapeutic intervention. Expert feedback on statistical interpretation is always valuable, but the lack of such a person does not prevent a careful exploratory analysis and a transparent, qualified release.

## Resource and honesty constraints

A normal CPU notebook/workstation and public software are plausible for the GSE51981 processed matrix. The **full raw CEL archive** is listed around 711.8 MB and is optional for the first reproduction; GEO provides already processed values. Actual runtime, memory and retrieval success must be measured rather than promised.

Use a fixed limited research effort: prioritize data provenance and one reproducible report before expanding functionality. Record failed download, metadata or sample-count barriers openly. Do not claim we analyzed real data until an actual controlled run finishes, its outputs are checked, and the source records are audited.

**The result we can reasonably hope for is a defensible independent computation or falsification, not an independently invented cure.**
