# Source audit — 12 November 2026

**Metadata only.** We have not downloaded or analyzed any human donor data.

### Tissue: GEO GSE266265

The peer-reviewed atlas reports 466,371 cells from 60 donors (35 cases/25 controls). GEO lists a ~38.9 GB RDS and a ~190.7 MB epithelial subset. Paper and GEO agree on GSE266265.

A significant provenance anomaly: the authors' GitHub README mentions **GSE111976** in prose while linking to **GSE266265**. The authoritative paper and GEO record use GSE266265. Upstream source: https://github.com/DuempelmannLea/endometriosis_endometrium_scRNA_atlas.

The authors' README explicitly reports **no formal code license**. We have inspected public metadata field names but have not copied their analysis code. The paper discloses patent applications involving clinical uses of gene signatures.

The original method uses carefully selected cycle-phase subsets and batch-aware modeling. Our current Welch logCPM method is a software baseline, not an independent reproduction of the paper.

### Genetics: 2026 multi-ancestry GWAS

Nature Genetics, 29 April 2026, reports 80 associated regions, ~1.4 million participants including 105,869 cases. Summary statistics are linked in the paper (DOI 10.5281/zenodo.18983492); exact file layouts, genome build and licenses still need verification.

### Menstrual-fluid pilot: 2026 study

Five cases and five controls. Processed data/code at DOI 10.5281/zenodo.17203106; individual sequencing at DOI 10.57745/ELWEGD is controlled and needs a data-transfer agreement. The pilot cannot independently establish clinical diagnostic performance.

### Next checks

1. Review file-level reuse terms and licenses.
2. On approved computing, inspect Seurat metadata/covariates, assay layers and counts.
3. Pre-specify external study eligibility; splitting cells or donors within one source is not external replication.
4. Engage a clinician-scientist and independent biostatistician before human-data interpretation.

Primary references: https://www.nature.com/articles/s41467-026-73020-4 ; https://www.nature.com/articles/s41588-026-02582-2 ; https://link.springer.com/article/10.1186/s12916-026-05020-6
