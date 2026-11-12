# Data source and access registry

Last metadata review: 12 November 2026. **No real participant data have been ingested into this repository.** License/terms review must happen before use or redistribution. Public metadata does not imply permission to rehost individual-level files.

| Role | Study / ID | Published evidence and access | Readiness |
| --- | --- | --- | --- |
| Primary tissue discovery | GEO GSE266265, PRJNA1106473 | 2026 60-donor atlas (35 cases, 25 controls); public raw and processed Seurat RDS; primary processed atlas ~38.9 GB and epithelial subset ~190.7 MB | Metadata verified; file format, software compatibility, licenses and donor metadata not audited |
| Genetic mechanism anchors | 2026 multi-ancestry GWAS (Nature Genetics), Zenodo 18983492 | ~1.4 million women, 105,869 cases and 80 genomic regions; summary statistics linked on Zenodo | Public summary-statistics link verified in paper; contents and licensing review pending |
| External exploratory diagnostic translation | 2026 menstrual-fluid single-nucleus pilot, Zenodo 17203106 | 5 affected + 5 control donors; processed data/code shared; donor-level sequencing files require committee approval + data transfer agreement at DOI 10.57745/ELWEGD | Pilot is underpowered and not valid for diagnostic claims; access review required |

Source links:

- NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265
- Published atlas: https://www.nature.com/articles/s41467-026-73020-4
- Atlas original code: https://github.com/DuempelmannLea/endometriosis_endometrium_scRNA_atlas
- GWAS paper: https://www.nature.com/articles/s41588-026-02582-2
- GWAS summary data record: https://doi.org/10.5281/zenodo.18983492
- Menstrual-fluid pilot paper: https://link.springer.com/article/10.1186/s12916-026-05020-6
- Processed pilot data/code: https://doi.org/10.5281/zenodo.17203106
- Controlled raw pilot data: https://doi.org/10.57745/ELWEGD

## Data intake checklist — mandatory for every real dataset

1. Record research question, immutable data version, access date, DOI/accession, publication, principal investigators, genome build, gene annotation release and provenance.
2. Confirm license and terms permit the contemplated local analysis; verify human subjects consent and whether IRB/data-use agreements apply. Never clone another lab's analysis code without reviewing its license.
3. Obtain donor-level **non-identifying** clinical covariates and verify study-level and donor-level independence; maintain access controls.
4. Hash each downloaded file with the local sha256 command and store checksums in a controlled research log. A checksum for unaccessed content must not be guessed or generated.
5. Document RNA count layers, cell filtering, hormone-exposure definitions, cycle phase, tissue source, covariate availability and cell-type mapping.
6. Export donor-level pseudobulk in the schema (SCHEMA.md). Check raw sequencing data handling under source agreements; never upload human gene-level counts or patient metadata to issues or the repository.
7. Obtain independent study/cohort for a genuine external evaluation. A pilot with ten donors cannot establish clinically useful accuracy.

### Important conflict and patent note

The atlas paper discloses co-inventor interests in patent applications relating to clinical uses of gene signatures. Open-source research is **not** a freedom-to-operate determination. Obtain qualified legal review before commercialization or deployment.

### Access and scale

A ~39 GB RDS is not an appropriate automatic CI download. The repository deliberately includes no real-data fetch automation until access, licensing, storage, metadata and checksum checks are complete.
