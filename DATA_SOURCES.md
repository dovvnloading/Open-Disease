# Data sources and reproducibility contract

**Last evidence audit: 2026-11-12.** Links are catalog entries, not authorizations to download, republish, or clinically use any material. Machine-readable metadata: [data/sources.json](data/sources.json).

## Executed and audited FOSB comparisons

- GSE25628 contributed one eutopic biopsy per person by the published Italian study design, with eight cases and six surgical controls. The single-probe analysis found a positive, imprecise estimate. [Results and input checksums](docs/research-notes/GSE25628_FOSB_RESULTS_2026-11-13.md)
- GSE141549 and the authors’ EndometDB archive support a donor-linked eligibility audit with recorded hormone status and histology. Platform assignment and probe provenance limit the proposed comparison. [Finnish cohort audit](docs/research-notes/GSE141549_ELIGIBILITY_2026-11-13.md)

## Primary: human endometrium atlas

**GSE266265 — Duempelmann et al., Nature Communications (2026).**

- [Authoritative GEO accession](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266265) · [peer-reviewed article](https://doi.org/10.1038/s41467-026-73020-4) · [authors' analysis repository](https://github.com/DuempelmannLea/endometriosis_endometrium_scRNA_atlas).
- 466,371 single cells, from **60 donor participants**: 35 endometriosis, 25 non-endometriosis; biopsied eutopic endometrium from participants without exogenous hormone treatment, and multiple cycle phases.
- GEO lists a **38.9 GB** processed Seurat object and a **190.7 MB** separate epithelial trajectory object. The smaller object may not provide a suitable raw donor × gene matrix; **do not substitute it without auditing assays and metadata**.
- The authors' public README sometimes mentions a different GEO accession, **GSE111976**, while its own URL and the paper identify **GSE266265**. Follow the actual GEO series and article, not the apparent README typo.
- The authors' analysis repository does **not** state a formal software license. Do **not** copy upstream source into this MIT repository without explicit permission. Independently implement documented concepts.
- Stage: metadata and rights audit pending; no human files downloaded, tested, or mirrored into GitHub.

## Genetic data (not participant genotype uploads)

**Koller et al., Nature Genetics (2026).**

- [Published paper](https://doi.org/10.1038/s41588-026-02582-2) · [GWAS summary-statistics release](https://doi.org/10.5281/zenodo.18983492).
- Approximately 1.4 million participants (105,869 cases); 80 risk regions, including 37 newly identified regions in this study. Study includes different ascertainment definitions and ancestry groups.
- The paper says its GWAS results/score weights are available on Zenodo, but some constituent cohorts (including 23andMe) have distinct restricted-access policies. **Public summary statistics are not raw individual genotypes.**
- Stage: verify exact released file manifest, assembly, allele definitions, mapping, restrictions, gene mapping/colocalization inputs, and source-version hashes before use. Do not naively map risk regions to nearest genes.

## Potential external expression cohorts — comparability is NOT assumed

**GSE213216 — Lawrenson et al. (2023).** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213216) · [paper](https://doi.org/10.1038/s41588-022-01254-1).

- 54 **specimens from multiple anatomical sites**: 10 eutopic endometrium samples (7 patients with disease, 3 controls), plus peritoneal lesions, ovarian endometriomas and unaffected ovary tissue; 6 participants in the 21-person cohort were taking exogenous hormones at surgery.
- Cannot be used as a blanket 54-donor uterine case/control validation group; tissue and cycle-phase restrictions may shrink eligible n below power thresholds. Stage: **conditional**, pending explicit sample audit.

**GSE179640 — independent endometriosis atlas.** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179640).

- Public studies have used a subset comparing 9 eutopic endometriosis and 3 control endometrium samples. This is a **reported subset**, not a verified count for our analysis; review original sample-level metadata, hormones, cycle and cohort overlap. Stage: conditional, likely underpowered.

**GSE51981 — uterine expression microarray study.** [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981).

- A 2026 secondary analysis describes **77 case / 34 healthy-control** samples after excluding uterine-pelvic-pathology controls. The cited inclusion exclusions and original study labels must be independently checked.
- Microarray expression is **not raw RNA-seq count data** and must never be run through the current count-based pseudobulk pipeline as if they were counts. Cross-platform gene-level direction could be a separate, pre-specified analysis.

## Menstrual fluid — separate biological endpoint

**Leap et al., BMC Medicine, June 2026.** [paper](https://doi.org/10.1186/s12916-026-05020-6), [processed data/code](https://doi.org/10.5281/zenodo.17203106), [individual sequencing archive](https://doi.org/10.57745/ELWEGD).

- 5 donors with disease and 5 without. This is a **small pilot**, not a clinical diagnostic validation cohort. Reported immune/epithelial and co-morbidity differences need replication.
- Individual-level sequencing files require Data Access Committee approval and a signed Data Transfer Agreement. **Do not attempt access through unapproved paths or redistribute them.**
- For example, previously reported candidate markers from this work must not be labeled as new Open-Disease discoveries.

**Schwalie et al. menstrual-fluid comparison data (2024).** [Public Zenodo release](https://doi.org/10.5281/zenodo.11105267), discussed in the 2026 BMC Medicine paper.

- Source reports 18 donors (7 disease, 11 controls), with an analysis subset of 5 cases and 7 controls after excluding hormonal contraceptive exposure; some donors contributed multiple samples/libraries. The library count is **not** independent donor n.
- Menstrual fluid is a different material from endometrial biopsy and should be used only for separately justified, clearly labeled cross-compartment work.

## Minimum input contract for exploratory pseudobulk

CSV columns, with exactly one row per study × donor × cell_type × gene:

```
study_id,donor_id,diagnosis,phase,hormone_treatment,cell_type,gene_id,counts,library_size,n_cells
```

- `diagnosis`: precisely `case` or `control`. Use the study's actual clinical definition in a separate provenance file.
- `hormone_treatment`: `none` only for baseline cohort v0; other exposure populations require a separate protocol (not automatic imputation).
- `phase`: an audited, nonmissing phase label consistent within each donor. Do not derive phase from disease label.
- `counts`: nonnegative **raw integer** gene counts aggregated from cells of that donor/cell type. Do not use log-normalized/integrated data or TPM as counts.
- `library_size`: sum of *full feature set* raw counts for that donor × cell type; constant within that sample. Include explicit zero rows for genes in the tested universe.
- `n_cells`: observed cell count, not independent-patient n.
- `study_id` and `donor_id`: stable keys. Confirm technical sample and participant overlap; same string in two cohorts is not sufficient to decide identity.

The code structurally validates data; **it cannot infer consent, sample independence, complete gene universe, batch/phase validity, anatomical comparability, or copyright licensing**. Record those checks by a human in a separately versioned study-level QC report.

## Data handling

1. Review publisher/Zenodo/GEO terms and individual consent/access limits. Seek institutional support where required.
2. Download source data into an access-controlled local environment, not this repository or CI.
3. Store original URL, version, size, SHA-256, retrieval date, local extraction script hash, reference genome/transcript IDs and QC exclusions. Use `python -m open_disease sha256 /path/to/file` locally.
4. Convert raw counts to sample-level tables; run `python scripts/audit_cohort.py --pseudobulk ... --out ...` before modeling.
5. Publish scripts, fixed non-identifying *aggregate statistics*, and provenance as permitted; **do not publish participant-level genotype, sequencing, clinical or identity metadata in GitHub**.
6. Maintain a discoverable statement of third-party licenses, patent disclosures, contributors' rights and permitted uses.

Do not interpret successful demonstration on artificial data as successful import of GEO/Zenodo human data.
