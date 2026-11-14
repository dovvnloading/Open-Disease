# GSE51981 cohort audit before expression analysis

**Decision on November 13, 2026 (US Eastern): NO-GO for the frozen donor-level model.** The 148 array records support a three-group, phase-filtered inventory of 144 candidate arrays. The reviewed sources do not establish the person-level key or individual recent-treatment status required for inclusion. Zero arrays were approved for fitting. This is a public-metadata assessment under [task #15](https://github.com/dovvnloading/Open-Disease/issues/15). It reports no expression outcome.

## Reconciled cohort

We matched every accession across the GEO Series Matrix header and the accession viewer’s sample-targeted SOFT quick view. Both representations agreed on platform GPL570, submitter source name, and all characteristic values for every record. They describe the same submission. Matching them checks transcription, alignment and field interpretation. The original article provides the separate published comparison below. [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981) · [Original paper, Table 1 and Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/)

| Original group | Proliferative | Early secretory | Mid-secretory | Late secretory in GEO | Unknown in GEO | Total arrays | Three-phase candidates |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Endometriosis | 29 | 18 | 28 | 0 | 2 | 77 | 75 |
| Other uterine/pelvic pathology, no endometriosis | 15 | 6 | 14 | 2 | 0 | 37 | 35 |
| No endometriosis or uterine/pelvic pathology | 20 | 6 | 8 | 0 | 0 | 34 | 34 |
| **Total** | **64** | **30** | **50** | **2** | **2** | **148** | **144** |

The nine cells for the three modeled phases match the published Table 1. Its footnotes call two additional cases and two additional other-pathology phases indeterminate. GEO calls the two out-of-scope case records unknown and two other-pathology records late secretory. The four fall outside the predeclared three phases under either source. This is an aggregate reconciliation; the paper's counts do not establish which individual records it designated indeterminate. Preserve both source descriptions. The primary comparison would have **109 candidate arrays** across the three modeled phases, before person-level and treatment review. All 37 other-pathology records retain their own group. [Table 1 with footnotes](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/)

| Decision field | Primary evidence and disposition |
| --- | --- |
| Array and platform | 148 distinct accessions on GPL570 matched across the two GEO formats. Each join used the accession, independent of file order. |
| Group and pathology | Used the explicit endometriosis characteristic for the 77 cases. For the 71 records explicitly labeled non-endometriosis, the separate pathology characteristic distinguishes the 34 pathology-free controls and 37 other-pathology samples. Specific comorbidities among cases cannot be recovered from this breakdown. |
| Phase and tissue | The submitter encodes phase in the **value of the tissue characteristic**. There is no separately named phase characteristic. The three modeled tissue labels explicitly describe endometrial tissue in 144 records. The study Methods describe archived endometrium across the cohort. We retained the two unknowns and the two out-of-scope late secretory labels as separate exclusions. |
| Disease severity | All 77 cases have a severity characteristic, but two also have a source name indicating that severity is unavailable. Do not use those two severity assignments for a stage comparison without reconciliation. Severity is outside the planned primary contrast. |
| Independent participant | The series design describes 148 women and the matrix contains 148 distinct arrays. This supports one sample per woman at the study-design level. Neither GEO representation supplies a named donor, subject, participant or patient characteristic for any array. Titles are specimen labels and descriptions are numbered sample labels. The source-to-person mapping required by the frozen model remains unresolved. |
| Hormonal treatment | The original study excluded women treated with hormones during the preceding three months. GEO has zero sample characteristics for hormone exposure or treatment. The exclusion criterion is relevant cohort evidence. Individual-level verification under the frozen manifest rule remains unresolved for all arrays. |
| Scale | All 148 individual GEO SOFT value descriptors identify base-2 logarithmic GC-RMA. The paper documents joint GCRMA normalization. No intensity values were decoded into an expression table, used to infer labels, or tested. |

[Original study, specimens and processing](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/) · [GEO sample-targeted metadata](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981&targ=gsm&view=quick&form=text) · [Frozen analysis rules](GSE51981_PREANALYSIS_PLAN_2026-11-12.md)

## Supplement and reuse review

The authors identify Supplemental Table 1 as the clinical-data supplement. The publisher indexes the workbook as `en-14-1490_stab1.xlsx`. The original link returned HTTP 403 in this environment. Its contents and any possible sample-to-participant crosswalk remain unverified. The later same-cohort analysis matches arrays by specimen title and applies additional severity, source and phase exclusions. Its code references a cleaned metadata CSV absent from the linked code revision. The ten public records in the associated Figshare project contain RDS artifacts. The cited secondary evidence does not resolve the missing participant key. [Original supplement link](https://press.endocrine.org/doi/suppl/10.1210/EN.2014-1490/suppl_file/en-14-1490_stab1.xlsx) · [Publisher index](https://academic.oup.com/endo/article/155/12/4986/2423578) · [Bunis et al. processing code](https://github.com/dtm2451/EndometrialDeconvolution/blob/37362434d5b2bbe569567d05a31bd0d72780cb48/microarray_processing.Rmd#L64-L105) · [Related Figshare project](https://figshare.com/projects/Whole-tissue_deconvolution_and_scRNAseq_analysis_identify_altered_endometrial_cellular_compositions_and_functionality_associated_with_endometriosis/127208)

EMBL-EBI's imported ArrayExpress study indexes 148 samples and 148 assays. Its indexed source characteristics cover organism, tissue, disease, severity and pathology. The index supplies no participant or treatment field. The SDRF download returned HTTP 503 here, so the contents of that file remain unverified. The public index response was 123,649 bytes, SHA-256 `6b058e6f61274a316199fff562abef3b87354349b493173929d0af7fdda21637` (2026-11-14 UTC). [EMBL-EBI study index](https://www.ebi.ac.uk/biostudies/api/v1/studies/E-GEOD-51981)

GEO describes public unrestricted access and does not impose use or distribution restrictions of its own. It also states that submitters and third parties may hold rights. The original study reports written informed consent for its collection. This audit used public GEO metadata locally and reports aggregate counts and source citations. It makes no determination of additional downstream rights to redistribute participant records, the clinical workbook, or expression files. [GEO data disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html) · [NCBI data-use policy](https://www.ncbi.nlm.nih.gov/home/about/policies/) · [Study consent statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/)

## Repeat the audit

The following opt-in acquisition transferred **677,247 bytes** of GEO data (response bodies) on November 14, 2026 UTC. The matrix server returned HTTP 206 and `Content-Range: bytes 0-131071/27099296`. Its modification timestamp was 2026-09-17 05:35:10 UTC. The bounded prefix contains the complete metadata header. Its digest identifies only this prefix. GEO's SOFT quick view embeds sample-table previews; the reconciliation skips each table body.

| Input | Bytes used | SHA-256 of those bytes |
| --- | ---: | --- |
| Official Series Matrix, compressed bytes 0–131071 | 131,072 | `fae8cdc59664d0a239c9d55120b0cd01154f8ec4d23cf3defc5ba099447d3321` |
| Official GSM-targeted SOFT quick response | 546,175 | `767499b710e33598cffcc7cdaf414925630ec241d90fcc04bcfec93b903bd5d8` |

Run from the repository root with Python 3.11 or newer. Keep both inputs in the ignored directory. Verify the digests before comparing results across runs.

```bash
mkdir -p data/private/gse51981-audit
cd data/private/gse51981-audit
curl --fail --silent --show-error --max-filesize 131072 --range 0-131071 \
  -o matrix-prefix.gz \
  'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51981/matrix/GSE51981_series_matrix.txt.gz'
curl --fail --silent --show-error --max-filesize 2400000 \
  -o gsm-quick.soft \
  'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981&targ=gsm&view=quick&form=text'
sha256sum matrix-prefix.gz gsm-quick.soft
cd ../../..
python scripts/reconcile_gse51981_sources.py \
  --matrix data/private/gse51981-audit/matrix-prefix.gz \
  --soft data/private/gse51981-audit/gsm-quick.soft
python -m pytest -q tests/test_reconcile_gse51981_sources.py
```

The deterministic reconciliation used CPython 3.12.14, 0.05 seconds wall time and 12.5 MiB peak child RSS on Linux. No random seed applies. The parser reuses the existing Series Matrix identity and platform guards from base commit `832b06a3a98f00c1e7293cf6d00a4440f63c6a17`. The corresponding script and synthetic regression tests accompany this note. The test data deliberately reorder records and fields, preserve unresolved categories, and reject duplicate identities, source disagreements, or malformed tables.

The remaining decision for [#18](https://github.com/dovvnloading/Open-Disease/issues/18) is whether a permitted primary record establishes the sample-to-person relationship and the individual treatment rule. A separately dated change to accept study-design evidence would change the protocol and the strength of the claim. This audit does not make that change.
