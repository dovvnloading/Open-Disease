# First GSE51981 reproduction attempt

**Outcome on November 13, 2026 (US Eastern) — blocked before model fitting.** Of 148 arrays, 144 fell within the three modeled phases. None could be assigned both an individually verified participant key and recent-treatment status under the frozen eligibility rules. No expression estimate, interval, p-value or permutation result exists for this attempt. [Cohort audit](GSE51981_COHORT_AUDIT_2026-11-13.md) · [Reproduction issue #18](https://github.com/dovvnloading/Open-Disease/issues/18)

## What is ready

| Requirement | Evidence |
| --- | --- |
| Published endpoint | The prospective commit selected FOSB probe set `202768_at`, its expected positive case–control direction, the other-pathology comparison, model, multiplicity and interpretation rules. [Frozen endpoint](https://github.com/dovvnloading/Open-Disease/commit/9cb999cf7bf3e714c6c53e398c48363d4c339fc8) · [Original Table 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/#T4) |
| Array annotation | The August 2016 GEO file has 54,675 distinct probe sets. The published audit gives its digest, release, complete classification and mapping policy. [GPL570 report](GPL570_ANNOTATION_RESULTS_2026-11-13.md) |
| Phenotype and phase | All 148 records reconcile across the two GEO representations. The three-phase inventory contains 75 cases, 34 pathology-free controls and 35 other-pathology arrays. [Source audit and individual-field limitations](GSE51981_COHORT_AUDIT_2026-11-13.md) |
| Processed data | Full GSE51981 compressed Series Matrix, 27,099,296 bytes, downloaded at 2026-11-14 02:46:06 UTC. SHA-256 `810e560c2b10dabb2f92dde03f64a2f9af4d72d51360f29c102569942dd84211`. Its gzip integrity check passed. All 148 sample value descriptions identify log2 GC-RMA. [Official source](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51981/matrix/GSE51981_series_matrix.txt.gz) |

## Reproduce the data-integrity and eligibility checks

After installing the project, download to an unused ignored path. Run the SOFT quick-view acquisition in the [cohort audit](GSE51981_COHORT_AUDIT_2026-11-13.md) before the reconciliation command below. The single complete download was within the 150 MB #18 budget. The header-only acquisition for #15 used a separate, bounded range request.

```bash
open-disease geo-download --out data/raw/GSE51981_series_matrix.txt.gz
gzip -t data/raw/GSE51981_series_matrix.txt.gz
echo '810e560c2b10dabb2f92dde03f64a2f9af4d72d51360f29c102569942dd84211  data/raw/GSE51981_series_matrix.txt.gz' | sha256sum --check
python scripts/reconcile_gse51981_sources.py \
  --matrix data/raw/GSE51981_series_matrix.txt.gz \
  --soft data/private/gse51981-audit/gsm-quick.soft
```

The full-file reconciliation reproduced the counts from the header prefix. The downloader and checksumming tools read bytes. Reconciliation read only metadata fields for the audit and used no expression values. All source files stayed local. The gzip digest covers the complete transfer, including expression bytes, and is distinct from the prefix digest in the cohort report. Software used for this check was CPython 3.12.14 and the repository downloader at baseline `832b06a`. No random seed applies to these operations. The synthetic model suite passed using R 4.5.0 and limma 3.62.2, including the exclusion of late-secretory records from the model.

## Evidence needed to resume

A permitted source connecting the intended array records to independent participants and supporting the recent-treatment rule would allow a fresh eligibility review. The paper documents a cohort-level treatment exclusion. Its clinical supplement was located but could not be verified. Reusing a specimen title as a person identifier would leave repeated donors undetectable. Keep the original study-design evidence and the unresolved individual fields visible. If evidence becomes available, record its version, reconcile eligible rows, and apply the prospectively fixed endpoint. Independent cohort validation remains outstanding.
