# Execution and reproducibility

The [first checkpoint](execution.json) records the original metadata-only audit. On 2026-11-16 the continuation used ordinary bounded public HTTP reads to get both primary workbooks (one per study), source metadata, reference sequence and aggregate regional statistics. All downloaded matrices and workbooks remained in ignored local results. No restricted data, participant records, enrollment, outreach, paid compute or changes to issue #49 were involved. There was no stochastic model and no seed.

| Computation | Input provenance and execution |
| --- | --- |
| Primary workbooks | [2023 source/size/digest and source cells](rahmioglu_workbook.json), [2026 workbook and source rows](koller_workbook.json). Python 3.12.14, openpyxl 3.1.5. The final UK scan took about 79 seconds, with peak RSS about 114 MiB. |
| Molecular and disease intervals | [GTEx v8 nominal data, official variant lookup, GRCh37/38 mapping and GRCh38 sequence](molecular_sources.json), [Catalog regional slice, 2026 manifest and README](disease_release_evidence.json). All summary-statistic slices were fetched by indexed HTTP ranges, with HTTP 206 extent and SHA-256 checks. No genome-wide file was downloaded. Per-gene QC was rerun byte-for-byte. |
| Pairwise tagging | [Original two-marker NIH LDlink response](pairwise_ld_input.json), [count-based calculation](pairwise_ld_result.json), [six Ensembl pairs](four_marker_ld.json) and [S31 proxy comparisons](s31_proxy_ld.json). The last two include the local freeze commit and the equivalent publicly archived design snapshot. The raw LD responses are aggregate counts or correlations. Public URLs may change output formatting across releases. |
| Accounting | Workbooks total 146,321,941 bytes. The three tabix indexes and six data ranges sum to 16,846,854 bytes (see source receipts for exact spans). The GRCh38 reference adds 2,000,003 bases; molecular APIs, manifests and targeted diagnostic reads add small ancillary transfers. The 250 MB task ceiling permitted the bounded inputs. These totals describe scientific inputs rather than all transport, dependencies, repository traffic or memory high-water marks across workers. |
| Inference and scope | Descriptive disease and QTL estimates, complete coverage by gene and reference LD were executed. There were zero fitted colocalization models, zero posterior prior-grid runs and zero multi-signal fits. The result on a changed disease input is exploratory. The source release and the original study share participants. |

## Replay

From the repository root, place the [2023](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01323-z/MediaObjects/41588_2023_1323_MOESM4_ESM.xlsx) and [2026](https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41588-026-02582-2/MediaObjects/41588_2026_2582_MOESM4_ESM.xlsx) publisher XLSX files under the two filenames used below. The readers enforce the original checksums. The two commands with an explicit fetch flag retrieve bounded public inputs when permitted by the operator and runtime. For pre-populated directories, omit the flag to verify local hashes without network access.

```sh
run_dir=results/gwas-3p21-colocalization-audit
evidence_dir=evidence/gwas-3p21-colocalization-audit
python scripts/gwas_3p21_rahmioglu_workbook.py "$run_dir/rahmioglu_2023.xlsx" --members-out "$run_dir/rahmioglu_s11_bsn_member_rsids.txt" | cmp - "$evidence_dir/rahmioglu_workbook.json"
python scripts/gwas_3p21_koller_workbook.py "$run_dir/koller-final.xlsx" | cmp - "$evidence_dir/koller_workbook.json"
python scripts/gwas_3p21_disease_slice.py "$run_dir/clinical-metadata" --fetch
python scripts/gwas_3p21_molecular_slice.py "$run_dir" --fetch
python scripts/gwas_3p21_sensitivity.py "$run_dir" | cmp - "$evidence_dir/sensitivity.json"
python scripts/gwas_3p21_pairwise_ld.py | cmp - "$evidence_dir/pairwise_ld_result.json"
python scripts/gwas_3p21_colocalization_audit.py | cmp - "$evidence_dir/verification.json"
python -m unittest discover -s tests -p 'test_gwas_3p21*.py'
python scripts/research_queue.py validate
python scripts/verify_catalog.py
python scripts/repo_safety.py
```

All 36 focused tests passed locally. Source hash verification succeeded for both workbooks, the regional disease and molecular inputs and the original coordinate/TWAS check. The next scientific comparison requires the exact clinical disease interval and independently measured, usable endometrial and myometrial QTLs. The [data inventory](inventory.md) identifies each unresolved input.
