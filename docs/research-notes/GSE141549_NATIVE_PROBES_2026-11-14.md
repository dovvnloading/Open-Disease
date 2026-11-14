# Physical FOSB probe provenance in the native EndometDB release

**November 14, 2026 UTC · Literature and metadata assessment · Issue #26 remains open.** The native EndometDB raw table preserves enough information to identify HumanWG-6 v2 on two complete assays and link both to the authors’ own clinical assay table. Neither assay belongs to the 38 donors with exactly one eligible analysis result. The Finnish case–control comparison remains unestimated. [Versioned author release](https://doi.org/10.6084/m9.figshare.12195975.v1) · [Clinical eligibility](GSE141549_ELIGIBILITY_2026-11-13.md)

## Native physical measurements

The dump explicitly declares the foreign key `microarray_data.sample_analysis_result_id → sample_analysis_result.id`. The earlier eligibility audit documents the remaining sample-to-patient joins within this release. The raw table uses integer `probeid` values. In each of two complete consecutive assays, the set of 48,701 addresses equals the full WG-6 manufacturer manifest, with exactly one record per address. Of those addresses, 11,543 appear in neither HumanHT-12 v4 manifest. A third assay is only partly present in the byte range and contributes no chip assignment. [Author release](https://api.figshare.com/v2/articles/12195975/versions/1) · [Illumina WG-6 manifest](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13376/suppl/GPL13376_HumanWG-6_V2_0_R4_11223189_A.txt.gz)

| Native integer address | WG-6 manufacturer probe | HT-12 manufacturer probe, both revisions |
| ---: | --- | --- |
| 6040086 | **ILMN_1730706 (FOSB)** | ILMN_1685043 (CYP3A7) |
| 7160239 | ILMN_1718671 (GJB6) | **ILMN_1751607 (FOSB)** |

Both addresses exist on both chips. Only the platform-specific pairing identifies FOSB. The WG-6 and HT-12 FOSB probes have distinct 50-base sequences. [HT-12 manifest R1](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/suppl/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz) · [HT-12 manifest R2](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/suppl/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz) · [Platform comparison](GSE141549_ELIGIBILITY_2026-11-13.md)

## Limits of the normalized export

The separate `normalized_microarray_data` COPY stream declares `sample_analysis_result_id, avg_signal, search_key, symbol, probeid, array_address_id`. Its first 70,503 complete records span 169 reported probe labels and 355 native assay IDs. In this prefix every search key and array address is 0. There are only 53,770 distinct combinations of the five nonmeasurement fields. The first probe label, ILMN_1802380, occurs 864 times across 334 assay IDs, with one to eight records per ID. FOSB is absent from this prefix. These observations establish that joining the normalized table on assay and reported probe can return several rows. FOSB-specific multiplicity and the reason for repetitions remain unknown. The audit counted annotations without converting or comparing signal values. [Versioned archive](https://doi.org/10.6084/m9.figshare.12195975.v1)

The previous complete annotation scan located FOSB in combined GEO workbook row 21,565 under the HT-12 label. The three members of the GEO RAW tar are the manufacturer manifests, as listed in the depositor’s file inventory. The authors’ GitHub repository supplies the web application, and their paper says that data-processing code is available upon request. We found no published per-chip derivation for that row or author-supplied GEO-to-native assay crosswalk in these materials. [Pooled workbook audit](GSE141549_FOSB_MEASUREMENT_2026-11-13.md) · [GEO tar inventory](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141549/suppl/filelist.txt) · [Author code](https://github.com/micawo/EndometDB/tree/5465ed005d0210256fb0e22e0c544f037568b460) · [Paper and code availability](https://doi.org/10.1038/s41597-020-00623-x)

## What remains to run a Finnish test

A native raw-data comparison can stay within the declared foreign keys. It requires verifying chip identity and address coverage for the **eligible** native assays, recounting donor overlap by run and phase on a single chip, and freezing raw-data preprocessing before examining FOSB signals. The earlier ceiling across chip generations was 12 cases and 13 controls in seven shared run-by-phase groups. The survey has no treatment lookback interval, and the control definition is healthy women undergoing tubal ligation. No cross-release matching was attempted. The raw PG table occupies approximately 512.6 MB of the dump before ZIP compression. A sequential download of this complete table exceeds the issue’s 200 MB allowance. The two complete assays in the source prefix establish a feasible identity check but do not establish Finnish outcome eligibility. [Clinical eligibility and control definitions](GSE141549_ELIGIBILITY_2026-11-13.md) · [Task budget](../../research/tasks/endometdb-fosb-measurement-audit.json)

## Reproduce the bounded audit

Use the existing pinned metadata head, tail and three manufacturer manifests listed in the [eligibility audit](GSE141549_ELIGIBILITY_2026-11-13.md). Add the two ranges below from Figshare file 22809185 to the same ignored local directory. The combined additional archive transfer was 5,500,001 bytes. The program verifies all seven input hashes, checks the native foreign key, inflates the release-specific PG data prefixes, and prints aggregate counts. It never parses quantitative columns. Python 3.12.14 and inflate64 1.0.4 were used. No random procedure was involved.

| Inclusive ZIP byte range | Local filename | SHA-256 |
| --- | --- | --- |
| 65,536–2,665,535 | native-raw-continuation.bin | `90cbe73c4fd51f102f73a7bfdf834597c0bf7f6b50eb9c2f8ec577d2abdd58e6` |
| 510,400,000–513,300,000 | native-normalized-head-range.bin | `a62f3374bf39d6b460994089298fec0a90dc8049757c4565504c782c063c129a` |

```bash
cd data/private/gse141549
url=https://ndownloader.figshare.com/files/22809185
curl --fail --location --range 65536-2665535 --max-filesize 2600000 "$url" -o native-raw-continuation.bin
curl --fail --location --range 510400000-513300000 --max-filesize 2900001 "$url" -o native-normalized-head-range.bin
cd ../../..
python scripts/audit_endometdb_probe_keys.py data/private/gse141549
python -m pytest -q tests/test_audit_endometdb_probe_keys.py
```
