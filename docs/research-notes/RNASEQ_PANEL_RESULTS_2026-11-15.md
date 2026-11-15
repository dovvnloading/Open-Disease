# Five-gene lesion panel in paired Peking ovarian RNA sequencing

**Exploratory human-data computation, 15 November 2026.** The same five genes selected using the ovarian microarrays changed in their locked direction in **every one of five new donor pairs**. The signed panel was **42.94 within-sample percentile points higher** in ovarian lesions than in paired uterine endometrium. Only one of the 32 possible within-donor tissue-label assignments yielded an equally large or larger score (prespecified one-sided p = 0.03125, two-sided p = 0.0625). The five donors determine the test's resolution. Removing any one pair gave a mean shift of +41.26 to +44.99 points. [Local follow-up specification](RNASEQ_PANEL_PLAN_2026-11-15.md) · [Original panel lock](../../evidence/lesion-transfer/panel.json) · [Aggregate results](../../evidence/lesion-transfer/peking-rnaseq.json)

| Gene | Discovery direction in lesions | Peking change, percentile points | Pairs in that direction | Deposited transcripts |
| --- | --- | ---: | ---: | ---: |
| TMEM255A | Higher | +29.98 | 5/5 | 5 |
| DIRAS3 | Higher | +65.93 | 5/5 | 2 |
| GABRP | Lower | −33.46 | 5/5 | 9 |
| SALL1 | Lower | −25.23 | 5/5 | 4 |
| JCHAIN | Higher | +60.09 | 5/5 | 6 |

Each gene's one-sided p was 1/32 and its five-gene Holm-adjusted p was 0.15625. The 25 gene-by-donor comparisons reuse the same five people. Ranking transcripts within each sample first, then averaging their percentiles within each gene, gave a smaller panel shift of **+27.31 points**, positive in all five pairs (one-sided p = 1/32). [All estimates](../../evidence/lesion-transfer/peking-rnaseq.json)

## Provenance and pairing

The GSE282532 study was submitted by Peking University First Hospital. Its BioProject documents surgery, five participants with ovarian endometriomas, a cyst-wall and uterine endometrium biopsy from each participant, regular menstrual cycles, no hormone therapy and proliferative-phase collection. The associated publication describes separate patient enrollment and institutional ethics approval. Recruitment at this hospital and the 2024 submission distinguish it from the earlier separately collected microarray cohorts and the older Chinese PLA General Hospital sequencing cohort. The ten GEO records here correspond to the five sequenced pairs. Other assays described in the paper involved different sample counts. [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1189111) · [Source publication](https://doi.org/10.1038/s41467-026-69551-5) · [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282532)

Each workbook column matched exactly one GEO sample description. Each tissue annotation agreed with both the title prefix and assay alias, and the numbered source titles and paired library suffixes independently gave the same five pairs. The source release contains 83,064 uniquely identified transcript rows mapping to 19,688 Ensembl genes and supplies FPKM, a length- and depth-normalized abundance measure. All five locked Entrez IDs mapped uniquely to approved HGNC entries, and the corresponding Ensembl gene IDs and symbols matched the workbook. The reference distribution retained all measured genes and zero values. [Series metadata](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282532/matrix/GSE282532_series_matrix.txt.gz) · [Original workbook](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE282nnn/GSE282532/suppl/GSE282532_mRNA_FPKM.xlsx)

We summed deposited transcript FPKMs within a gene, assigned percentile midranks across genes separately within each specimen, subtracted the same donor's uterine rank from the lesion rank, and averaged the five genes with the frozen signs. A second implementation read the workbook XML directly, independently summed the original decimal values and counted genes below or tied with each target to verify the primary percentiles. It checked every gene's mean shift and the panel's exact directional tail. [Analysis](../../scripts/rnaseq_lesion_panel.py) · [Independent calculation](../../scripts/verify_rnaseq_lesion_panel.py) · [Verification receipt](../../evidence/lesion-transfer/peking-verification.json)

## Interpretation

The panel now has independent positive aggregate contrasts in the Hull series (7/9 people), CodeLink series (11/11) and this sequencing series (5/5). All three test the same fixed gene identities and directions. The Peking data extend the observation to sequenced, proliferative-phase ovarian endometriomas. The different tissue compositions of lesions and uterine biopsies, the distinction between ovarian and extraovarian lesions and the use of different platforms limit comparison of effect sizes. Source recruits had established disease. No disease-free group, diagnostic threshold, clinical prognosis or treatment outcome was tested. [Discovery and Hull transfer](LESION_TRANSFER_RESULTS_2026-11-15.md) · [CodeLink test](CODELINK_PANEL_RESULTS_2026-11-15.md)

## Reproduction

Place the linked source workbook and series matrix in the ignored data directory shown below. Use the October 2026 source-pinned HGNC complete set from the repository's existing Swiss workflow. All input paths and checksums are validated. The local plan had SHA-256 `8517398e766c2b1ba6d1056ca7135716ef6e2cdd7904d7262f47952fc753756d` before expression analysis. Existing validations had been seen when this separate follow-up was chosen. [Source pins](../../evidence/lesion-transfer/peking-rnaseq.json)

```sh
mkdir -p results/lesion-transfer
python scripts/rnaseq_lesion_panel.py data/private/gse282532 \
  data/private/gse266265/hgnc_complete_set_2026-10-06.tsv \
  > results/lesion-transfer/peking-rnaseq.json
cmp results/lesion-transfer/peking-rnaseq.json evidence/lesion-transfer/peking-rnaseq.json
python scripts/verify_rnaseq_lesion_panel.py data/private/gse282532/GSE282532_mRNA_FPKM.xlsx \
  results/lesion-transfer/peking-rnaseq.json
```
