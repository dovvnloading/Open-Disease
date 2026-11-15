# A paired ovarian lesion expression profile transfers to extraovarian lesions

**Exploratory human-data computation, 15 November 2026.** We learned the direction and size of expression differences between ovarian lesions and the same patient's uterine endometrium in ten people. Applied unchanged to a separately collected series, the profile shifted in the predicted direction in **eight of nine** people with peritoneal or broad-ligament lesions. The mean signed score difference was **+6.10 percentile points**. Only six of 512 possible within-person label assignments produced a score at least as large (one-sided p = 0.0117). A smaller, discovery-selected five-gene panel shifted by +22.18 points in the second series (7/9 pairs, p = 0.0371). [Local analysis specification](LESION_TRANSFER_PLAN_2026-11-15.md) · [Panel lock](../../evidence/lesion-transfer/panel.json) · [Complete results](../../evidence/lesion-transfer/results.json)

## What transferred across anatomical sites

All 19 participants had endometriosis. Both biopsies in each comparison came from the same person. Hever et al. compared ovarian lesions with uterine endometrium in GSE7305. The human series from Hull et al. comprised seven broad-ligament and two peritoneal lesions and their paired uterine biopsies in GSE11691. We used the native tissue and donor fields to establish pairing and required a complete set of distinct feature identifiers on each platform. [Hever original study](https://pmc.ncbi.nlm.nih.gov/articles/PMC1941489/) · [Discovery deposit](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7305) · [Validation deposit](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11691)

For each sample, rank all array features, assign midranks to ties and convert the ranks to percentiles. Retain each uniquely annotated Entrez target present on both platforms, excluding the five genes already examined in the genetic-anchor analysis. Average eligible probes within a gene. For each of the resulting 12,498 genes, calculate its average within-donor lesion-minus-uterus difference in the discovery series. The signed weight of each gene is that difference divided by the sum of all absolute differences. Apply those weights to the held-out nine pairs. The resulting score measures the extent to which genes higher or lower in ovarian lesions move in the same directions in extraovarian tissue. Its unit is the signed weighted change in within-array percentile points. It cannot be converted to a fold change. [Methods and parser](../../scripts/lesion_transfer.py) · [Aggregate results for all genes](../../evidence/lesion-transfer/genes.tsv)

| Validation analysis | Matched pairs in predicted direction | Mean score change | Exact one-sided p |
| --- | ---: | ---: | ---: |
| All 12,498 genes, primary | 8/9 | +6.10 points | 6/512 = 0.0117 |
| Same physical probes and rank reference on both platforms | 8/9 | +6.56 points | 6/512 = 0.0117 |
| Remove 125 most heavily weighted genes, renormalize | 8/9 | +5.14 points | 5/512 = 0.00977 |
| Score only the sign of each validation gene change | 8/9 | +0.282 on a −1 to +1 scale | 3/512 = 0.00586 |
| Frozen five-gene panel | 7/9 | +22.18 points | 19/512 = 0.0371 |

Swapping every gene together within each donor preserves cross-gene correlations. The exact test assumes the paired differences are sign-exchangeable under its null; tissue was observed rather than randomly assigned. Discarding magnitude and testing only whether each donor's primary score rose gave two-sided sign-test p = 0.0391. The primary effect ranged from +5.60 to +7.92 points after dropping each validation pair. Relearning the weights after omitting each discovery pair gave +5.78 to +6.16 points. [Aggregate sensitivity analyses](../../evidence/lesion-transfer/results.json)

The cohort mean effects had Pearson correlation 0.476 and Spearman correlation 0.415. Of 12,498 genes, 8,090 (64.7%) had matching signs. Restriction to 19,917 identically mapped shared probe sets left all 12,498 genes, raising the respective correlations to 0.505 and 0.440. These correlations and the 8,090 signs describe the two data sets. Correlated genes do not supply 12,498 independent biological replicates.

## Where the five-gene projection succeeded and failed

We locked two largest increases and two largest decreases in the discovery contrast, then the remaining largest absolute shift. Every gene received the same weight after orienting its sign to discovery. The lock was written before this analysis read validation expression. These public cohorts and prior analyses were accessible. No prospective clinical validation took place.

| Gene | Ovarian lesions minus uterus, mean points | Extraovarian lesions minus uterus, mean points | Extraovarian pairs in discovery direction | Directional p | Five-gene Holm p |
| --- | ---: | ---: | ---: | ---: | ---: |
| TMEM255A | +78.91 | +12.38 | 7/9 | 0.1523 | 0.3047 |
| DIRAS3 | +75.68 | +3.62 | 7/9 | 0.3184 | 0.3184 |
| GABRP | −69.29 | −29.46 | 8/9 | 0.0586 | 0.1758 |
| SALL1 | −61.92 | −30.14 | 7/9 | 0.0313 | 0.1250 |
| JCHAIN | +71.52 | +35.27 | 8/9 | 0.0117 | 0.0586 |

All five cohort means agreed in direction; all five shrank in magnitude. Any biological difference between ovarian and extraovarian lesions is confounded here with study and platform. Selecting the largest discovery changes also favors overestimation. The sign test for the smaller panel was p = 0.180 (two-sided). Its exact magnitude-weighted directional p is secondary to the single genome-wide primary endpoint. The five individual tests each missed their Holm threshold. [Selection rule and results](../../evidence/lesion-transfer/results.json)

In the independent CodeLink study, the same signed panel rose in 11 of 11 pairs (+24.91 percentile points, one-sided p = 1/2048). The ten pairs whose native title suffixes match exactly gave +26.66 points, 10/10 positive. The remaining pair was linked by the curated GEO DataSet. Four of the five directional gene tests passed within-panel Holm correction in this cohort. [CodeLink study and source audit](CODELINK_PANEL_RESULTS_2026-11-15.md)

In a third external study, all five genes moved in their locked direction in each of five Peking ovarian lesion–uterus pairs measured by RNA sequencing. The panel rose by +42.94 points, one-sided p = 1/32. [RNA-seq study and verification](RNASEQ_PANEL_RESULTS_2026-11-15.md)

For symmetrical evidence about single genes, we also tested a same-direction association jointly across the two cohorts. Each gene's joint p is twice the smaller of the two directional maxima across studies, capped at one. The minimum possible joint p was 2/512 = 0.00391. With 12,498 hypotheses, the first Holm threshold was 0.00000400. Even complete directional consistency cannot identify a genome-wide Holm-significant gene in this design. [All joint and Holm p-values](../../evidence/lesion-transfer/genes.tsv)

## Prior evidence and interpretation

Cross-study lesion expression is well-established. A 2009 analysis included these studies and GSE5108 and identified pathways shared between lesion data sets. A 2021 analysis overlapped differentially expressed genes across these and a fourth series. A later study reported higher JCHAIN and lower GABRP in lesions. Independent laboratory work localized elevated DIRAS3 predominantly to lesion epithelium. The original Hever report documented plasma cells in ovarian lesions, which offers one explanation for immunoglobulin-related expression. [2009 pathway analysis](https://pmc.ncbi.nlm.nih.gov/articles/PMC2752458/) · [2021 expression overlap](https://pmc.ncbi.nlm.nih.gov/articles/PMC8381947/) · [JCHAIN and GABRP](https://pmc.ncbi.nlm.nih.gov/articles/PMC10616016/) · [DIRAS3 tissue and cell study](https://pubmed.ncbi.nlm.nih.gov/37598541/) · [Hever study](https://pmc.ncbi.nlm.nih.gov/articles/PMC1941489/)

This contribution estimates the size and reproducibility of a profile learned in one anatomical site and fixed before testing in another with the donor as the inference unit. Different cell mixtures, local organ tissue, inflammation, fibrosis, tissue handling and menstrual phase can contribute. A matched uterine biopsy controls for person-level variation but cannot distinguish these explanations. There are no disease-free lesion-site tissue controls here. The bulk comparison cannot identify a cause of disease or establish whether a uterine sample can diagnose endometriosis. A cellular follow-up can test how much of the signed pattern persists within epithelial and stromal compartments.

## Reproduction

Use the source-pinned GSE7305 and GSE11691 series matrices plus the GPL570 and GPL96 annotation files under the private data directory. Their SHA-256 values are enforced by the data reader. These are GEO's series matrix and platform annotation downloads, not raw CEL files. The training command reads validation platform annotations only. The scoring command refuses a lock whose discovery genes, directions or effects cannot be reproduced. Input series and platform files remained outside version control. The analysis is deterministic and was run with Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. [Pinned input checksums](../../scripts/genetic_anchor_tissue.py) · [Script](../../scripts/lesion_transfer.py) · [Tests](../../tests/test_lesion_transfer.py)

```bash
mkdir -p results/lesion-transfer
export OPENBLAS_NUM_THREADS=1
python scripts/lesion_transfer.py lock data/private/genetics-expression \
  --lock results/lesion-transfer/panel.json
cmp results/lesion-transfer/panel.json evidence/lesion-transfer/panel.json
python scripts/lesion_transfer.py validate data/private/genetics-expression \
  --lock evidence/lesion-transfer/panel.json --genes results/lesion-transfer/genes.tsv \
  > results/lesion-transfer/results.json
cmp results/lesion-transfer/results.json evidence/lesion-transfer/results.json
cmp results/lesion-transfer/genes.tsv evidence/lesion-transfer/genes.tsv
python -m pytest -q tests/test_lesion_transfer.py
```

The local specification had SHA-256 `9b3b11f1045a19da30502a633c0525039deb32fe5b096132a570eb66c39c181a`. The five-gene panel lock had SHA-256 `4c634a65d648a83af8418e7ad1bf50e03d1a83982d0265d47d9ccbd06aba0114`. The primary endpoint had a prespecified directional alpha of 0.05 within this analysis. Estimates for the panel and all sensitivities share donors with the primary endpoint.
