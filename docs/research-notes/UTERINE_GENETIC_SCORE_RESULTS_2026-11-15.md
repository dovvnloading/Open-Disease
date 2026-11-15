# Does the direction predicted by uterine genetics appear in endometrial biopsies?

**Exploratory human-data computation, November 15, 2026.** The score fixed from published uterine expression genetics moved in opposite directions in the Finnish and Mumbai biopsies. Finland gave −0.60 percentile points and Mumbai +2.80. Neither passed the two-cohort decision rule. The result limits the claim that the genetically predicted directions can be read out as a common case–control expression shift in these studies. The inherited association itself remains a separate observation. [Pre-outcome specification](https://github.com/dovvnloading/Open-Disease/commit/ba37b1ae3518d0750a644de31f8fd779fd1434d8) · [Aggregate results](../../evidence/uterine-genetic-direction-2026-11-15/results.json)

## Which genetic prediction did the biopsies cover?

Koller et al. published 2,533 genes with a uterus-specific model in Supplementary Table 13. Seven passed the project’s uterus-only threshold of 0.05/2,533 (approximately 0.00001974) for the European-ancestry combined endometriosis definition. We aligned their signed z-scores to the published Ensembl gene IDs and the October 2026 HGNC snapshot. Under the frozen common-coverage rule, the Rome array left just three genes, all pointing toward lower expression with higher genetically predicted susceptibility. All three fell in the same cytoband, 3p21.31. The primary score was therefore the negative mean of their measured expression ranks. [Primary genetic source and Supplementary Table 13](https://doi.org/10.1038/s41588-026-02582-2)

| HGNC band | Gene | Uterus z | Gene included in the common score |
| --- | --- | ---: | --- |
| 3p21.31 | RBM6 | −5.19 | Yes |
| 3p21.31 | TMEM115 | −4.42 | Yes |
| 3p21.31 | HYAL3 | −4.41 | Yes |
| 3p21.31 | GMPPB | +5.21 | Absent from Rome mapping |
| 4q24 | CISD2 | +4.78 | Absent from Rome mapping |
| 8q21.13 | CRISPLD1 | +6.04 | Absent from Rome mapping |
| 15q15.1 | EIF2AK4 | −4.29 | Absent from Rome mapping |

HYAL3 has a separate gene-specific result. In the European clinical-diagnosis analysis, the authors reported colocalization with rs13100173 in uterus (full-colocalization posterior 0.817, candidate SNP contribution 0.721) and five additional tissues. It was the only uterus row in their positive expression-colocalization table. This analysis used a different case definition from the signed TWAS. The colocalization probabilities do not orient a risk allele or distinguish the causal tissue. GTEx sought uterine samples containing myometrium and endometrium. When evaluable, the endometrium was typically inactive basalis. This differs from the cycling endometrial biopsies measured here. [Supplementary Table 8](https://doi.org/10.1038/s41588-026-02582-2) · [GTEx histopathology](https://doi.org/10.5858/arpa.2023-0468-OA)

## Effects in whole biopsies

Positive values denote a measured shift toward the genetic risk direction. The unit is one point on the signed within-biopsy percentile scale. Each cohort retains its own rank reference. The 95% intervals invert the exact label-allocation test under a constant additive shift and exchangeability in the stated design.

| Cohort | Cases / controls | Case minus control, points | 95% interval | One-sided tail in the genetic direction | Holm p |
| --- | ---: | ---: | --- | ---: | ---: |
| Finland, seven matched run-by-phase strata | 12 / 13 | −0.599 | −1.200 to +0.108 | 11,522 / 12,096 | 0.953 |
| Mumbai, mid-proliferative | 4 / 3 | +2.804 | −1.602 to +11.891 | 3 / 35 | 0.171 |
| Rome, descriptive, two mixed scan dates | 5 / 6 | −5.416 | −7.548 to −2.817 | 199 / 200 | — |

All three raw gene contrasts in Mumbai were negative. In Finland, RBM6 and HYAL3 ranked higher in cases and TMEM115 ranked lower. The Finnish proliferative-only contrast was −0.723 points (8 cases, 4 controls). In Rome, all controls were scanned before any case on each mixed date, making scan chronology a competing explanation. A two-sided Rome tail of 2/200 under date-only reallocation does not resolve the source of that difference. Participant overlap with the genetic meta-analysis is unknown, and genetic liability and expression after disease develops describe different processes. Recruitment, tissue composition, menstrual timing and technical handling may all contribute to the measured contrasts. [Rome scan audit](GENOME_TRANSFER_RESULTS_2026-11-15.md)

## Secondary analyses after the coverage result

Finland, Mumbai and the Swiss glandular-cell dataset measured all seven genes. After seeing the shared-panel results, we retained all seven conditional on complete coverage, averaged signed ranks within cytobands, then gave each of the four bands one-quarter weight. Equal band weights stop the four genes at 3p21.31 from dominating the score. Cytobands provide a physical grouping only. [Secondary analysis specification](UTERINE_GENETIC_SCORE_COVERAGE_2026-11-15.md)

| Cohort | Seven-gene score, points | 95% interval | One-sided allocation p | Omit-one-band range |
| --- | ---: | --- | ---: | --- |
| Finland | −1.398 | −3.655 to +0.600 | 0.904 | −2.297 to +0.298 |
| Mumbai | +1.793 | −7.127 to +7.679 | 0.229 | +1.042 to +2.538 |
| Swiss proliferative glandular epithelium | +0.601 | −0.444 to +2.213 | 0.125 | −0.149 to +1.126 |

The whole-tissue pair still disagreed in sign. Secondary Holm p-values were 0.904 in Finland and 0.457 in Mumbai. The Swiss value came from six affected and seven comparator donors in four capture batches, using SCT-corrected aggregates. The original three-gene score was +0.730 points there (95% interval −4.798 to +8.481). A study that can distinguish a genetic tissue mechanism would measure the risk allele and allele-aligned expression in the same, cycle-dated uterine compartment, with batch-balanced samples. Disease-status contrasts alone do not supply that test.

## Reproduction

The [implementation](../../scripts/uterine_genetic_direction.py) verifies the genetic workbook, HGNC release, original Finnish native assay and archived matrices through pinned loaders. Its gene-first calculation reproduces the original Finnish per-gene differences for all three primary genes. An independent rational enumeration checks the directional tail, two-sided tail and observed contrast. Rereading the complete original RDS reproduced the Swiss five-pathway sentinel before either new score was evaluated. The repeat run with the secondary panels reproduced all fields of the first primary result. No random seed applies. Source tables are public aggregates. The genetic workbook SHA-256 is `c021b6c20c668644939fd83381307a2c8cf08f18d95156b33e65fb0301997815`.

With the [Finnish archive and manufacturer manifests](GSE141549_CHIP_CENSUS_2026-11-14.md), the [Rome and Mumbai sources](GENOME_TRANSFER_RESULTS_2026-11-15.md), the [Swiss RDS, workbook and HGNC release](SWISS_GLANDULAR_RESULTS_2026-11-14.md), and the [Koller workbook](../../evidence/locked-genetics-2026-11-15/genetics_sources.json) installed locally, run the following from the project root. The Python environment also needs openpyxl and inflate64 1.0.4. The cellular comparison uses base R.

```sh
d=data/private/gse266265
OPENBLAS_NUM_THREADS=1 python scripts/uterine_genetic_direction.py \
  data/private/genetics-lookup/koller.xlsx data/private \
  --swiss "$d/GSE266265_epithelial_cells_monocle3.rds" "$d/41467_2026_73020_MOESM3_ESM.xlsx" \
  > results/uterine-genetic-direction.json
```
