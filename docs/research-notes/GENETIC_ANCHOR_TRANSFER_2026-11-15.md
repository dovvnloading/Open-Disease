# Do inherited endometriosis signals identify a reproducible uterine RNA signature?

**Exploratory analysis of human-cohort measurements and published aggregate genetics, November 15, 2026.** An externally defined genetic panel produced a strong case–control difference in Rome. The other 12,029 assayed genes showed a larger average difference. Conditioning on physical scan date preserved the panel association and widened the gap in favor of the background. These data provide little evidence that disease genetics specifically explains the Rome transcriptome contrast.

## Inherited susceptibility and measured RNA

Rahmioglu et al. integrated disease GWAS with expression and methylation QTLs. We took all 23 distinct named genes in the middle column of their Table 1, across both tissues and molecular traits. Fifteen were measurable in the repository’s shared Finnish–Rome universe. Koller et al. separately predicted genetically regulated uterine expression using TWAS. Their publisher workbook explicitly labels Table S13 as the European-ancestry combined endometriosis definition. Applying its stated p < 3.72 × 10⁻⁶ threshold to all 2,533 nonmissing uterus models among 13,457 rows identified exactly four genes. RBM6 expression was genetically associated with lower disease risk and GMPPB, CRISPLD1 and CISD2 with higher risk. A GWAS/TWAS is a susceptibility analysis. A difference in tissue from people who have disease can incorporate its consequences and differences in tissue composition. Overlap among the genetic consortia precludes independent-genetics replication claims. [2023 primary study, Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10042257/#T1) · [2026 primary study](https://doi.org/10.1038/s41588-026-02582-2) · [2026 publisher workbook, S13](https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41588-026-02582-2/MediaObjects/41588_2026_2582_MOESM4_ESM.xlsx)

The panel was defined from published studies. This analysis was designed retrospectively after the Finnish/Rome exploratory results, including the known ABO result. We recorded a local method before calculating these panel outcomes. The competitive, scan-day and additional cohort checks were designed after seeing the initial Rome results. No prospective external holdout was registered. [Local analysis design](../../evidence/genetic-anchors-2026-11-15/design.md) · [Follow-up order](../../evidence/genetic-anchors-2026-11-15/followup.md)

## Why the Rome result looks convincing in isolation

For each measured genetic gene we squared its case–control difference in tissue RNA percentile, divided by its variance under cohort label permutation and averaged across genes. This tests for changes in either direction. Every row receives equal weight. Every permutation moves whole donor labels, preserving dependence between genes. The competitive check subtracts the corresponding average for all other eligible genes in the same permutation. The conditional Rome analysis gives equal weight to the two scan dates containing both diagnoses (three cases and three controls on one, two cases and three controls on the other). Three additional cases shared a day without controls.

| Cohort and design | Measured genetic genes | Genetic mean z² | Other genes' mean z² | Genetic panel vs no disease shift | Genetic panel vs other genes |
| --- | ---: | ---: | ---: | ---: | ---: |
| Rome, 8 cases and 6 controls | 15 | 3.353 | 3.549 | 2/3,003 (0.000666) | 2,221/3,003 (0.740) |
| Rome, within two shared scan dates | 15 | 2.652 | 3.217 | 2/200 (0.010) | 198/200 (0.990) |
| Mumbai, 4 cases and 3 controls | 14 | 1.008 | 0.766 | 14/35 (0.400) | 8/35 (0.229) |

Mumbai lacked a qualifying FSHB measurement. Its background comprised 11,450 genes, against 12,029 in either Rome comparison. In unblocked Rome, removing ABO gave 1/3,003 and the four measured endometrium-eQTL genes alone also gave 1/3,003. The four endometrium-eQTL gene measurements and the nearby gene detections are correlated. These are descriptive sensitivity estimates. Six of 15 panel genes had matching directions in Finland and Rome. The remaining nine reversed. For each of the 23 genes we required the same sign and evidence in both cohorts; all Holm-adjusted conjunction p-values were 1 for assayed genes. Eight were unassessed. Eleven of 23 had at least one exact source-symbol probe in Melbourne's published high-detection panel. [All 23 genes and coverage](../../evidence/genetic-anchors-2026-11-15/rahmioglu_expression.tsv) · [Specificity and Mumbai](../../evidence/genetic-anchors-2026-11-15/transfer.json) · [Shared-day findings](../../evidence/genetic-anchors-2026-11-15/scan_day.json)

ABO, Finland's leading gene, was previously implicated by a blood eQTL in the 2023 Table 1 (rs550057). In 2026 the clinical-definition GWAS colocalized with four blood methylation sites upstream of ABO at rs635634 (both reported posterior probabilities 1). These prior results establish a literature link missed in the project’s original blood-type discussion. Finland's +18.36 percentile points and Rome's +0.79 points remained inconclusive on external comparison. This work does not establish a uterine ABO RNA mechanism. [2023 primary study](https://doi.org/10.1038/s41588-023-01323-z) · [2026 results](https://pmc.ncbi.nlm.nih.gov/articles/PMC13195585/)

## An inherited direction for RBM6, and conflicting biopsy results

Only RBM6 of the four TWAS genes passed the project’s common chip mapping. Although every unblocked Rome disease label aligned with lower RBM6 rank (1/3,003 one-sided), it ranked 2,982nd of 12,044 genes by absolute standardized contrast. The exact permutation tail is determined by ordering relative to each gene’s own null. Its genome-wide standardized rank addresses a different question. Within mixed scan days RBM6 was −10.17 percentile points (1/200 one-sided). The Holm-adjusted primary pair of p-values was 0.000666 each in unrestricted Rome, and 0.010 each within shared dates. Discovery and technical follow-ups reuse the same women.

| Cohort | Measurement | Case minus control | Interpretation |
| --- | --- | ---: | --- |
| Finland | Shared-chip RBM6 RNA rank | +1.20 pp | Opposite the lower-RBM6 susceptibility direction |
| Rome | RBM6 RNA rank | −7.83 pp | Within mixed scan dates −10.17 pp |
| Mumbai | RBM6 rank of summed gene FPKM | −2.85 pp | One-sided 10/35 = 0.286 |
| GSE315857 | Mean of five measured RBM6 transcript ranks | +8.35 pp | Opposite the Rome estimate |
| Melbourne | Published cycle-adjusted high-detection probe | Direction unresolved | Two-sided p = 0.456, Holm over four genes = 1 |

All four primary gene-total ranks in Mumbai pointed in the genetic directions. The direction-weighted mean case–control shift was +0.428 within-cohort gene SD units (6/35 one-sided, p = 0.171, Holm for the two primary Mumbai tests = 0.343). Ranking transcripts before averaging within genes reversed GMPPB and CRISPLD1 and reduced the score to +0.027 (16/35, p = 0.457). In GSE315857 three of four targets mapped unambiguously via the frozen GENCODE v19 / HGNC annotation, with no qualifying CISD2 transcript. The signed score was −0.397. It stayed below zero after any one withdrawal (−0.503 to −0.270, fixed gene scaling). Library-name prefixes were uneven across diagnoses, and the RNA abundance units are unresolved. We report no inferential p there. Melbourne covered three of the four targets; none passed Holm. Adult uterine expression differs with cycle, treatment and cell composition. Technical and biological heterogeneity, timing, and the small cohorts prevent a general claim of disease-induced change.

[Four genetic genes and complete coverage](../../evidence/genetic-anchors-2026-11-15/koller_uterus.tsv) · [Eight-biopsy result](../../evidence/genetic-anchors-2026-11-15/additional_rnaseq.json) · [Melbourne study](https://doi.org/10.1038/s41598-018-29462-y) · [Mumbai study](https://doi.org/10.1093/humrep/deaa255) · [GSE315857 primary study](https://doi.org/10.1093/reprod/xaag038)

## Reproduction

The [input manifest and commands](../../evidence/genetic-anchors-2026-11-15/README.md) use the existing verified cohort readers. Separate enumeration via raw per-stratum group means reproduced the six primary extreme counts in Rome, shared scan days, and Mumbai. Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and openpyxl 3.1.5. Enumeration uses no seed. Public aggregate outputs contain no sample-level rows.
