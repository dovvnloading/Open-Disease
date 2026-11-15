# The fixed five-gene score in 11 paired CodeLink biopsies

The signed score rose in every lesion relative to the same woman's uterine biopsy. Its mean increase was **24.91 within-specimen expression percentile points**. Exactly one of 2,048 within-pair sign assignments had a mean as high as the observed one (one-sided p = 0.000488). The two-sided sign-only p was 0.000977. Each patient contributed one difference. Removing one pair at a time gave means from +23.67 to +26.66 points. This is exploratory secondary analysis of real human biopsies. [Frozen gene list](../../evidence/lesion-transfer/panel.json) · [Source study](https://doi.org/10.1016/j.fertnstert.2007.01.056)

The five genes and their directions came from a separate study of 10 ovarian lesion–uterus pairs. The selection did not use CodeLink or Hull outcomes. In the nine-pair Hull study, the signed score averaged +22.18 points and seven pairs were positive (one-sided p = 0.0371). The different research teams, platforms and recruitment sites support study-level separation. All three studies used tissue from patients with endometriosis. Healthy tissue from the same anatomical site would answer a separate question about disease specificity. [Discovery and first external cohort](LESION_TRANSFER_RESULTS_2026-11-15.md) · [CodeLink analysis plan](CODELINK_PANEL_PLAN_2026-11-15.md)

| Gene | Frozen lesion direction | Mean lesion minus uterine rank, points | Pairs in frozen direction | Directional p | Holm across five |
| --- | --- | ---: | ---: | ---: | ---: |
| TMEM255A | Higher | +8.73 | 8/11 | 0.1880 | 0.1880 |
| DIRAS3 | Higher | +33.37 | 11/11 | 0.000488 | 0.00244 |
| GABRP | Lower | −38.53 | 11/11 | 0.000488 | 0.00244 |
| SALL1 | Lower | −4.17 | 9/11 | 0.01563 | 0.03125 |
| JCHAIN | Higher | +39.73 | 11/11 | 0.000488 | 0.00244 |

The matrix contains 54,359 distinct finite probe rows from GPL2895, which has 56,448 annotated rows. Each complete deposited feature supplied the within-specimen reference. Only uniquely mapped probes supplied candidate values. SALL1 used three eligible probes and the other genes each used one. We multiplied the five gene differences by their frozen signs and averaged them with equal weights. The p-values assume tissue labels can be exchanged within each patient under the sharp null. Different cell composition, tissue of origin or a tissue-specific processing effect could each make a tissue contrast reproducible. This design cannot estimate diagnostic performance or establish gene function. [Aggregate results](../../evidence/lesion-transfer/codelink.json) · [Implementation](../../scripts/lesion_codelink.py)

GEO's curated individual subsets place one uterine biopsy and one lesion in each pair. All submitted numeric donor identifiers agree with that pairing. One pair differs in a trailing letter in the original titles. Requiring an exact suffix match leaves ten pairs, all positive (mean +26.66). Six explicitly labeled ovarian endometriomas averaged +29.21 points, three explicitly labeled peritoneal lesions +18.70, and two without a specific lesion site +21.31. The small subgroups serve as descriptive checks. [GEO curated design](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3092) · [Published series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5108)

An independent reader matched the 22 array identifiers across the differently ordered source tables, recovered the 11 curated individual and two tissue subsets, and enumerated the tests with integer subset sums. All five gene estimates, the panel, subgroup summaries, sensitivities and adjusted p-values agreed. [Independent verification](../../evidence/lesion-transfer/codelink-verification.json)

## Reproduce

Keep the following original files in the private directory. All are checksum-verified by the offline script. The raw samples and donor pair assignments remain outside Git.

| File | Public source |
| --- | --- |
| GSE5108_series_matrix.txt.gz | [GEO](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE5nnn/GSE5108/matrix/GSE5108_series_matrix.txt.gz) |
| GPL2895.annot.gz | [GEO platform](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL2nnn/GPL2895/annot/GPL2895.annot.gz) |
| GDS3092.soft.gz | [GEO curated dataset](https://ftp.ncbi.nlm.nih.gov/geo/datasets/GDS3nnn/GDS3092/soft/GDS3092.soft.gz) |

```sh
mkdir -p results/lesion-transfer
OPENBLAS_NUM_THREADS=1 python scripts/lesion_codelink.py data/private/genetics-expression > results/lesion-transfer/codelink.json
cmp results/lesion-transfer/codelink.json evidence/lesion-transfer/codelink.json
```
