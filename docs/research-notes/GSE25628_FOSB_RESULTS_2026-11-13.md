# FOSB in a separate proliferative-phase endometriosis cohort

**Exploratory human-data computation, November 13, 2026.** FOSB expression was higher on average in the eight eutopic case biopsies than in the six surgical control biopsies. The estimated difference was **+1.382 log2 units** (95% CI **−0.647 to +3.411**, two-sided Welch **p = 0.156**). The prespecified decision is **inconclusive**. The interval admits both lower expression and a substantial increase. [Frozen protocol](https://github.com/dovvnloading/Open-Disease/commit/456455c5b27437fd001256f9edf65188da68377a) · [GEO GSE25628](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25628)

| Prespecified quantity | Result |
| --- | ---: |
| Women included according to the published sampling design | 8 cases, 6 controls |
| Paired lesion arrays excluded | 8 |
| Mean RMA log2 intensity, case / control | 6.540 / 5.158 |
| Case minus control, log2 units | +1.382 |
| Two-sided 95% Welch interval, log2 units | −0.647 to +3.411 |
| Ratio of geometric means (95% interval, back-transformed) | 2.61 (0.64 to 10.63) |
| Welch degrees of freedom | 8.362 |
| Two-sided p-value, one prespecified comparison | 0.1560 |
| Range after omitting each included biopsy in turn, log2 units | +0.667 to +1.651 |
| Probes present / tested | 22,277 / 1 |

## Interpretation

The Italian estimate has the same sign as the FOSB increase reported by Tamaresis et al. in the UCSF cohort. All 14 leave-one-out estimates were positive. Those correlated estimates describe sensitivity to individual biopsies and do not supply 14 independent confirmations. With six controls and eight cases, the interval is wide. [Original Table 4](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/#T4)

Crispi et al. recruited surgical cases in Rome. Each affected participant provided a uterine biopsy and a lesion biopsy, and each surgical control provided a uterine biopsy. This sampling description supports the independent observational units after excluding lesions. The GEO field called `donor` encodes clinical group. It supplies no individual participant key. The paper confirms proliferative phase for all array specimens. It describes the controls as women undergoing surgery for benign conditions and does not establish individual recent hormone exposure. Their indications, possible differences in treatment and the surgical case mix limit interpretation and comparability with the UCSF pathology-free contrast. Cross-study independence is assessed from the separate recruitment sources. This result supports no clinical use and leaves the original GSE51981 donor-level reproduction unresolved. [Crispi et al., Methods](https://doi.org/10.1002/jcp.24358) · [GSE51981 status](GSE51981_REPRODUCTION_STATUS_2026-11-13.md)

## Execution and source provenance

The expression matrix came directly from GEO. The 22 unique sample accessions and their column order matched the matrix header. The source characteristics split the arrays into eight eutopic cases, eight lesions and six controls. Both clinical groups provided complete intensities for the endpoint. All 22 source sample descriptions specify log2 RMA intensities, and GEO reports quantile normalization. The August 2016 GPL571 annotation mapped exactly one probe to the specified FOSB symbol and Entrez ID 2354. The analysis quantified only `202768_at`, without screening other probe outcomes. The normalizations differ between studies (RMA in GSE25628 and GC-RMA in GSE51981).

| Public input, retrieved November 14 UTC | Compressed or response bytes | SHA-256 |
| --- | ---: | --- |
| [GSE25628 Series Matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25628/matrix/GSE25628_series_matrix.txt.gz) | 2,795,955 | `3d77aae9769afb1f4bb1d03e2bc73fbab4ccc52ac5cd0f0945e4c0224fcca089` |
| [GPL571 annotation](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL571/annot/GPL571.annot.gz) | 4,504,226 | `9f5a482f9f3d478b47fc80124600319e569f1e9ab447f3f256d61b1fbacce4fa` |
| [GEO sample metadata, quick view](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25628&targ=gsm&view=quick&form=text) | 57,395 | `b57b522cb17114a9ed412e09b0d4fc532e92ae4dfd49fa3c05e15778ef7c8df6` |

These transfers total 7,357,576 bytes. Public access and local exploratory computation are consistent with the GEO disclaimer, subject to original submitters' rights. This report publishes only aggregates, original code and source references. [GEO terms](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)

## Reproduce

From a checkout containing this report, install the project as in the README. Use a fresh ignored input directory. The command checks both byte hashes before computation.

```bash
mkdir -p data/private/gse25628
cd data/private/gse25628
curl --fail --silent --show-error --max-filesize 4000000 --output GSE25628_series_matrix.txt.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE25nnn/GSE25628/matrix/GSE25628_series_matrix.txt.gz
curl --fail --silent --show-error --max-filesize 5000000 --output GPL571.annot.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL571/annot/GPL571.annot.gz
gzip -t GSE25628_series_matrix.txt.gz GPL571.annot.gz
cd ../../..
python scripts/analyze_gse25628_fosb.py \
  --matrix data/private/gse25628/GSE25628_series_matrix.txt.gz \
  --annotation data/private/gse25628/GPL571.annot.gz
```

The source processing record lists oneChannelGUI 1.16.2 and Bioconductor 2.7. The deterministic Python computation used CPython 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. An independent implementation parsed the complete matrix in base R 4.5.0 and reproduced the group means, contrast, Welch interval and p-value, geometric-mean ratio and all leave-one-out bounds. A timed repeat of the Python analysis took 3.69 seconds and used 81.8 MiB peak RSS. No random seed applies.
