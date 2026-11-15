# Five-gene lesion signature in matched spatial compartments

**Exploratory human-data computation, November 15, 2026.** A score learned from ovarian lesions shifted in the prespecified direction in all five stromal pairs and four of five epithelial pairs from superficial peritoneal lesions. No compartment crossed the three-test Holm threshold. Three genes rarely exceeded the spatial assay’s limit of quantification. Two genes, **GABRP** and **SALL1**, had measurable, lower epithelial expression in the authors’ filtered analysis. This experiment narrows the question to a possible epithelial contribution, with insufficient evidence to localize the entire five-gene signal. [Full results](../../evidence/spatial-lesion-panel/results.json) · [Frozen ovarian gene panel](../../evidence/lesion-transfer/panel.json) · [Local testing specification](SPATIAL_PANEL_PLAN_2026-11-15.md)

## Five independent anatomical pairs

Burns et al. used archived secretory-phase tissue from five women undergoing lesion excision and hysterectomy in a study approved by Ponce Health Sciences University. Each had a matched superficial peritoneal lesion and eutopic uterine endometrium. We checked that all five deposited patient keys link both tissue types, each with two regions per compartment, for 60 region measurements. Pan-cytokeratin marked epithelial regions, CD68 marked macrophage regions, and the stromal regions were negative for pan-cytokeratin, ACTA2 and CD68. Segments contain mixtures of neighboring cells. The source reports patient-matched acquisition and an independent clinical setting from the original ovarian discovery study. The study has no disease-free or normal peritoneum arm. [Primary study and methods](https://doi.org/10.1016/j.isci.2025.111790) · [GSE263897](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263897)

For each of the 60 raw digital count files we mapped all 18,676 single-probe endogenous genes through the deposited probe configuration. We ranked all genes within each region, including zeros tied at their midpoint, and calculated each panel gene’s position on the 0–100 percentile scale. The five genes were equally weighted, with positive signs for TMEM255A, DIRAS3 and JCHAIN and negative signs for GABRP and SALL1. We averaged region ranks within each patient, tissue and compartment, took paired lesion-minus-uterus differences, and averaged over patients. Sign flips operate on complete patient differences. At five pairs, even five positive signs only achieve 1/32 one-sided; the first of three Holm thresholds is 0.0167. Neither a second region nor an additional gene adds an independent patient. [Computation](../../scripts/spatial_lesion_panel.py)

| Compartment | Positive / paired patients | Score change, pp | Exact one-sided tail | Holm p | Detection-filtered regions / pairs | Filtered score, pp |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Stroma | 5 / 5 | +11.72 | 1 / 32 | 0.094 | 17 / 4 | +7.85 |
| Epithelium | 4 / 5 | +11.44 | 2 / 32 | 0.125 | 18 / 5 | +9.47 |
| Macrophages | 3 / 5 | +1.08 | 10 / 32 | 0.313 | 14 / 3 | +2.32 |

For the sensitivity analysis, the local detection rule retained regions with at least 10% of endogenous counts above the threshold estimated from 139 negative probes. It retained 49 regions. The source used outlier-probe and tissue/sequence QC, reported retaining 50 and then applied Q3 normalization. The local rule does not reproduce the complete published pipeline. Filtered immune and stromal analyses lost matched patients. The filtered Holm p-values were 0.188, 0.188 and 0.250 in table order. All raw region rank scores retain low-count signals; the following audit is essential to interpreting them.

## Which genes were quantifiable?

The authors’ three supplementary differential-expression tables contain the same 7,945 filtered biological genes plus one negative control. Of the five candidate genes, only GABRP and SALL1 passed through to these tables. We confirmed the lesion-minus-uterus orientation using complement C3, explicitly described as increased in lesion epithelium. Their source model used log2 Q3-normalized counts and patient as a covariate, followed by limma empirical Bayes inference. The raw-region reanalysis above first averaged replicate regions and used patient-paired exact tests.

| Gene | Stromal log2FC | Epithelial log2FC | Macrophage log2FC | Epithelial nominal p / original BH | Epithelial Holm across 15 planned tests |
| --- | ---: | ---: | ---: | ---: | ---: |
| GABRP | −0.702 | −1.082 | +0.068 | 0.0176 / 0.268 | 0.247 |
| SALL1 | −0.178 | −0.853 | −0.320 | 0.0158 / 0.254 | 0.237 |

GABRP exceeded the local limit of quantification in 16 of 20 epithelial regions but in no stromal region. SALL1 exceeded it in four of 20 epithelial regions and three of 20 stromal regions. Across all 60 regions, TMEM255A exceeded the threshold once, DIRAS3 twice, and JCHAIN four times. All five target IDs map uniquely in the physical probe panel. For the three absent from the filtered source tables, the author-model effects are unmeasured. The 15-hypothesis adjustment retains those nine tests with p = 1 for multiplicity bookkeeping. A rank shift below the limit of quantification can reflect counting background. The five positive stromal score differences therefore offer weak evidence about expression intrinsic to stromal cells. Control for cell composition within each segment, lesion environment, sample quality and a larger independent patient series would be needed to identify the source of the tissue difference. [Spatial study supplementary tables S1–S3](https://pmc.ncbi.nlm.nih.gov/articles/PMC11810701/#supplementary-material) · [Source data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263897)

## Reproduction

The SHA-256-verified inputs are the locked five-gene Entrez list, GEO family SOFT, native DCC archive, matching WTA PKC probe configuration and published workbooks S1–S3. The offline analysis uses Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and openpyxl 3.1.5. The source enumeration, all six SHA-256 digests, exact URLs, bounded opt-in acquisition and transformation rules are in the [analysis script](../../scripts/spatial_lesion_panel.py). Sparse missing DCC probe entries represent zero counts, consistent with the manufacturer’s [reader implementation](https://github.com/Nanostring-Biostats/GeomxTools/blob/master/R/readNanoStringGeoMxSet.R). Only donor-aggregated outputs are published.

```bash
export PYTHONPATH=src:scripts OPENBLAS_NUM_THREADS=1
mkdir -p results/spatial-lesion-panel
python scripts/spatial_lesion_panel.py --fetch > results/spatial-lesion-panel/results.json
cmp results/spatial-lesion-panel/results.json evidence/spatial-lesion-panel/results.json
python -m pytest -q tests/test_spatial_lesion_panel.py
```

Omit `--fetch` for a cached, offline run. The complete panel was fixed before evaluating its compartment-specific results. Local diagnostic plans lack an external preregistration timestamp.
