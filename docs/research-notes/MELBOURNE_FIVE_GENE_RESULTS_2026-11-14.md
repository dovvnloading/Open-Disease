# Melbourne follow-up of the five locked genes

**November 14, 2026 · Computational analysis of published human cohort summaries.** One of the five nominated genes had an eligible probe in the complete high-detection table. Its two-sided test for GLTP gave p = **0.2116723**, with Holm-adjusted p = **1.0** across the original five targets. The other four genes were unmeasured under the registered rule. Melbourne supplies no additional statistically supported gene association in this test. [Public pre-analysis commit](https://github.com/dovvnloading/Open-Disease/blob/95cae294edd57367735d2033d34e71dd374ed4ce/docs/research-notes/MELBOURNE_FIVE_GENE_PLAN_2026-11-14.md) · [Original five-gene nomination](https://github.com/dovvnloading/Open-Disease/commit/f2a3755ff055cde241e0eaaaedc11a1471270beb) · [Issue #37](https://github.com/dovvnloading/Open-Disease/issues/37)

| Locked gene | Finnish direction | Eligible S20 probes | Two-sided probe p | Within-gene adjusted p | Five-target Holm p |
| --- | --- | ---: | ---: | ---: | ---: |
| ABO | Higher | 0 | — | — | — |
| KIF3C | Lower | 0 | — | — | — |
| ZNF835 | Lower | 0 | — | — | — |
| GLTP | Higher | 1 | 0.2116723 | 0.2116723 | 1.0 |
| F11 | Lower | 0 | — | — | — |

We checked all 15,262 distinct probe rows in Table S20, used exact trimmed ILMN_GENE symbols and retained every matching probe. GLTP probe ILMN_1764380 has a published `logFC` of −0.04005445. The table title and paper identify a comparison between cases and controls but do not explicitly identify the coefficient's reference level. The sign of disease-associated change is unassessable from this evidence. Each uncovered target contributed p = 1 only to the five-test Holm calculation. Its biological effect remains unknown. We did not evaluate the differently defined expression-detection endpoint in Table S21. [All probe results](../../evidence/cross-cohort-2026-11-14/melbourne_probes.tsv) · [All gene results and multiplicity bookkeeping](../../evidence/cross-cohort-2026-11-14/melbourne_genes.tsv)

## What this comparison contains

Fung et al. measured uterine endometrial expression in 229 women from two Melbourne clinics. Of 165 surgical patients, 112 had endometriosis at laparoscopy, 48 had a negative recorded diagnosis and five had unknown status. In the IVF group the 64 participants provided 21 positive, 41 negative and two unknown self-reports. The combined enrollment counts are 133 positive, 89 negative and seven unknown. Participants in the surgical group had a three-month hormone washout. Some IVF participants underwent ovarian stimulation in the preceding cycle. Fung et al. adjusted the disease comparison for menstrual cycle stage. The deposited summary does not give the model sample size or the treatment of the seven unknown diagnoses. The recovered table contains probe-level statistics already computed by the authors. It has no donor-level expression values or covariances for re-estimating the Finland–Rome pathway scores. [Primary publication, Methods, Results and Table 7](https://doi.org/10.1038/s41598-018-29462-y) · [Supplementary release audit](GSE266265_SUMMARY_FEASIBILITY_2026-11-14.md)

## Reproduction and checks

The source was a 4,476,677-byte CC BY 4.0 supplementary workbook. SHA-256 was `b83a3d32a5882f750d978f4de1a6473fdd5d7627dbccf7d7c1774162e104e8e3`. The high-detection worksheet TableS20 has a title, column labels and 15,262 distinct Illumina probe identifiers. The paper reports the same probe count for 12,321 high-detection genes. [Fung et al. supplementary workbook](https://pmc-oa-opendata.s3.amazonaws.com/PMC6065421.1/41598_2018_29462_MOESM2_ESM.xlsx)

```bash
p=evidence/cross-cohort-2026-11-14
book=data/private/fung-2018/Fung_supplementary_tables.xlsx
python scripts/analyze_melbourne_genes.py --workbook "$book" --nominations "$p/nominations.json" --out "$p"
python scripts/verify_melbourne_s20.py "$book" "$p/nominations.json"
python -m pytest -q tests/test_analyze_melbourne_genes.py
```

The primary analysis used Python 3.12.14, openpyxl 3.1.5 and NumPy 2.3.5. The cached-input run took 2.01 seconds and used 100,584 KiB peak RSS. A separate standard-library implementation decoded the workbook's OOXML relations and checked all S20 rows using decimal arithmetic. It agreed on the single probe, its coefficient, nominal and adjusted probabilities and absence of the other four nominations. A separate base R Holm calculation also returned 1.0 for the measured gene. Nine synthetic tests passed, including ambiguous inputs, within-gene Bonferroni, missing-target handling and an end-to-end streamed workbook. No random procedure applies. [Primary implementation](../../scripts/analyze_melbourne_genes.py) · [Independent verification](../../scripts/verify_melbourne_s20.py) · [Regression tests](../../tests/test_analyze_melbourne_genes.py)
