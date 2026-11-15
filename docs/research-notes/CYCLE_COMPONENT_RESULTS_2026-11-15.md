# How much of the Rome reversal follows the Melbourne cycle patterns?

**Exploratory secondary analysis, November 15, 2026.** On 10,809 identically mapped genes, the Finnish score was 4.525 points lower in Rome cases. The component aligned with the seven Melbourne histologic stages contributed −0.543 points. The remaining component contributed −3.981 points. The observed reversal persists after mathematically removing the six measured stage contrasts. The result narrows an explanation for the cross-cohort disagreement while leaving the acquisition-order problem unresolved. [Gene and source checks](../../evidence/genome-transfer-2026-11-15/cycle-components.json) · [Diagnostic specification](CYCLE_COMPONENT_PLAN_2026-11-15.md)

## Two different quantities

In 257 Melbourne biopsies, seven histologic stages explained 71.5% of variation in the full Finnish score after accounting for sequencing batch. That percentage depends on how genes vary together across Melbourne biopsies. To examine what happens in other studies, we used the genes mapped in all five settings (Finnish training, Melbourne, Rome, Mumbai and Swiss glandular cells), retaining 90.7% of the original absolute Finnish weight. We estimated six independent batch-adjusted vectors of between-stage expression differences in Melbourne. We split the Finnish weights into the portion along these vectors and the remainder perpendicular to them. The two parts use the same normalization so their effects add to the original common-gene contrast in each study. No parameters were fitted to Rome, Mumbai or Swiss outcomes. The diagnostic was chosen after seeing the earlier contrasts. [Melbourne phase analysis](MELBOURNE_PHASE_RESULTS_2026-11-15.md) · [Implementation](../../scripts/cycle_component_transfer.py)

| Contrast on the common genes | Full score | Along the six stage contrasts | Remaining score |
| --- | ---: | ---: | ---: |
| Rome, cases minus controls on two mixed scan dates | −4.525 | −0.543 | −3.981 |
| Mumbai, cases minus controls in mid-proliferative tissue | −0.029 | −0.159 | +0.130 |
| Swiss glandular cells, cases minus controls within capture batch | +0.034 | +0.100 | −0.066 |
| Melbourne, mid-secretory minus mid-proliferative, fitted stage reference | +2.080 | +2.080 | 0.000 |

The Melbourne row is an arithmetic identity from fitting these stages, using all seven stages and the reduced common feature set. Its value differs from the primary +1.939-point contrast, which used just the two specified stages and 11,995 genes. The data used to define the stage patterns cannot validate a phase-independent biomarker. In Rome, the additive split describes how expression directions line up in gene space. A larger reversal in the remainder does not estimate a percentage caused by the scanner, disease or any other exposure. Restricted sampling, disease-dependent cycle effects, histology errors, cell mixture, processing and assay differences can all leave structure outside the fitted patterns. In particular, every control preceded every case on both mixed Rome dates. [Rome native scan times and technical control results](GENOME_TRANSFER_RESULTS_2026-11-15.md)

## Reproduction

Use the same local inputs as the [Rome and Mumbai analysis](GENOME_TRANSFER_RESULTS_2026-11-15.md), the [Melbourne analysis](MELBOURNE_PHASE_RESULTS_2026-11-15.md) and the [Swiss study](SWISS_GENOME_TRANSFER_2026-11-15.md). The implementation checks source digests and the Swiss sentinel. Changing the reference stage and using QR in place of a singular value decomposition gave the same projected weights within 10⁻⁹. Computing external differences by donor and by gene agreed within 10⁻¹⁰ points. The residual was orthogonal to every fitted Melbourne stage contrast under the same tolerance. The signed components reconstruct the original common-gene contrast in each cohort. Results do not rely on sampling a random reference distribution.

```sh
d=data/private/gse266265
OPENBLAS_NUM_THREADS=1 python scripts/cycle_component_transfer.py data/private \
  "$d/GSE266265_epithelial_cells_monocle3.rds" "$d/41467_2026_73020_MOESM3_ESM.xlsx"
```
