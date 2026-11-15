# Can the Fonseca atlas test lesion specificity against peritoneum from participants without endometriosis?

**Literature and metadata assessment, 2026-11-15.** GSE213216 contains no peritoneal specimen from any of the four participants without confirmed endometriosis. Its four histology-negative peritoneal specimens all came from two people with confirmed disease. Consequently, this atlas leaves open whether the [five-gene panel](LESION_TRANSFER_RESULTS_2026-11-15.md) distinguishes a peritoneal lesion from peritoneum in a person without endometriosis. The [paired core and adjacent peritoneum result](ADJACENT_PERITONEUM_RESULTS_2026-11-15.md) addresses a different comparison.

| Candidate comparator | Published specimen and donor evidence | Consequence for the panel |
|---|---|---|
| Peritoneum from participants without confirmed endometriosis | Four participants provided three uterine endometrium specimens and two unaffected ovary specimens. None provided peritoneum. | The requested same-organ comparison has zero control donors. |
| Peritoneum with negative histology in participants with confirmed disease | Four biopsies from two donors. Both donors also had positive peritoneal biopsies. No positive and negative biopsy shared the same exact anatomic location label within either donor. Distance from a positive lesion was undocumented. | A secondary within-person comparison could describe differences across the sampled locations. Even if both donor differences favored the panel, the minimum one-sided exact sign-flip P would be 1/4. |

The negative biopsies appeared suspicious at surgery. The authors found no endometriosis on histology, including five deeper sections. They detected 28–248 endometrial-type epithelial cells per sequenced negative specimen and proposed that the pieces used for sequencing may have included microscopic lesions. That explanation remains unconfirmed. Following pooling and quality control the study reported 23 pathology-positive and four pathology-negative peritoneal libraries. Their large cell counts do not add independent donors. We closed this comparison at metadata review without downloading expression data. [Fonseca et al., *Nature Genetics* (2023), Results and Discussion](https://pmc.ncbi.nlm.nih.gov/articles/PMC10950360/) · [GSE213216](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213216).

## Reproduction

The audit reads explicit person-level diagnosis in S2, carries the merged donor key down across specimen rows in S1 and checks the S1 diagnosis against S2. It interprets the four peritoneal S1 pathology entries labeled NA† using the authors’ description of the four pathology-negative specimens. Exact location overlap uses the full anatomic label recorded in S1. The remaining tissue and diagnosis counts reproduce the primary article. The workbook contains 1,638,632 bytes and the script checks its SHA-256 before reading it. [Primary supplement, Tables S1–S25](https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41588-022-01254-1/MediaObjects/41588_2022_1254_MOESM4_ESM.xlsx) · [Aggregate output](../../evidence/lesion-transfer/fonseca-comparator-audit.json) · [Audit script](../../scripts/audit_gse213216_peritoneum.py).

Save the linked workbook under the ignored private input directory, install the repository’s dev dependencies, and run from the repository root.

```sh
python scripts/audit_gse213216_peritoneum.py data/private/gse213216/supp.xlsx
```
