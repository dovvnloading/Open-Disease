# Finnish genome-wide expression profile in Swiss glandular cells

**Exploratory human-data computation, November 15, 2026.** In 6 cases and 7 controls, the full Finnish signed expression score differed by **+0.023 percentile points** between groups after matching capture batch. Among 144 assignments preserving case counts in each batch, 127 were at least as extreme (two-sided p = 0.882). The exact 95% interval under a constant additive shift was −0.204 to +0.413 points. It is compatible with both directions. The Swiss epithelial comparison addresses a different tissue compartment from the whole-biopsy analyses. [Swiss source study](https://doi.org/10.1038/s41467-026-73020-4) · [Aggregate output](../../evidence/genome-transfer-2026-11-15/swiss.json)

We reused the published Swiss donor selection and the bulk study's continuous Finnish weights. The clinical workbook, 199,996,542-byte epithelial trajectory RDS, ordered feature axis, dated HGNC mapping, and Finnish weights passed their pinned SHA-256 checks. Each donor contributed at least 20 qualifying proliferative glandular cells. The 3,564 selected cells supplied one pseudobulk profile per donor. We ranked all 29,034 original features, including zero ties, after summing SCT-corrected counts. Exactly 10,945 of 12,044 Finnish genes mapped uniquely, retaining 91.73% of absolute Finnish weight. The score is the sum of percentile ranks times signed weights, divided by the sum of absolute measured weights. [Swiss eligibility and annotation](SWISS_GLANDULAR_PLAN_2026-11-14.md) · [Source receipts](GSE266265_EPITHELIAL_FEASIBILITY_2026-11-14.md)

Each of four mixed capture batches contributed the same weight to the case-minus-control contrast. Omitting one batch at a time gave differences ranging from −0.056 to +0.099 points. Removing all measured genes in the earlier nominations left 10,830 genes and an effect of +0.020 points (129/144 assignments at least as extreme). The same ranks reproduced the earlier five-pathway average at +0.709243 points and its exact tail of 24/144. Separate matrix and rational-arithmetic implementations agreed on the full-score tail count. Small donor numbers and SCT correction limit biological interpretation. This comparison cannot identify which cell types account for differences between whole-biopsy studies. [Earlier Swiss results](SWISS_GLANDULAR_RESULTS_2026-11-14.md)

The mapping, whole-profile score, mixed-batch test, previous-nomination exclusion and capture-batch omission were written down locally before calculating Swiss outcomes. The earlier pathway result was already known. The compatibility interval was added after the new score was computed. The p-value is descriptive and was outside the Rome–Mumbai two-cohort Holm family.

## Reproduce

Use base R, the project Python dependencies and openpyxl. Place the two linked human source files in a private working directory and obtain the pinned HGNC archive. The runner uses local input files, writes only aggregate JSON to stdout, and deletes transient cell-level data on exit. The verified run used R 4.5.0, Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. [Analysis code](../../scripts/swiss_genome_signature.py) · [HGNC archive](https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2026-10-06.tsv)

```sh
d=data/private/gse266265
OPENBLAS_NUM_THREADS=1 python scripts/swiss_genome_signature.py \
  "$d/GSE266265_epithelial_cells_monocle3.rds" \
  "$d/41467_2026_73020_MOESM3_ESM.xlsx" \
  "$d/hgnc_complete_set_2026-10-06.tsv"
```
