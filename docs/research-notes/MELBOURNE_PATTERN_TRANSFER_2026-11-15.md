# Does the Melbourne expression pattern recur in Finnish or Rome uterine biopsies?

**15 November 2026 · Exploratory human-data computation.** The Melbourne disease coefficients showed little alignment with either case–control comparison across 6,908 shared genes. Both exact two-sided tests were inconclusive. An independently published early-secretory versus proliferative pattern showed much stronger alignment with Finland on its 872 measured phase-responsive genes. A strong gene-wise correlation was also plausible for chance relabelings of these small donor groups.

| External Melbourne pattern | Shared genes | Finland, 12 cases and 13 controls | Rome, 8 cases and 6 controls | Rome restricted to two mixed scan dates |
| --- | ---: | --- | --- | --- |
| Cycle-adjusted disease coefficient, Table S20 | 6,908 | r = −0.064, 5,083/12,096, p = 0.420 | r = +0.048, 1,534/3,003, p = 0.511 | r = +0.063, 43/200, p = 0.215 |
| Early secretory minus proliferative, Table S3 | 872 | r = +0.777, 754/12,096, p = 0.062 | r = −0.370, 1,280/3,003, p = 0.426 | r = −0.218, 101/200, p = 0.505 |

Counts are assignments at least as extreme in either direction. Holm p-values for disease and phase, respectively, were 0.420 and 0.125 in Finland, 0.852 and 0.852 in full Rome, and 0.430 and 0.505 in scan-date-restricted Rome. The Melbourne paper and Table S20 do not fix the disease coefficient's reference level. Reversing that reference would reverse both S20 correlations and leave the tests unchanged. Table S3 explicitly labels up/down in early secretory phase. [Melbourne primary study](https://doi.org/10.1038/s41598-018-29462-y) · [Workbook](https://pmc-oa-opendata.s3.amazonaws.com/PMC6065421.1/41598_2018_29462_MOESM2_ESM.xlsx)

## Measurements and tests

Fung et al. recruited 229 women, including 165 surgical participants and 64 IVF participants. Of those, 133 had positive, 89 negative and seven unknown disease status across the enrollment groups. Some diagnoses were self-reported, cycle stage entered the disease comparison, and the exact disease model sample size is unavailable. S20 contains all 15,262 high-detection probes. S3 contains 1,555 externally selected differential probes and omits nonsignificant probes. We use S3 as a targeted phase pattern. It cannot estimate whole-transcriptome similarity to a complete cycle contrast. [Source audit and diagnosis counts](MELBOURNE_FIVE_GENE_RESULTS_2026-11-14.md)

Both comparisons start with the 12,044-Entrez Finnish screening universe. Trimmed Melbourne symbols must match exactly one GPL571 gene with exactly one valid symbol across unambiguously annotated non-control probes. We averaged the published moderated t over eligible probes per gene. S20 contributed 8,921 mapped probes across 6,908 genes. S3 contributed 1,047 probes across 872 genes. Historical annotation, restricting to genes mapped by both original arrays, and the Melbourne high-detection filter limit generalization. The case–control expression profiles are donor-level contrasts in within-array probe percentiles, divided gene by gene by the label-invariant null SD. Across genes we calculated Pearson correlation between the external reference and the case–control profile. Rome's unstratified 14-array SD is proportional to that contrast's null SD, and Pearson correlation is invariant to this positive constant. [Earlier Finland–Rome measurement and design](CROSS_COHORT_DISCOVERY_RESULTS_2026-11-14.md)

The exact Finnish analysis permutes complete donor label vectors within seven original run-by-histologic-phase strata and weights those strata equally. Rome uses all whole-donor allocations of eight cases across 14 proliferative biopsies. The post-hoc Rome day analysis holds the 14-array gene scales fixed, then calculates the mean within-day contrast on the two acquisition dates containing both groups, with five cases and six controls. It excludes three cases on other dates. Day matching leaves within-day acquisition order, treatment and indication unresolved. Gene-gene and probe-probe dependence are preserved in each allocation. All p-values are conditional on exchangeability within their analysis strata. Null distribution quantiles in the JSON summarize permissible allocations, not confidence intervals on population correlations.

## What the phase pattern can explain

Within Finnish proliferative histology the correlation with the Melbourne phase pattern was +0.628 (13/72, p = 0.181). Within secretory histology it was +0.674 (47/168, p = 0.280). Both phase-specific Holm p-values were 0.361. As a positive control, we used actual histologic phase as the exposure and kept only Finnish run-by-diagnosis strata containing both phases. Ten of the 25 donors qualified in four strata, five in each phase. The correlation was +0.507 (8/36, p = 0.222). These data offer limited resolution on individual cycle timing and tissue composition. The Finnish direction persisted within each broad phase.

On exactly the same 872 selected genes the Melbourne disease-coefficient correlations were −0.121 in Finland and +0.066 in Rome. The Finnish and Rome disease contrasts correlated at −0.477 on that panel. Projecting both contrast vectors off the external phase coefficient vector gave partial r = −0.324. This descriptive calculation removes one linear component across genes. Individual donors' timing and cell mixtures remain unknown, and their causal contribution cannot be quantified from it. The difference between these correlations cannot be interpreted as a fraction of original pathway reversals caused by phase.

## Reproducibility

Using native source data, we reconstructed every previously published Finnish statistic on the 6,908-gene panel exactly. A second implementation independently read the Melbourne workbook's OOXML, repeated symbol mapping, obtained the null SD directly from all possible whole-donor contrasts, and reproduced the extreme counts for both full-cohort tests in both cohorts, both phase subsets and the positive control. For scan-date conditioning, direct vector correlations agreed with the projection for all 200 assignments in each reference. Four synthetic checks cover correlated genes, orientation reversal, exhaustive stratified standardization, eligibility, partial correlations and constant genes.

The two initial Rome tests and the Finnish descriptive comparison were specified locally before this analysis. Exact Finnish testing, its positive control and phase splits, the same-gene sensitivity and Rome's mixed-date comparison followed the observed descriptive results. Earlier targeted analyses had already seen these cohorts. There was no separate blinded or public prospective registration. Two-sided tests and both references are retained. All comparisons use deterministic enumeration without a random seed.

Use the original native Finnish input recipe, the original Rome matrix and platform annotation, the scan-date header cache and the Fung CC BY 4.0 workbook. The checksum for Table S20's workbook is b83a3d32a5882f750d978f4de1a6473fdd5d7627dbccf7d7c1774162e104e8e3. The checksum for the Finnish table is 9df65b29615ec56767e9abaaf3844f8ebde5e21baf8cdab193ac24857a170da5. Matching SHA-256 checks on the native Finnish archive ranges, three platform manifests, Rome matrix and GPL571 annotation are mandatory in the existing repository loaders. Both native GEO and Figshare studies supply single-platform cohort definitions as previously audited. [Finnish archive, ranges and manifests](GSE141549_CHIP_CENSUS_2026-11-14.md) · [Clinical range instructions](GSE141549_ELIGIBILITY_2026-11-13.md) · [Rome input files](GSE25628_FOSB_RESULTS_2026-11-13.md)

The newly acquired research inputs for this analysis totaled 564,825,600 bytes, before the separate scan audit. The Finnish source loader took 75.1 seconds and 1,170,332 KiB peak RSS. The independent cross-check took 6.7 seconds once donor-level inputs were available. Runtime was Python 3.12.14, NumPy 2.3.5, openpyxl 3.1.5 and inflate64 1.0.4. The workbook and all array-level inputs remain local. The two result files contain only aggregates. [Melbourne transfer results](../../evidence/melbourne-transfer-2026-11-15/transfer.json) · [Finnish null tests, phase splits and controls](../../evidence/melbourne-transfer-2026-11-15/finland.json) · [Independent verification program](../../scripts/verify_melbourne_pattern_transfer.py)

From the repository root, pass the actual local input paths and an existing aggregate output directory.

```bash
export PYTHONPATH=src:scripts
export OPENBLAS_NUM_THREADS=1
data=data/private
book=$data/fung-2018/Fung_supplementary_tables.xlsx
table=evidence/cross-cohort-2026-11-14/discovery_gene.tsv
out=results/melbourne-transfer
mkdir -p "$out"
python scripts/melbourne_pattern_transfer.py --rome "$data/gse25628" --workbook "$book" \
  --finnish "$table" --scan-audit "$data/matrisome-audit" > "$out/transfer.json"
python scripts/finnish_cycle_alignment.py --rome "$data/gse25628" --workbook "$book" \
  --finnish "$table" --finland "$data/gse141549" > "$out/finland.json"
python scripts/verify_melbourne_pattern_transfer.py --rome "$data/gse25628" --workbook "$book" \
  --finnish "$table" --finland "$data/gse141549" --rome-results "$out/transfer.json" \
  --finnish-results "$out/finland.json"
python -m pytest -q tests/test_melbourne_pattern_transfer.py
```
