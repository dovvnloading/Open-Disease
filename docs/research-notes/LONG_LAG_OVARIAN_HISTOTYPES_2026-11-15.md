# Ovarian cancer histotypes a decade after endometriosis diagnosis

**Exploratory computation from published aggregate human data, 15 November 2026.** Clear-cell and endometrioid cancers together made up an unusually large share of the ovarian cancers diagnosed at least a decade after an endometriosis diagnosis in three national registries. The other histotypes still accounted for most cases in each cohort. This comparison measures the distribution of cancers among diagnosed cases. Absolute cancer incidence in an individual requires separate information.

| Study and lag from diagnosis | Clear-cell or endometrioid cancers / all ovarian cancers | Comparison within that registry | Estimate (95% interval) |
| --- | ---: | --- | ---: |
| Finland, ≥10 years | 35/86 (40.7%) | Pair represented 19.8% of expected cancers after age and period standardization | Ratio of SIRs 2.78 (1.76–4.36) |
| Sweden, ≥10 calendar years | 20/68 (29.4%) | Pair represented 13.4% of expected cancers after age, period, education and region standardization | Ratio of SIRs 2.70 (1.52–4.63) |
| Netherlands, 10–24 years | 30/69 (43.5%) | 41/241 (17.0%) after a benign dermal nevus, stratified into three 5-year follow-up windows | Conditional case-mix odds ratio 3.73 (1.99–7.02) |

The excess can be counted in the Nordic cohorts. Finland observed 35 cancers in the two histotypes, against an estimated 8.96 expected, and 51 in the remaining histotypes, against 36.31 expected. Sweden observed 20 against 4.79 expected, and 48 against 31.00 expected. The two histotypes accounted for approximately 64% and 47%, respectively, of each study’s total observed-minus-expected case count at long lag. These descriptive shares use registry reference rates and carry sampling uncertainty. They are distinct from the fraction of cancers caused by endometriosis.

**Overlap with existing research.** PR #46 also analyzed the Finnish, Swedish and Dutch registry sources. Those calculations are parallel examinations of overlapping participant populations, not three additional independent cohorts. [Earlier registry comparison](OVARIAN_HISTOTYPES_LONG_LAG_2026-11-15.md).

## What this adds

Previous studies established preferential associations with clear-cell and endometrioid ovarian cancer. Here the question was whether the distribution remains shifted after the diagnostic episode is more than a decade in the past, against the full remaining ovarian cancer denominator, and whether the effect survives in a different comparator design. We chose the histotypes and long-lag window during a targeted literature review. The comparisons are retrospective. An independent, prespecified replication would require a fourth population with accessible long-lag data. [Danish primary study](https://pubmed.ncbi.nlm.nih.gov/27430396/)

A same-operation explanation alone cannot account for cancers registered ten or more years after the exposure record. Longitudinal differences in surgery, age, ascertainment and classification can still change which cancers appear. There is considerable geographic and design variation in the risk estimate for any ovarian cancer. For example, the Dutch investigators reported an age-adjusted incidence rate ratio of 1.01 (0.75–1.35) after excluding the first year in their endometriosis cohort. Even with that near-null overall rate ratio, the distribution of long-lag cancers differed in our analysis. [Primary study](https://pubmed.ncbi.nlm.nih.gov/34266690/) · [Full tables, chapter 3](https://cris.maastrichtuniversity.nl/ws/files/93704393/c7507.pdf)

## Calculation

For Finland and Sweden, let K denote the sum of observed cases in the two histotypes and N all observed ovarian cancers. For each histotype reconstruct its expected count as observed divided by the published standardized incidence ratio (SIR). Sum the two expected counts to obtain E_pair. Obtain E_total using N divided by the SIR for all ovarian cancers and E_other by subtraction. Divide K/E_pair by (N−K)/E_other. The null value of 1 corresponds to the same standardized incidence ratio for the two groups. With fixed expected counts and independent Poisson counts for mutually exclusive cancer types, K conditional on N is binomial with null share E_pair/E_total. We transformed an exact two-sided Clopper–Pearson interval on that share to an interval on the ratio of SIRs. Interval estimates condition on the registry reference rates and assume additive counts across mutually exclusive histotypes. The two-sided binomial p-values are exploratory and unadjusted for analysis choices.

For the Netherlands, compare the odds of a cancer falling in the two histotypes among people with histological endometriosis and people with a benign nevus. Condition on both margins of each 2×2 table and estimate a common odds ratio across 10–14, 15–19 and 20–24 years. Convolving the Fisher noncentral hypergeometric distributions gives the conditional maximum likelihood estimate, an equal-tailed exact 95% interval and a central two-sided p-value. These data permit stratification by lag only. A separate unmeasured age effect within lag windows is possible. The calculations use the observed-number columns in supplementary Table S3. Its displayed incidence rate ratio for the zero exposed clear-cell events at 20–24 years appears inconsistent. The observed counts for clear-cell exposed cases sum to the 38 in the main table across all five time bands. The observed control count for clear-cell cancer across the five bands is 34, matching Table 4 but differing from 32 in Table 3 of the author's thesis. We used the time-since-diagnosis breakdown from Table S3 for the long-lag analysis. [Author's tables, chapter 3](https://cris.maastrichtuniversity.nl/ws/files/93704393/c7507.pdf)

## Checks against alternative interpretations

- **Rounding.** Allowing each long-lag SIR to range over half a unit of its displayed precision gave point estimates of 2.77–2.80 in Finland and 2.57–2.83 in Sweden. The most conservative lower confidence bounds across those combinations were 1.75 and 1.45.
- **Unspecified histology.** As a deliberately adverse Dutch sensitivity, we assigned all 25 control diagnoses recorded as adenocarcinoma of unspecified type to the histotype pair and retained the 7 unspecified exposed diagnoses in the remainder. The conditional odds ratio fell to 2.00 (1.10–3.63). This checks one classification error in one direction. Other definitively coded cancers could also have been misclassified.
- **Adenomyosis.** In the Dutch adenomyosis arm, the histotype pair comprised 33 of 136 cancers at 10–24 years. Against the same nevus controls the lag-stratified odds ratio was 1.63 (0.93–2.84). The interval includes 1. Concurrent endometriosis and adenomyosis were possible, and the two exposure arms share controls. These estimates cannot establish a difference between the two conditions.
- **Earlier intervals and individual histotypes.** The code reports every earlier Nordic time band available in the source tables. Finland’s first six months included five cancers and the very wide interval spans 1. The two histotypes are not interchangeable. In the long-lag stratum the published clear-cell SIR was 5.50 (3.08–9.06) in Finland and 6.2 (3.3–10.6) in Sweden. For endometrioid disease the corresponding SIRs were 3.21 (1.96–4.95) and 2.6 (1.0–5.3).
- **Absolute risk and scope.** The Finnish investigators estimated about two additional ovarian cancers per 1,000 women with surgically diagnosed ovarian endometriosis over ten years of follow-up. That estimate concerns a specific subtype of endometriosis and a different time horizon. This analysis cannot establish screening benefit, malignant transformation of individual lesions or treatment effects. [Finnish paper](https://pubmed.ncbi.nlm.nih.gov/29742675/)

## Cohorts and limitations

| Primary source | Exposure, comparator and relevant selection |
| --- | --- |
| [Saavalainen et al., 2018, Table 5](https://www.endometriozisdernegi.org/konu/dosyalar/pdf/makale_ozetleri/Haziran2018/makale21.pdf) | Finnish surgical endometriosis compared with national rates. Adenomyosis alone excluded. Ovarian analysis comprised 36,524 women from the 49,933-person cohort, excluded index oophorectomy and ended follow-up at either unilateral or bilateral oophorectomy. General-population reference rates lack equivalent surgery adjustment. |
| [Moberg et al., 2026, Table 2](https://obgyn.onlinelibrary.wiley.com/doi/full/10.1111/aogs.70223) | Swedish endometriosis or adenomyosis records compared with women without either recorded diagnosis during the study. The 49,133 exposed women include prior or undifferentiated disease. Most diagnoses used the broad N80 code. No comparable exclusion for prior organ removal. Calendar-year differences can differ from exact elapsed time. |
| [Hermens et al., 2021](https://pubmed.ncbi.nlm.nih.gov/34266690/) and [Table S3 of chapter 3 in the author's thesis](https://cris.maastrichtuniversity.nl/ws/files/93704393/c7507.pdf) | Dutch pathology cohort with 51,544 histologically recorded endometriosis diagnoses and a 132,654-person benign dermal nevus control group in the source study. Women could also have adenomyosis. Ovarian outcomes included tubal and peritoneal primaries. We used the available 10–24-year case counts and treated only this cohort as the Dutch comparison. Related ENOCA publications share participants. |

The three populations provide separate observations with different methods. We retained separate estimates because the mix of histotypes, exposure definitions, surgery and standardization differs. Study selection was targeted and the interval estimates do not account for literature or design selection. Individuals without a recorded endometriosis diagnosis may still have had the disease. A decisive follow-up would test a locked long-lag contrast in a separately sourced registry while accounting for age, diagnostic period, ovarian surgery, adenomyosis and central histotype classification.

## Reproduce

The two small CSV files contain selected numeric facts transcribed from the cited tables. No participant data or source PDFs are included. The standard input checksums below refer to these transcriptions. The study uses deterministic calculations in Python 3.12.14, NumPy 2.3.5 and SciPy 1.17.0. Run from the repository root after installing the project dependencies. Minor floating-point differences are possible with other SciPy releases. The test verifies the primary estimates and compares the exact-distribution implementation against independent combinatorial weights. [Script](../../scripts/ovarian_histotypes.py) · [All estimates and sensitivity analyses](../../evidence/long-lag-histotypes/results.json)

```sh
mkdir -p results/long-lag-histotypes
python scripts/ovarian_histotypes.py > results/long-lag-histotypes/reproduced.json
diff -u evidence/long-lag-histotypes/results.json results/long-lag-histotypes/reproduced.json
python -m pytest -q tests/test_ovarian_histotypes.py
```


| Input | SHA-256 |
| --- | --- |
| [dutch_counts.csv](../../evidence/long-lag-histotypes/dutch_counts.csv) | `3f1325fe88e3082c53a08f64907c3419a046393f57de6127d7822409db086e37` |
| [nordic_counts.csv](../../evidence/long-lag-histotypes/nordic_counts.csv) | `b9d8f9c355592db5a7968752d727137fef537bf8cdcf97ce9cbf4e21e57996c3` |
