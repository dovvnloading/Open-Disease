# Which ovarian cancer histotypes remain overrepresented a decade after an endometriosis diagnosis?

**Exploratory human-data computation · 2026-11-15**

Finland and Sweden each recorded clear-cell and endometrioid cancers more often than the expected mix of ovarian cancers at least ten years after the index diagnosis. Of 86 Finnish ovarian cancers, 35 (41%) had either histotype. Population rates implied a 20% share. In Sweden the figures were 20 of 68 (29%) and an expected 13%. A Dutch pathology cohort showed the same direction using an observed comparison group. The studies already reported associations for these histotypes. Pearce et al. had shown in seven case–control studies that associations persisted after excluding diagnoses within ten years, including an association for low-grade serous cancer. The contribution here is the numerical contrast between two histotypes combined and the other cancers in these registries. [Pearce et al., 2012, Table 4](https://discovery.ucl.ac.uk/id/eprint/1343061/1/1343061.pdf)

| Cohort and reference | Observed cases with either histotype / all ovarian cancers | Comparison | Contrast (95% CI) |
| --- | ---: | --- | ---: |
| Finland, surgically verified endometriosis, ≥10 years | 35/86 | 8.96 expected out of 45.26, using Finnish population rates | Ratio of standardized incidence ratios 2.78 (1.76–4.36) |
| Sweden, endometriosis or adenomyosis, ≥10 years | 20/68 | 4.79 expected out of 35.79, using rates among women without either recorded diagnosis | Ratio of standardized incidence ratios 2.70 (1.52–4.63) |
| Netherlands, histologically confirmed endometriosis, 10–24 years | 30/69 | 41/241 following a benign dermal nevus diagnosis | Case-mix odds ratio 3.73 (2.00–6.98) |

In Finland, the combined histotypes occurred 3.91 times as often as expected, and the remaining cancers occurred 1.40 times as often. The corresponding Swedish figures were 4.18 and 1.55. The reference distributions and exposure definitions differed, so we report country-specific estimates. The Dutch contrast conditions on already having developed cancer. None of these contrasts estimates a person's chance of developing ovarian cancer.

A related Finnish study checked whether clear-cell cancers also occurred in excess outside the ovary. It found 23 cases (about 26 expected, SIR 0.89, 95% CI 0.56–1.33). Nineteen were coded as kidney cancers, leaving little information on rare pelvic sites. The studies draw on the same surgical registry. Extraovarian follow-up began in 2007 when detailed morphology coding became available and used a ≥5-year stratum (SIR 0.92, 0.58–1.37). The ovarian study followed women from surgery and censored at oophorectomy. These eligibility and timing differences prevent a matched ≥10-year organ comparison. [Saavalainen et al., 2021, Table 2](https://doi.org/10.1371/journal.pone.0253270)

## Sources and dates

- **Finland.** Saavalainen et al. (2018), Table 5, all endometriosis, follow-up ≥10 years. The ovarian cancer cohort had 36,524 women. Those undergoing unilateral or bilateral oophorectomy at index surgery were excluded, and follow-up ended at oophorectomy. The reference rates account for age and calendar period. [Paper](https://doi.org/10.1097/AOG.0000000000002624) · [Full text and Table 5](https://www.endometriozisdernegi.org/konu/dosyalar/pdf/makale_ozetleri/Haziran2018/makale21.pdf#page=6)
- **Sweden.** Moberg et al. (2026), Table 2, ≥10 years of follow-up. The 49,133-person exposed cohort used the full N80 code range for endometriosis or adenomyosis. The comparator excluded both recorded diagnoses. Standardization included age, calendar period, education and region. Follow-up ended at first gynecological or breast cancer, death, emigration or study end. [Paper and Table 2](https://doi.org/10.1111/aogs.70223)
- **Netherlands.** Hermens et al. (2021), Supplementary Table S3, also printed in chapter 3 of the first author's thesis (printed page 53). The endometriosis cohort included 51,544 women. We used all published late bins (10–14, 15–19, 20–24 years), and used the nevus cohort's observed cancers as the reference. The study included cancers of the ovary, fallopian tube and peritoneum. It shares a registry with the authors' 2020 paper, which does not supply an additional population here. [Paper](https://doi.org/10.1016/j.ygyno.2021.07.006) · [Author's supplement](https://cris.maastrichtuniversity.nl/files/93704393/c7507.pdf)

We consulted public article tables on 2026-11-15. The input files contain selected numerical facts with new field names. The original tables, figures, papers and registry data remain with their sources. The Finnish paper carries ACOG copyright, the Swedish article uses Creative Commons Attribution–NonCommercial, and the Maastricht portal reserves rights to the Dutch publication.

## Computation and checks

For Finland and Sweden, recover approximate expected counts by dividing each observed count by its printed standardized incidence ratio (SIR). Sum the expected clear-cell and endometrioid counts, then subtract that sum from the expected all-cancer count to obtain the complement. Apply the same subtraction to observed counts. The reported contrast is the SIR for the two specified histotypes combined divided by the SIR for the complement.

Conditioning on the observed total yields a binomial model under independent Poisson counts with fixed expected counts. The exact two-sided binomial test uses the expected target share under equal SIRs. Transform the Clopper–Pearson 95% interval for the observed share to odds and divide by the reference odds. This conditions on national reference rates. It omits uncertainty in those rates. Allowing every printed SIR to range independently across its rounding interval gave point estimates from 2.77 to 2.80 in Finland and 2.57 to 2.83 in Sweden. Combining rounding uncertainty with the sampling interval gave envelopes of 1.75–4.38 and 1.45–4.85.

For the Netherlands, subtract the two histotypes from the total in each group, sum the three bins, and compute the conditional maximum-likelihood odds ratio with an exact 95% interval. Stratifying by follow-up bin using the Mantel–Haenszel estimator gave 3.69. Neither Dutch estimate adjusts for age. The original endometriosis and nevus cohorts had median entry ages of 39 and 45. A sensitivity analysis using the adenomyosis group yielded 33 of 136 cancers and a case-mix odds ratio of 1.56 (0.90–2.70) against the same nevus cohort. Diagnoses can overlap and the two comparisons share controls. The difference between the two ratios cannot establish specificity to endometriosis.

The ten-year lag places the cancers outside the immediate diagnostic workup. Ongoing surveillance, age, selection for surgery, removal of ovaries, reproductive history and tumor classification may still influence the contrast. There was no prospective registration or systematic cohort search. Intervals and p-values describe the chosen exploratory contrasts. A decisive follow-up would use a previously unexamined registry with age-specific denominators, separate endometriosis and adenomyosis coding, comparable pelvic assessment, and organ-removal histories. Estimate the absolute rate difference for each histotype as well as the ratio between histotypes.

## Reproduce

From the repository root with its Python dependencies installed, run:

    python scripts/ovarian_histotype_contrast.py
    python -m unittest discover -s tests -p test_ovarian_histotype_contrast.py

The calculation runs offline and is deterministic. The [input counts and saved output](../../evidence/ovarian-histotypes-2026-11-15) refer only to published aggregates. The recorded execution used Python 3.12.14 and SciPy 1.17.0. Input SHA-256 values are below.

| File | SHA-256 |
| --- | --- |
| standardized.tsv | aebe1e83785bb104ba9643b7b03a64a477d7519cb33275b473db3adc6213866e |
| dutch_counts.tsv | ece37382c1b87bad340a6e7eb86fe82dba741aea1d4713fe72dcbce851af4b3f |
