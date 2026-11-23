# Menstrual outflow obstruction and recorded endometriosis

**Exploratory computation from published aggregate human data, 23 November 2026.** Endometriosis was more common in the anatomical group expected to retain more menstrual flow in both Beijing and London. The observed differences were 25.1 and 12.6 percentage points, respectively. Differential detection remains a plausible explanation, especially in London, where horn excision was 48.3 points more frequent in women with a functional horn. Neither cohort measured menstrual reflux or the onset of endometriosis.

## Published counts and project calculations

Tong and colleagues had already reported the Beijing association. Tellum and colleagues had reported an association with a functional rudimentary horn and published counts by anatomical subtype. We used only the two **noncommunicating** horn groups in London, compared three diagnoses in those same women, and calculated how many missed cases would erase each endometriosis difference. The source reports and their directions were known when the project specified the analysis. [Project protocol and decision rules](https://github.com/dovvnloading/Open-Disease/issues/60)

| Cohort and outcome | Higher expected retention, events / women | Reference, events / women | Difference in percentage points (95% interval) | Additional reference diagnoses to erase excess |
| --- | ---: | ---: | ---: | --- |
| Beijing · pelvic endometriosis | 10/27 (37.0%) | 8/67 (11.9%) | +25.1 (+6.7 to +44.7) | 17 of 59 recorded noncases (28.8%) |
| London · endometriosis | 23/88 (26.1%) | 17/126 (13.5%) | +12.6 (+1.9 to +23.8) | 16 of 109 recorded noncases (14.7%) |
| London · adenomyosis | 24/88 (27.3%) | 36/126 (28.6%) | −1.3 (−13.1 to +11.1) | — |
| London · fibroids | 15/88 (17.0%) | 22/126 (17.5%) | −0.4 (−10.3 to +10.4) | — |
| London · horn excision | 46/88 (52.3%) | 5/126 (4.0%) | +48.3 (+36.9 to +58.7) | — |

**Beijing.** Tong et al. reviewed clinical and operative records and used follow-up contact for a series of 94 women with Herlyn–Werner–Wunderlich syndrome at Peking Union Medical College Hospital. The groups had complete and incomplete hemivaginal obstruction, respectively. The abstract gives both counts, but leaves the recruitment period and the intensity of pelvic investigation within each group unspecified. An earlier 70-woman series from the same team and hospital could overlap and contributes no additional observations here. [Tong et al., 2014, abstract, Results](https://pubmed.ncbi.nlm.nih.gov/25044087/) · [Related 2013 series, abstract, Methods](https://pubmed.ncbi.nlm.nih.gov/23465857/)

**London.** Tellum et al. reviewed consecutive women aged 16 or older with a unicornuate uterus who attended University College London Hospital between January 2008 and September 2021. Table 3 gives the group sizes in the column headers and the counts in the Endometriosis, Adenomyosis, Fibroids and Excision of rudimentary horn rows. We compared functional with nonfunctional noncommunicating horns. Women with an absent or communicating horn and the separately matched women with a normally shaped uterus were outside this comparison. The study used systematic specialist ultrasound and surgical histories. The ultrasound criteria describe nodules and ovarian endometriomas, which can leave other lesions undetected. The authors also excluded women whose fibroids obscured the uterine cavity. [Tellum et al., Methods and Table 3](https://obgyn.onlinelibrary.wiley.com/doi/full/10.1002/uog.26073)

## Decision and sensitivity to missed diagnoses

Both endometriosis differences were positive and the London endometriosis difference exceeded the point estimates for both comparison diagnoses. The study therefore meets the protocol's descriptive predictions. The comparator intervals allow meaningful differences in either direction. Adenomyosis and fibroids differ in biology and detectability, and counts of joint diagnoses are unavailable. The three London outcomes cannot support a formal test of differences between the contrasts from these margins alone.

The last column fixes both group sizes and treats recorded cases as true. It adds missed endometriosis diagnoses only in the reference group. In Beijing, adding 16 cases leaves a **+1.22** point difference and adding 17 gives **−0.28** points. In London, adding 15 leaves **+0.74** points and adding 16 gives **−0.05** points. These are the smallest counts that erase the respective point estimates under those assumptions. How many cases went undetected in either group remains unknown.

Horn excision is evidence of different care pathways. It can give clinicians an opportunity to inspect the pelvis, but its frequency alone cannot reveal which diagnoses surgery detected. Selection into a surgical or ultrasound clinic, symptoms prompting referral, time since menarche, and time to anatomical correction could all influence the contrast. These two targeted cohorts characterize different anomalies and diagnostic practices. The associations remain compatible with a role for retained flow and with differential detection or selection. An independent comparison with the same prespecified endometriosis assessment in both anatomical groups, and a record of who received that assessment before surgery, would directly address the detection explanation.

## Offline reproduction

The calculation subtracts the reference prevalence from the higher-retention prevalence. The two-sided 95% intervals use Newcombe's method 10 for independent binomial proportions, combining uncorrected Wilson score limits. They quantify sampling uncertainty within each comparison. They do not account for selection, measurement error or multiple comparisons. The threshold calculation uses exact integer arithmetic. All five intervals were checked using independently calculated SciPy Wilson limits, and both thresholds by enumerating every feasible number of missed cases. [Interval method](https://pubmed.ncbi.nlm.nih.gov/9595617/) · [Source-attributed counts](../../evidence/outflow-obstruction-endometriosis/counts.json) · [Calculated results](../../evidence/outflow-obstruction-endometriosis/results.json)

Run from the repository root with Python 3.12.14 or another compatible Python 3 installation. The calculation uses only the standard library, runs in under a second and is deterministic. The aggregate input SHA-256 is `e87ec4c3ba1b18bd055e9f79c137bdc212a6635ecb644c5d9cc3aa6e402c4911`.

```sh
python evidence/outflow-obstruction-endometriosis/calculate.py | diff -u evidence/outflow-obstruction-endometriosis/results.json -
```
