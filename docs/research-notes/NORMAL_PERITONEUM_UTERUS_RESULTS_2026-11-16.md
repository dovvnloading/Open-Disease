# The five-gene peritoneum–uterus difference occurs in surgically negative controls

**Exploratory human-data computation, November 16, 2026.** The frozen five-gene score was higher in peritoneum than in the same woman’s uterine endometrium in **all eight** Finnish controls whose laparoscopy excluded endometriosis. The mean paired difference was **+29.22 within-array percentile points**. The exact two-sided paired sign-flip p-value was **2/256 = 0.0078125**. The 95% compatibility interval under a common additive shift was **+14.21 to +41.57 points**. This meets the prospectively defined criterion for a positive organ difference in controls. Visible endometriosis is unnecessary for the score to rise in this direction. [Frozen protocol at d4c6ba1](https://github.com/dovvnloading/Open-Disease/blob/d4c6ba11ddad7c5a9abf7fbaedf22d43de4e1184/docs/research-notes/NORMAL_PERITONEUM_UTERUS_PLAN_2026-11-16.md) · [Primary study](https://doi.org/10.1038/s41597-020-00623-x) · [Aggregate results](../../evidence/normal-peritoneum-uterus/results.json)

| Prespecified quantity | Observed result |
| --- | ---: |
| Independent donors / assays | 8 / 16 |
| Peritoneum-minus-uterus mean | +29.22 percentile points |
| Positive / zero / negative donor differences | 8 / 0 / 0 |
| Two-sided exact test | 2 of 256 assignments |
| 95% common-shift compatibility interval | +14.21 to +41.57 points |
| Mean after omitting each donor in turn | +27.21 to +31.99 points |
| Prespecified phase sensitivity, 6 donors | +29.53 points |
| Repeated specimens remaining for technical comparison | 0 |

The phase sensitivity includes three proliferative and three secretory uterine specimens with stored false medication status. One donor with atrophic uterine histology and one without a phase classification remain in the primary analysis. No additional test was run on this subset. Each selected specimen had one verified assay, so there was no within-specimen repeat difference to estimate. [Aggregate results and eligibility](../../evidence/normal-peritoneum-uterus/results.json)

## Who contributed the eight pairs

The native recruitment and confirmed-diagnosis fields had to agree. Among 46 healthy-control donors with at least one relevant assay, 22 supplied only uterine tissue and ten supplied only peritoneum. Fourteen supplied both. Six had no common assay run. The other eight had exactly one common run and one distinct specimen per tissue within that run. No case donors entered the selection. We froze native donor, specimen and assay keys privately before physical-chip identification. [Eligibility audit](PERITONEAL_CONTROL_AUDIT_2026-11-16.md) · [Source release](https://doi.org/10.6084/m9.figshare.12195975.v1)

The blinded chip census examined address keys across 25,355,096 source records. All 16 selected assays contained each of the 48,701 WG-6 v2 physical addresses once. None matched either HT-12 reference. No assays or donors were removed at this step. We locked the verified pairs before opening target signals. The software rechecked both private locks and the pinned archive before scoring. All 16 verified arrays had one finite submitted average signal at each expected address. The six-donor minimum was satisfied. [Census and scoring code committed before outcome access](https://github.com/dovvnloading/Open-Disease/commit/6f69a4aabc160d460338ea5cef071acbf1944124)

## What the comparison measures

On each array we ranked the full physical address set with midranks for ties, scaled as 100 × (rank − 0.5) / 48,701. We averaged three positive contributions (TMEM255A, DIRAS3, JCHAIN) and two negative contributions (GABRP, SALL1) using the original equal-magnitude signed weights. Each frozen Entrez gene mapped to exactly one verified address. The unit is the published five-gene score in percentile points, independent of the absolute microarray intensity scale. No gene, probe, weight or cutoff was selected from the controls. [Frozen panel](../../evidence/lesion-transfer/panel.json) · [Manifest and addresses](../../evidence/peritoneal-control-audit/geo-platforms.json)

The two-sided test exhausts all 256 within-donor sign assignments and includes ties. Its model assumes independent donors and exchangeability of the signs under the null. Tissue was observed rather than randomized. Inverting the same test gives the 95% compatibility interval for a common shift. Donors can differ in their actual tissue responses. A single primary endpoint was tested. The phase-restricted and deletion estimates are descriptive and were not used in the decision.

The previous lesion–uterus replication and the core–adjacent comparison now have a disease-negative anatomical reference. A whole-tissue organ or cell-composition difference is sufficient to produce the direction in these controls. The Finnish result does not estimate how much of a lesion–uterus difference is anatomical, nor distinguish cell proportions from within-cell changes. Subtracting means across these cohorts would mix platforms, body sites and sampling protocols. Positive lesion–uterus pairs therefore cannot by themselves establish lesion specificity. [Prior lesion synthesis](LESION_PANEL_SYNTHESIS_2026-11-15.md) · [Prior core–adjacent comparison](ADJACENT_PERITONEUM_RESULTS_2026-11-15.md)

The source classifies controls by inspection at laparoscopic sterilization. It lacks biopsy-specific negative pathology for their peritoneum, pelvic subsite, and a visit-date key that verifies same-day collection of paired tissues. A common run limits assay variation, but does not establish identical processing order, cycle timing or exposure. The native medication flag defaults to false. Published control-peritoneum samples include a medication category. None of the eight can be called untreated based on that flag. Only eight of 14 paired donors had same-run material, limiting transfer to other women and surgical settings. [Sources and covariate discrepancies](PERITONEAL_CONTROL_AUDIT_2026-11-16.md)

**Decision on the original biological question.** The disease-negative organ contrast is positive. The full three-arm test of an *additional* lesion effect remains a NO-GO with currently audited sources. It would require source-linked biopsy pathology, comparable pelvic subsites and sufficient donor/run support for lesions, negative case peritoneum and noncase peritoneum. A documented distant case biopsy would help distinguish a local tissue response. No diagnostic, causal or therapeutic conclusion follows. [Ranked cohort inventory and reopening criteria](PERITONEAL_CONTROL_AUDIT_2026-11-16.md)

## Reproduction and verification

The separate [task and budget](../../research/tasks/normal-peritoneum-uterus.json) were committed at [b711516](https://github.com/dovvnloading/Open-Disease/commit/b7115162bf9168b163b2252ab9012f0dd5641725). Source data are the CC BY 4.0 EndometDB v1 archive from Gabriel et al. and the pinned Illumina manufacturer annotations from GEO. This run reused the clinical ranges and WG-6 manifest. It acquired the two missing HT-12 manifests and assembled only the required archive prefix using range requests. The new source transfer was 539,966,947 bytes within the 600 MB cap. The full archive is 2,955,844,826 bytes. Private keys, patient rows, signals and run files remain outside Git. [Figshare v1 and terms](https://api.figshare.com/v2/articles/12195975/versions/1) · [Manufacturer manifest links and SHA-256 pins](GSE141549_ELIGIBILITY_2026-11-13.md) · [Exact native prefix and SHA-256](GSE141549_CHIP_CENSUS_2026-11-14.md)

From the repository root, install the repo dev dependencies and inflate64 1.0.4. In the ignored input directory below, use the filenames from the linked recipes: archive-start.bin, archive-tail.bin, GPL13376-platform.gz, GPL10558-R1.gz, GPL10558-R2.gz and raw-prefix.zip.part. The latter must be the exact 512,595,211-byte prefix with SHA-256 `74316094156c43ebb693c060178238b46c0b3f5fa032f8333488063da9e5ce7a`. The audit and scoring CLIs are offline. Use the first command to create candidate and verified private locks. On subsequent runs disagreement stops the computation.

```sh
inputs=data/private/normal-peritoneum-uterus
out=results/normal-peritoneum-uterus
mkdir -p "$out"
python scripts/normal_peritoneum_preflight.py "$inputs" > "$out/chip.json"
python scripts/normal_peritoneum_score.py "$inputs" > "$out/results.json"
diff -u evidence/normal-peritoneum-uterus/results.json "$out/results.json"
python scripts/verify_normal_peritoneum_score.py "$inputs" "$out/results.json"
python -m pytest -q
```

The tests cover donor/assay attribution, off-run and repeated specimens, chip failure, ties, the minimum-pair gate, absent outcomes, sign reversal and interval endpoints. Synthetic all-tied and identical arrays returned zero paired difference. Swapping tissue labels negated synthetic differences and preserved the two-sided tail count. Full-reference ranks agreed with independent counts of lower and tied features for every target on every retained array.

A second implementation re-read the native records selected by the same verified design. It recomputed ranks using insertion points in independently sorted signals, calculated the exact tail using subset sums, and checked the compatibility set at every distinct rational breakpoint and between successive breakpoints. It reproduced the mean, sign counts, test, interval, deletion range, phase-restricted mean and technical repeat count. This cross-check shares the cohort definition, archive decompression and physical address lock with the primary calculation. It is a software verification within the cohort. [Independent verification code](../../scripts/verify_normal_peritoneum_score.py) · [Verification result](../../evidence/normal-peritoneum-uterus/verification.json)

The blinded preflight took 19.99 s (peak RSS 589,868 KiB). The primary command, which reverifies eligibility, took 38.43 s (1,111,596 KiB). Independent verification took 35.11 s (1,118,544 KiB). CPython 3.12.14, NumPy 2.3.5, SciPy 1.17.0, inflate64 1.0.4 and pytest 9.2.0 were used. All 333 Python tests passed, with one environment-plugin deprecation warning. The procedure is deterministic and uses no random seed.
