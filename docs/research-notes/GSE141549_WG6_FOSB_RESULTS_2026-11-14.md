# Finnish FOSB comparison across matched run and cycle phase

**November 14, 2026 UTC · Exploratory human-data computation · Positive observed direction, inconclusive test.** Among 12 Finnish women with confirmed endometriosis and 13 controls, the case group’s FOSB signal ranked **10.89 percentile points higher** on average across seven shared run-by-phase strata. The exact two-sided permutation p-value was **0.1591**. The comparison did not meet the prospectively specified 0.05 threshold for directional support. [Frozen plan at commit 44dd06b](https://github.com/dovvnloading/Open-Disease/blob/44dd06b2526540b29dcb54801d61a9f9fae68e94/docs/research-notes/GSE141549_WG6_FOSB_PLAN_2026-11-14.md) · [Issue #32](https://github.com/dovvnloading/Open-Disease/issues/32)

| Frozen aggregate | Result |
| --- | ---: |
| Donors with one clinically eligible assay | 38 |
| Assays with full, uniquely verified WG-6 address sets | 32 (15 cases, 17 controls) |
| Verified donors in shared strata | 25 (12 cases, 13 controls) |
| Exact shared analysis-run and histologic-phase strata | 7 |
| Complete finite physical addresses on each included array | 48,701 |
| Equally weighted case-minus-control difference | +10.891045 percentile points |
| Assignments at least as extreme, both tails | 1,925 of 12,096 |
| Exact p-value | 275/1,728 (0.15914352) |
| Effect after omitting each stratum in turn | +6.379347 to +14.961135 percentile points |

Six of the 38 clinically eligible assays had repeated physical addresses. Seven additional donors on verified WG-6 assays fell outside shared run-by-phase strata. The native joins connected each included assay to a distinct patient. All 25 included arrays had exactly one finite submitted signal for every WG-6 address, including FOSB probe ILMN_1730706 at address 6040086. No quantitative failures or additional exclusions occurred. The target percentile used the complete physical probe reference and included half of all ties. [Clinical eligibility](GSE141549_ELIGIBILITY_2026-11-13.md) · [Complete chip census and source pins](GSE141549_CHIP_CENSUS_2026-11-14.md) · [Native release](https://doi.org/10.6084/m9.figshare.12195975.v1)

Every leave-one-stratum-out estimate retained the positive sign. This deletion range describes sensitivity to including each stratum. The exact sharp-null test assumes disease labels are exchangeable within run and phase. The small, observational surgical cohort supplies limited precision, and several strata contained only one donor in an arm. The frozen test supplies no interval for a heterogeneous mean effect. The signal percentile depends on the rest of the chip and has no fold-change units. Cases and controls had different recruitment pathways. The survey field records false for hormone medication, defaults to false, and supplies no washout interval. Site, precise cycle timing and technical effects remain potential explanations. The Finnish result alone establishes neither a disease-wide association nor equivalence. [Original study](https://doi.org/10.1038/s41597-020-00623-x) · [Prespecified interpretation](GSE141549_WG6_FOSB_PLAN_2026-11-14.md)

## Independent calculation and reproduction

The primary program repeated the blinded chip census, checked all selected signal records, counted values strictly below and tied with FOSB, and enumerated integer-scaled label assignments. The separate verifier sorted the signal arrays to recover tied ranks and convolved exact stratum-specific contrast distributions. Both programs read the pinned native archive, and share its release decoding and the native clinical joins. They agreed on the estimated effect, both ends of the deletion range, and all 1,925 extreme assignments. This arithmetic cross-check uses the same Finnish cohort and source. [Primary analysis](../../scripts/analyze_endometdb_wg6_fosb.py) · [Independent verifier](../../scripts/verify_endometdb_wg6_fosb.py)

From the repository root, reuse the ignored, previously verified source files from the census. Both scripts accept the same optional single-file ZIP prefix in place of the cached segments. The aggregate source inventory is below. No individual data or quantitative table is needed in Git.

| Input | SHA-256 |
| --- | --- |
| Native ZIP prefix, inclusive bytes 0–512,595,210 | `74316094156c43ebb693c060178238b46c0b3f5fa032f8333488063da9e5ce7a` |
| Native ZIP clinical head | `023a9ee95ecb6901cac61f1fdf82ed699e1c8d2159e61665d2ef499ebd8708fa` |
| Native ZIP clinical tail | `75f0e4ce7f600af1716c16b5b551331725d437f0608b92dc4544024f2a351fc8` |
| HumanWG-6 v2 manifest | `c6f2c87ab7bc094b22270c907e1f9a6d74e9990b62726b53bf44660b0db700c1` |
| HumanHT-12 v4 manifest R1 | `770fb74bb90e548ee731924310d8f5c48589e63ee530dc11c7686d28a3ba5ba6` |
| HumanHT-12 v4 manifest R2 | `932e2f110bd30ac70ef4d7dbac2f87f8088f982bf9aa596732c7477ab50d43c2` |

```bash
python scripts/analyze_endometdb_wg6_fosb.py data/private/gse141549
python scripts/verify_endometdb_wg6_fosb.py data/private/gse141549
python -m pytest -q tests/test_analyze_endometdb_wg6_fosb.py
```

The runs used CPython 3.12.14, inflate64 1.0.4, NumPy 2.3.5 (verifier) and pytest 9.2.0. No random seed applies. The primary calculation took 34.27 seconds with 1,189,876 KiB peak resident memory, and the independent check took 15.34 seconds with 683,444 KiB peak. Both ran offline using the census cache. The 15 dedicated synthetic tests passed, including tied ranks, finite zero and negative signals, missing or duplicate values, singleton strata, exact tail ties, both directional decisions, and frozen cohort counts. [Synthetic tests](../../tests/test_analyze_endometdb_wg6_fosb.py) · [Reviewable source and results](https://github.com/dovvnloading/Open-Disease/pull/29)
