# How much of the U.S. matrisome ranking comes from one cycle phase?

**Retrospective exploratory computation, November 23, 2026.** In the scan-disjoint U.S. benchmark, cases scored higher in 77 of 120 within-phase comparisons (64.2%). The mid-secretory phase supplied 72 of the 120 comparisons. Giving the three phases equal weight yields **53.0%**. The other two phases together yield **21/48 (43.8%)**. In Rome, the direction survives a same-day comparison, but the released counts allow several different results for its statistical test. [Original benchmark, cohort definitions and sources](MATRISOME_HOLDOUT_AUDIT_2026-11-15.md) · [Complete numerical results](../../evidence/matrisome-ranking-uncertainty/results.json)

## Evidence across the three targets

The score uses 688 matrisome genes and fixed weights trained on 56 case and 27 comparator GSE51981 arrays after physical overlaps with the U.S. test series were removed. For each target, we counted how often a case scored higher than a comparator. Multiple pairs share an array. The tests below reassign **whole array labels**, retaining case and comparator counts within each U.S. phase. Rome and Mumbai each contribute one phase. All three rankings were public before this analysis. [Score specification](MATRISOME_BENCHMARK_PLAN_2026-11-15.md)

| Target | Case / comparator arrays | Higher / eligible pairs | As extreme / possible label assignments | Two-sided p | Three-target Holm p |
| --- | ---: | ---: | ---: | ---: | ---: |
| U.S., within phase | 21 / 16 | 77/120 | 174,582,066 / 943,422,480 | 0.1851 | 0.3701 |
| Rome, proliferative | 8 / 6 | 44/48 | 24 / 3,003 | 0.0080 | 0.0240 |
| Mumbai, mid-proliferative | 4 / 3 | 2/12 | 8 / 35 | 0.2286 | 0.3701 |

These probabilities assume exchangeable array labels within the indicated groups. Historical control reuse and uncertain donor independence limit the U.S. interpretation. Mumbai has only 35 allocations and unresolved technical batch overlap. The small nominal Rome probability assumes that the scanner schedule can be ignored. Holm adjustment covers these three target comparisons. It cannot account for the retrospective choice of analysis or correct confounding.

## What changes with the mix of U.S. phases?

| Phase | Cases ranked higher | Share of all eligible pairs | Ranking fraction with this phase omitted |
| --- | ---: | ---: | ---: |
| Proliferative | 16/30 | 25% | 61/90 (67.8%) |
| Early secretory | 5/18 | 15% | 72/102 (70.6%) |
| Mid-secretory | 56/72 | 60% | 21/48 (43.8%) |

The original fraction describes a uniformly selected eligible pair. For equal phase weighting, first select a phase with probability one third, then select a pair within it. This gives (16/30 + 5/18 + 56/72)/3 = 143/270. The exact two-sided sensitivity test gives p = 0.7865. Neither weighting estimates performance in a population with an unknown phase distribution.

## What can the Rome scanner dates establish?

Only 15 of Rome's 48 case–comparator pairs shared a scan date. There were three cases and three comparators on one date, two and three on another, and three cases alone on the third. Four comparisons were losses overall. Even if all four occurred within dates, cases still ranked higher in at least **11/15** same-day pairs. We enumerated every attainable same-day win total consistent with the pooled 44 wins and the scanner margins. [Native scanner audit](GENOME_TRANSPORT_2026-11-15.md) · [Scanner counts](../../evidence/genome-transport-2026-11-15/rome_scan_transport.json)

| Possible same-day wins out of 15 | 11 | 12 | 13 | 14 | 15 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Two-sided p under within-day exchangeability (200 allocations) | 0.31 | 0.17 | 0.08 | 0.03 | 0.01 |

The published totals guarantee the direction and leave the actual same-day count unknown. A nominal 0.05 decision depends on which count occurred. Every case scan followed every comparator on both mixed dates. Even an exact same-day count would leave disease status confounded with within-day processing order.

## Sensitivity to omitting individual arrays

Holding scores fixed, we bounded the ranking fraction after deleting arrays, keeping at least one case and one comparator per stratum. Each endpoint is attainable in some ordering consistent with the released counts. These are limits on what an omission could change. Identifying the influence of a particular array requires its score.

| Target | Attainable range after one omission | Fewest omissions at which the original direction could disappear |
| --- | ---: | ---: |
| U.S., within phase | 61.3%–69.4% | 4 (lower bound 41/84) |
| Rome | 90.0%–100% | 5 (lower bound 1/2) |
| Mumbai | 0%–25.0% | 2 (upper bound 1/2) |

Thus any three individual U.S. omissions preserve the above-half ranking, while removing the entire mid-secretory phase reverses it. Omitting a phase changes which population the statistic describes. These bounds are distinct from sampling confidence intervals.

## Calculation and reproduction

Within each stratum we counted strict score orderings by the Mann–Whitney win statistic, then convolved their integer frequency distributions. Tails include both directions at least as far from the permutation mean as the observation. Define equal-score order by a fixed, diagnosis-blind key such as array accession. Every reported case–comparator tie count is zero, so every such ordering gives the same observed wins. This also specifies a permutation statistic when two arrays within one group originally tied. The calculations use only the scan-disjoint fit. [Mann and Whitney](https://doi.org/10.1214/aoms/1177730491) · [Holm](https://www.jstor.org/stable/4615733)

For omissions within a stratum with H wins among P pairs and R surviving pairs, the sharp lower and upper surviving win counts are max(0, H − (P − R)) and min(H, R). The implementation combines feasible omissions across phases and recomputes the surviving denominator. The Rome scanner calculation separately enumerates compatible counts by assigning each successive rank a diagnosis and scanner block. It uses only group sizes and scanner order from the earlier genome-wide audit. Its expression outcomes belong to a different score.

The two inputs were committed in `1859b40` and are pinned by SHA-256. The [task record](../../research/tasks/matrisome-ranking-uncertainty.json) references issue #45. All computations are deterministic. The commands below reproduce the result with the project development dependencies installed.

```sh
python scripts/matrisome_ranking_uncertainty.py | diff - evidence/matrisome-ranking-uncertainty/results.json
python -m pytest -q tests/test_matrisome_{rank_inference,scan_bounds,omission,ranking_uncertainty}.py
```

The tests independently enumerate small rankings, deletions and block assignments. The analysis required two committed aggregate files and no human-level inputs. Evaluating diagnostic performance would require a decision threshold and an independently sampled cohort with outcome-blind preprocessing and a processing schedule that mixes the groups.
