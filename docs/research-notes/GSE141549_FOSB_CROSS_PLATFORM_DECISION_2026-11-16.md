# EndometDB FOSB: source-level cross-platform decision

**2026-11-16 · Issue #26 · Task endometdb-fosb-measurement-audit · Literature/metadata assessment**

**Decision: NO-GO for a new donor-linked analysis of the authors' pooled normalized FOSB row.** The combined GSE141549 workbook labels a single FOSB row using the HumanHT-12 v4 probe ID, but WG-6 v2 measures FOSB with a different physical probe. The inspected public materials do not identify how WG-6 measurements entered that pooled row or provide a source-backed crosswalk between GEO columns and native EndometDB assay IDs. This is a measurement-provenance blocker, **not** a null biological finding. We opened **no new Finnish FOSB expression values** in this work.

## Source-backed trace

| Finding | Source | Consequence |
| --- | --- | --- |
| HumanWG-6 v2, GPL13376, FOSB probe **ILMN_1730706**, physical address **6040086** | [Illumina WG-6 manifest](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL13nnn/GPL13376/suppl/GPL13376_HumanWG-6_V2_0_R4_11223189_A.txt.gz) and [existing physical-probe audit](GSE141549_NATIVE_PROBES_2026-11-14.md) | Its address must be interpreted on WG-6 only |
| HumanHT-12 v4, GPL10558, FOSB probe **ILMN_1751607**, physical address **7160239** | [HT-12 R1 manifest](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/suppl/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt.gz), [R2](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/suppl/GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz), [probe audit](GSE141549_NATIVE_PROBES_2026-11-14.md) | The two FOSB 50-mer oligos differ; using the same integer address across platforms would even select different genes |
| Combined workbook has **28,247** data rows; **one** FOSB row at Excel row **21,565**, labeled **ILMN_1751607**; no ILMN_1730706 row | [Annotation-only complete workbook scan and SHA-256](GSE141549_FOSB_MEASUREMENT_2026-11-13.md), [GEO supplement](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141549/suppl/GSE141549_batchCorrectednormalizedArrayscombined.xlsx) | The output annotation does **not** document the older chip's contribution to this row |
| Published processing log-transformed and quantile-normalized by chip generation, remapped probes and used ComBat to combine chips | [Gabriel et al. 2020, Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC7455745/) and [GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141549) | A symbol-level match cannot authenticate physical probe identity, log scale or batch-correction provenance |
| Native normalized-table **prefix**, not the entire table, shows repeated nonmeasurement keys, zero search/address keys, and **no FOSB record within the inspected prefix** | [Native normalized-prefix audit](GSE141549_NATIVE_PROBES_2026-11-14.md), [versioned Figshare v1 SQL archive](https://doi.org/10.6084/m9.figshare.12195975.v1) | FOSB's normalized multiplicity and a unique normalized expression value per native assay are **not verified** |
| Paper says preprocessing code is available **upon request**; inspected [author application snapshot](https://github.com/micawo/EndometDB/tree/5465ed005d0210256fb0e22e0c544f037568b460) lacks an explicit published FOSB chip-pair rule or native-to-GEO column crosswalk | [Paper Data Availability](https://pmc.ncbi.nlm.nih.gov/articles/PMC7455745/) and [public GEO inventory](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE141549) | No supported pairing was found **in the inspected public sources**; absence from every possible source is not claimed |

The project had already pinned the combined workbook SHA-256 to **9919fc5d6d6c09dbdae5a815da798b8250edf17a781fdb86dc25c2f9c9b21d50**, WG-6 manufacturer SHA-256 to **c6f2c87ab7bc094b22270c907e1f9a6d74e9990b62726b53bf44660b0db700c1**, and HT-12 R1/R2 to **770fb74bb90e548ee731924310d8f5c48589e63ee530dc11c7686d28a3ba5ba6** / **932e2f110bd30ac70ef4d7dbac2f87f8088f982bf9aa596732c7477ab50d43c2**. These are **inherited previously verified source pins**, not new downloads or newly computed hashes. They are frozen in the [aggregate machine-readable ledger](../../evidence/endometdb-fosb-measurement-audit/measurement-provenance.json).

## Distinct source namespaces and actual eligible groups

From the [native clinical eligibility audit](GSE141549_ELIGIBILITY_2026-11-13.md), **38** donors had exactly one clinically eligible assay with recorded hormone-medication false and proliferative or secretory histology. The [full native chip census](GSE141549_CHIP_CENSUS_2026-11-14.md) verified **32** unique, complete WG-6 assays (15 cases, 17 controls), after excluding six duplicate-address assays without choosing between technical records. Only these seven **native run × histologic-phase** cells contain cases and controls:

| Native run | Histology | Cases | Controls |
| --- | --- | ---: | ---: |
| 10 | Secretory | 2 | 6 |
| 11 | Proliferative | 2 | 1 |
| 12 | Proliferative | 1 | 1 |
| 12 | Secretory | 1 | 2 |
| 13 | Proliferative | 2 | 1 |
| 13 | Secretory | 1 | 1 |
| 14 | Proliferative | 3 | 1 |
| **Total (7 cells)** | | **12** | **13** |

The remaining **seven** complete WG-6 donors have no within-cell case-control comparison. Each shared cell has fewer than three donors in one or both arms. Phase is only proliferative/secretory, not fine-grained. The hormone survey **defaults false** and has **no known lookback**; it does not verify washout. Controls underwent tubal-ligation surgery and were negative for endometriosis at laparoscopy, but this does not guarantee absence of every pathology; recruitment differs from cases. [Primary methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC7455745/) · [clinical audit](GSE141549_ELIGIBILITY_2026-11-13.md).

Separately, in **GEO sample-title** metadata, the previously audited primary HT-12 endometrium grouping contains **39 case and zero control arrays**. This is not a donor count and cannot support a platform-specific case/control test. **Do not connect GEO columns to native donors by phenotype, title, age, or other quasi-identifiers.** [GEO platform audit](GSE141549_ELIGIBILITY_2026-11-13.md).

## Locked decision and discrimination rule

1. **Blocked estimand:** no case-minus-control estimate of the pooled, ComBat-corrected FOSB row for the 25 native run-by-phase-matched WG-6 donors. The row origin and donor-to-column key are unresolved. No invented mapping, apparent probe alias or cross-platform gene symbol may bypass this stop.
2. **Separate experiment:** the native raw WG-6 **whole-array rank-percentile** question in [issue #32](https://github.com/dovvnloading/Open-Disease/issues/32) was frozen at [44dd06b](https://github.com/dovvnloading/Open-Disease/commit/44dd06b2526540b29dcb54801d61a9f9fae68e94) before its own expression analysis and has [its own report](GSE141549_WG6_FOSB_RESULTS_2026-11-14.md). Its measured direction must not be promoted to validation of pooled normalization. This new decision cannot be described as prospectively preceding that **earlier** result.
3. **Future GO prerequisites:** an inspectable author or source-derived, platform-specific WG-6/HT-12 FOSB pairing and transformation, uniquely identified normalized values by native assay ID (with authoritative mapping of any GEO columns), plus auditable phase/run counts and a prespecified independent comparison. Document treatment-window and recruitment caveats. Any new source requires a fresh, *separately committed* design before opening a new outcome.
4. **Fail-closed guard:** the companion Python program only validates existing aggregate, pinned observations. If a new provenance URL is added, it returns **RE_REVIEW_REQUIRED**, **never GO** based merely on that string.

## Offline reproduction

From the repository root:

    python scripts/verify_endometdb_fosb_provenance.py evidence/endometdb-fosb-measurement-audit/measurement-provenance.json --require-no-go
    python -m pytest -q tests/test_verify_endometdb_fosb_provenance.py
    python scripts/research_queue.py validate
    python scripts/research_queue.py ready

**Execution:** CPython 3.13.5, pytest 9.0.2, **19/19 offline synthetic aggregate tests passed**; deterministic, no random seed, no participant records, no raw/derived expression files, no external downloads or new source-data inspection. Whole-repository CI and tests must be checked on the PR. No clinical result, treatment inference, outreach, paid resource or new Finnish test was produced.
