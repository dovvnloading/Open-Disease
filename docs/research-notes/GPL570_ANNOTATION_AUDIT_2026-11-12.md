# GPL570 probe annotation policy — first pass

**Date:** 2026-11-12  
**Task:** `gpl570-annotation-audit` / [#16](https://github.com/dovvnloading/Open-Disease/issues/16)  
**Evidence level:** source-audit and synthetic software test. **No full GPL570 platform file retrieved or counted in this session, no patient expression matrix analyzed.**

The [completed real-source audit](GPL570_ANNOTATION_RESULTS_2026-11-13.md) records the subsequent release choice, checksum and results. The remainder of this note preserves the policy and status at the time it was written.

## Platform facts and mapping hazard

- [GEO GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) is the Affymetrix Human Genome U133 Plus 2.0 array. Its table includes distinct `ID` (probe set identifier), `Gene Symbol`, `ENTREZ_GENE_ID`, and source/annotation-date fields. The observed platform page includes an **October 6, 2014** example annotation date. This is *not* a verified download timestamp or statement that all annotations have a single release date.
- The [displayed GPL570 example row `1007_s_at`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) shows `DDR1 /// MIR4640` in the Gene Symbol field and multiple Entrez IDs. A probe set is not always a unique gene, and a historical gene symbol is not a stable modern causal target.
- The original [GSE51981 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/) used 54,675 probe sets. That is a **publication/platform context**, not the count produced by Open-Disease's current tool.
- [NCBI platform table guidelines](https://www.ncbi.nlm.nih.gov/geo/info/platform.html) describe `ID` and sequence/annotation provenance. Probe set IDs, Entrez Gene identifiers and gene symbols have different semantics; none establish causality.

## Conservative, pre-results policy

1. Preserve probe-level statistics as the primary analysis and multiplicity family. Do not discard significant or nonsignificant unmapped probes silently.
2. Use a *single pinned copy* of GPL570 platform annotation, record official download URL, retrieval UTC, SHA-256, file size, column schema and source annotation date(s) before evaluating gene-level results.
3. A unique syntactically valid symbol **and** unique positive integer Entrez ID is provisionally a one-to-one identifier link, **not a verified stable gene build or evidence of causality**. Multi-symbol, multi-Entrez, partial, missing and malformed entries remain separate categories.
4. No gene-level aggregations, lowest-p-value-per-gene selection, or candidate ranking until a separate documented and frozen gene-aggregation policy is selected. Multiple probes per same gene must not be mistaken for independent gene-level support.
5. If annotation is insufficient or ambiguous, report the probe accession, uncertainty and source rather than assigning a convenient nearest gene. No biological conclusions based on annotation alone.

## Offline audit tool

`python scripts/audit_gpl570.py --input data/raw/GPL570_platform_table.tsv --out results/gpl570/annotation_audit.json`

The input must already have been acquired with appropriate source checks. The command **does not download anything**. It accepts local TSV or gzip GEO SOFT platform tables with `ID`, `Gene Symbol` and `ENTREZ_GENE_ID`, rejects duplicate/malformed probe rows, and emits aggregate-only mapping category counts, compressed/source-file SHA-256, bytes and an explicit `source_verified: false` field. This flag is intentionally **not upgraded by a parser**: the operator/agent still has to verify the platform release. Output is local/ignored and refuses overwrite. Synthetic tests use fake probes and fake symbols only.

**Status / next gate:** policy and parser written; *actual GPL570 platform download, terms, complete row/annotation-date audit and versioned provenance are pending*. Do not close #16 or unblock the human-data reproduction simply because synthetic tests pass.

## Incremental offline provenance checks (2026-11-12)

The offline parser now rejects a GEO SOFT input that declares a platform other than `GPL570`, repeats `^PLATFORM`, or declares a platform without `!platform_table_begin`. The aggregate report includes `declared_platform_id` (`null` for a bare TSV or undeclared SOFT fragment) and `soft_platform_declaration_checked`. These are **internal file consistency checks**, not authentication of a public download.

After the source and its archived bytes have been independently verified, the operator can pin the SHA-256:

```bash
python scripts/audit_gpl570.py \
  --input data/raw/GPL570_platform_table.tsv.gz \
  --expected-sha256 "$VERIFIED_SHA256" \
  --out results/gpl570/annotation_audit.json
```

The optional expected digest must contain 64 hexadecimal digits for the **exact original input bytes** (compressed bytes when `.gz`). An invalid or mismatched pin causes an error before a report is written. A matching digest does **not** establish data rights, original GEO source URL, official annotation release or genome build; `source_verified` still remains `false`. Record the independently checked accession, download URL, retrieval date, annotation release/date, reuse terms and checksum separately.

**Execution evidence:** `python -m pytest -q tests/test_audit_gpl570.py` → 13 passed on synthetic-only fixtures with Python 3.13; `python -m compileall -q scripts/audit_gpl570.py` → passed. No real GPL570 annotation was retrieved, parsed or counted. **Issue #16 remains open** and reproduction #18 is still gated until actual source/release/rights and full annotation verification.
