# GPL570 annotation audit of all 54,675 probe sets

**Outcome for [task #16](https://github.com/dovvnloading/Open-Disease/issues/16): completed platform metadata audit.** Across the entire NCBI annotation file, 42,901 probe sets have one usable symbol and one Entrez Gene ID under the project's [pre-results mapping rules](GPL570_ANNOTATION_AUDIT_2026-11-12.md). These represent 20,845 distinct Entrez IDs. A probe set remains the reporting unit for the GSE51981 analysis. This audit used platform metadata only.

## Counts and interpretation

| Exclusive category | Probe sets | Rule |
| --- | ---: | --- |
| Single symbol and single Entrez ID | 42,901 | One distinct valid token in each field |
| Multiple mappings | 2,206 | All tokens pass syntax checks, with multiple distinct symbols or IDs |
| Unmapped | 9,557 | Both fields empty or designated missing |
| Incomplete | 0 | One field missing and a single valid token in the other |
| Outside the accepted syntax | 11 | At least one field fails the conservative token rules |
| **Total** | **54,675** | **Distinct nonempty probe set IDs, including platform controls** |

Eight of the 11 rows outside the accepted syntax also have multiple valid Entrez IDs. Across *all* rows, 2,214 probe sets carry multiple valid Entrez IDs. This overlapping count is separate from the five exclusive categories above. Historical symbols with punctuation can fall outside the parser's grammar. The label records an interpretation limit, not a finding of an error at GEO.

Apply multiplicity correction across all tested complete-case probes. Retain excluded probes in the output with an explicit reason as specified in the [analysis plan](GSE51981_PREANALYSIS_PLAN_2026-11-12.md). Use the single-mapping subset for descriptive labels only. Withhold gene assignments for the other categories. Many probe sets share an Entrez ID even among the single-mapping rows. No pooling of probe estimates, selection by significance, or gene-level test is part of this policy.

## Version and source

| Field | Observed value |
| --- | --- |
| Platform | GPL570, Affymetrix Human Genome U133 Plus 2.0, Homo sapiens |
| File | [NCBI GPL570.annot.gz](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz) |
| Format | GEO-generated annotated platform, 21 table columns between the SOFT table markers |
| Fields used | `ID` (platform probe set), `Gene symbol` (Entrez symbol), `Gene ID` (Entrez identifier) |
| Annotation date in file | 2016-08-09 |
| HTTP Last-Modified | 2016-08-22 22:09:22 UTC |
| Retrieval completed | 2026-11-14 02:38:45 UTC (November 13 in New York), HTTPS 200 |
| Exact compressed bytes | 8,471,521 |
| SHA-256 of compressed bytes | `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394` |

NCBI derives these annotations by using sequence identifiers from submitted platforms to query its gene resources. The submitter's separate GPL570 record reports a June 2016 NetAffx build 35 update and uses different columns (`Gene Symbol`, `ENTREZ_GENE_ID`). The chosen archive carries the August 2016 GEO annotation. Its header supplies no genome assembly release, and this audit performed no sequence remapping. Its identifiers can differ from those used in the 2014 paper or from current reference databases. [GEO platform record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) · [NCBI annotation method](https://www.ncbi.nlm.nih.gov/geo/info/profiles.html) · [Original GSE51981 publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/)

GEO provides unrestricted public access and places no additional restrictions on use of GEO data. Submitters may retain third-party rights. The audit publishes the source, aggregate counts and original analysis code. The upstream table stays local and retains its own terms. [GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html) · [GEO FAQ](https://www.ncbi.nlm.nih.gov/geo/info/faq.html)

## Reproduction and verification

From a fresh checkout of the commit containing this report, use an empty output path. The transfer check is in bytes and applies to the compressed archive. The file name, URL and byte hash together distinguish this snapshot from the separate GPL570 family archive.

```bash
mkdir -p data/raw results/gpl570
curl --fail --silent --show-error --max-filesize 8500000 --max-time 100 \
  --output data/raw/GPL570.annot.gz \
  https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz
gzip -t data/raw/GPL570.annot.gz
python scripts/audit_gpl570.py \
  --input data/raw/GPL570.annot.gz \
  --expected-sha256 d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394 \
  --out results/gpl570/geo_annotation.json
python -m pytest -q tests/test_audit_gpl570.py
```

The scan used Python 3.12.14, with deterministic parsing and no random seed. The download used curl 8.14.1 and took 5.46 seconds. The pinned scan took under one second and used under 30 MiB peak RSS on Linux. All 24 targeted synthetic tests passed with pytest 9.2.0. An independent pass using Python's tab-delimited dictionary reader with quoting disabled checked the terminal table marker, unique probe count, all five category counts, multiple-ID count and distinct IDs among single mappings.

The parser at baseline commit `832b06a` rejected the actual annotated archive because it expected the submitter's column names. The revised parser recognizes the annotated format from its entity, checks the platform accession and annotation date, and uses the declared format's fields. The hash above was computed from the HTTPS payload on acquisition. A later match establishes byte identity with this audit. Neither the header check nor the digest authenticates individual gene annotations. Sample eligibility, confounding, and any biological reproduction require separate evidence.
