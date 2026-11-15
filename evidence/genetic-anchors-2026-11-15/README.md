# Genetic-anchor evidence and execution

[Study results and interpretation](../../docs/research-notes/GENETIC_ANCHOR_TRANSFER_2026-11-15.md) · [Initial design](design.md) · [Follow-up order](followup.md)

The primary molecular panels come from Rahmioglu (2023) Table 1 and Koller (2026) Table S13. The former definition file condenses duplicate gene rows and preserves every middle-column gene; no locus labels or chromatin-interacting neighbors were added. All 27 externally defined gene identifiers were checked against the 2026-10-06 HGNC snapshot. HCG4P6 maps to HCG4B (80868), and SRP14-AS1 to SRP14-DT (100131089). For exact-symbol Melbourne matching we used current and paper names only. Other historical array aliases may be missed. [Rahmioglu Table 1](https://pmc.ncbi.nlm.nih.gov/articles/PMC10042257/) · [Koller et al.](https://doi.org/10.1038/s41588-026-02582-2)

| File | Contents |
| --- | --- |
| rahmioglu_table1.tsv | Source set and NCBI identifiers for all 23 genes |
| koller_crosswalk.tsv | NCBI identifier mapping for the four S13 hits |
| rahmioglu_expression.tsv | Finnish/Rome point estimates, per-gene exact tests, same-direction partial-conjunction tests, Holm over all 23 (eight unassayed), Melbourne coverage |
| koller_uterus.tsv | All significant uterus TWAS genes, signs, all relevant measurement gaps, and Melbourne two-sided tests |
| omnibus.json | Primary Rome self-contained panel tests with Holm correction over two, Finnish direction summary, diagnostic deletion/restriction |
| scan_day.json | Within-scan-date test on the eleven Rome biopsies on case-control shared scan dates, Holm over two primary tests, and competitive sensitivity |
| transfer.json | Genome-background check, Mumbai RNA-seq projection, gene-total and transcript-first measurements |
| additional_rnaseq.json | GSE315857 descriptive transfer; same source as GSE272606, counted once |

The independent statistical checker reuses the source cohort readers and constructs each possible contrast directly from group means within dates, without using the primary permutation-design/statistic implementation. It checks six exact tails. Gene-specific missing values have blank estimates. The 2023 and 2026 gene hypotheses and the later diagnostic checks were not one prespecified experiment. Holm correction applies only to the stated families; there is no study-wide prospective discovery claim.

## Inputs

Use the two [Rome source files](../../docs/research-notes/GSE25628_FOSB_RESULTS_2026-11-13.md), [Melbourne Table S20 workbook](../../docs/research-notes/MELBOURNE_FIVE_GENE_RESULTS_2026-11-14.md), the first three [GSE315857 source files](../../docs/research-notes/GSE315857_FOSB_RESULTS_2026-11-14.md), and [Mumbai and HGNC receipts](../external-bulk-2026-11-14/mumbai-source.json). Place them under the local paths shown below. The native Rome scan audit supplies the GEO brief SOFT and 22 two-kilobyte CEL prefixes under data/private/matrisome-audit. The date reader checks the SOFT digest and the full ordered prefix-series digest. [Scan-date reconstruction](../../scripts/rome_scan_dates.py)

Download the 25,743,840-byte Koller supplement to a local ignored directory. The script checks its SHA-256 c021b6c20c668644939fd83381307a2c8cf08f18d95156b33e65fb0301997815, table title, thresholds, and source coverage. The publisher workbook has 34 sheets; 13,457 rows contained results and 2,533 had measured uterine z scores in sheet 13. The repository’s frozen Finnish screen digest is 9df65b29615ec56767e9abaaf3844f8ebde5e21baf8cdab193ac24857a170da5. On this run ENA served the same seven Mumbai data records in a different sequence. The reconstruction below permutes the seven intact lines and accepts only the original source SHA-256. It rejects any content change. [Publisher supplement](https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41588-026-02582-2/MediaObjects/41588_2026_2582_MOESM4_ESM.xlsx)

```bash
p=data/private
mkdir -p "$p/koller-2026" "$p/gse153739"
curl -fLsS --max-filesize 26000000 'https://media.springernature.com/full/springer-static/esm/art%3A10.1038%2Fs41588-026-02582-2/MediaObjects/41588_2026_2582_MOESM4_ESM.xlsx' -o "$p/koller-2026/tables.xlsx"
curl -fLsS --max-filesize 5000 'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA643835&result=read_run&fields=run_accession,experiment_accession,sample_accession,secondary_sample_accession,instrument_model,library_layout,submitted_ftp,fastq_ftp,first_public,library_name&format=tsv' -o "$p/gse153739/current-ena.tsv"
python scripts/restore_mumbai_ena_order.py "$p/gse153739/current-ena.tsv" "$p/gse153739/ena-run-metadata.tsv"
python scripts/fetch_external_bulk_sources.py --fetch
python scripts/audit_matrisome_holdout.py "$p/matrisome-audit" --fetch --include-rome
```

## Offline execution

Run from the repository root using the project environment. Outputs go into an ignored results directory. Only the fixed source panels need to be copied there to run staged checks. The report provides source links for the study-specific human-data conditions and sample-selection details.

```bash
export OPENBLAS_NUM_THREADS=1
p=data/private
src=evidence/genetic-anchors-2026-11-15
out=results/genetic-anchors
mkdir -p "$out"
cp "$src/rahmioglu_table1.tsv" "$src/koller_crosswalk.tsv" "$out/"
fin=evidence/cross-cohort-2026-11-14/discovery_gene.tsv
hgnc="$p/gse266265/hgnc_complete_set_2026-10-06.tsv"
python scripts/evaluate_genetic_anchors.py --finland "$fin" --rome "$p/gse25628" --koller "$p/koller-2026/tables.xlsx" --melbourne "$p/fung-2018/Fung_supplementary_tables.xlsx" --definitions "$out" --out "$out"
python scripts/check_genetic_anchor_transfer.py --finland "$fin" --rome "$p/gse25628" --mumbai "$p/gse153739" --hgnc "$hgnc" --panel "$out" --out "$out/transfer.json"
python scripts/check_scan_day_genetic_anchors.py --rome "$p/gse25628" --scans "$p/matrisome-audit" --finland "$fin" --panel "$out" --out "$out/scan_day.json"
python scripts/describe_twas_transcripts.py --directory "$p/gse315857" --gencode "$p/gse153739" --hgnc "$hgnc" --panel "$out" --out "$out/additional_rnaseq.json"
python scripts/verify_genetic_anchors.py --rome "$p/gse25628" --scans "$p/matrisome-audit" --mumbai "$p/gse153739" --hgnc "$hgnc" --panel "$out"
```
