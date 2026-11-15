# Does inherited susceptibility support the frozen endometriosis candidates?

**Exploratory secondary analysis of published human aggregate results, November 15, 2026.** Among the 125 unique genes in the five nominated genes and five frozen Reactome sets, ABO alone passed the published genome-wide multi-tissue transcriptome-wide association study (TWAS) threshold. Genetically predicted ABO expression had opposite signs in whole blood and mammary tissue. The uterus model gave weak evidence. None of the 125 genes appeared in the study's reported endometriosis–expression colocalization signals. In a separate Danish hospital cohort with 6,640 diagnosed endometriosis events, measured blood group A versus O had an incidence rate ratio of 1.05 (study-wide adjusted 95% interval 0.92–1.20). These endpoints locate the main evidence at the inherited ABO locus. They leave its mechanism and relevance to uterine expression unresolved.

## Frozen genes and inherited susceptibility

Open-Disease selected the five genes and five projected Reactome pathways using Finland before opening the Rome expression results. The genetic and registry analyses documented here were designed after the earlier results, and some ABO literature was already known. Koller et al. pooled eight cohorts (105,869 cases and 1,282,731 controls overall). We intersected the original Entrez IDs with published expression-QTL and multi-tissue TWAS tables by version-stripped Ensembl gene ID, using unique HGNC mappings for all 125. The official symbols also agreed with NCBI Gene. Koller included the earlier International Endometriosis Genetics Consortium and UK Biobank results. The 2023 and 2026 GWAS findings should not be counted as independent replications. [Original nominations](https://github.com/dovvnloading/Open-Disease/blob/f2a3755ff055cde241e0eaaaedc11a1471270beb/evidence/cross-cohort-2026-11-14/nominations.json) · [Koller et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC13195585/) · [Rahmioglu et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10042257/) · [HGNC](https://hgnc.genenames.org/download/)

| Original set | Locked genes | Multi-tissue models | Uterus models | Smallest multi-tissue p |
|---|---:|---:|---:|---:|
| Interleukin-37 signaling | 19 | 12 | 2 | 0.0653 |
| Signaling by FGFR3 | 38 | 25 | 5 | 0.00516 |
| Triglyceride catabolism | 21 | 13 | 3 | 0.0496 |
| PI-3K cascade downstream of FGFR3 | 18 | 11 | 2 | 0.0139 |
| Golgi-associated vesicle biogenesis | 43 | 32 | 4 | 0.0639 |
| Five individual genes | 5 | 5 | 2 | 1.68 × 10⁻¹² (ABO) |

The pathway rows overlap. They contain 120 distinct genes, 82 of which had a multi-tissue model. Table 8 lists only positive colocalizations, with 29 rows covering 15 genes and 13 candidate variants (full-colocalization posterior >0.7 and SNP posterior >0.5). Zero matched the locked genes. A gene absent from that hit table might have lacked informative cis-eQTLs, expression or power. The table does not provide a complete tested-gene denominator, and its 29 rows include only one uterus row. Table 13 supplies 13,457 multi-tissue models, including 87 locked genes. At the authors' threshold of 3.72 × 10⁻⁶, 41 genes passed across the full table, including ABO. All five nominated genes had models. The next smallest p among the five was 0.0198 for KIF3C. This is a gene lookup, without an omnibus pathway-enrichment test or a fitted causal model.

The omnibus ABO p does not encode direction. The seven marginal tissue-specific z-scores include mammary tissue +6.02 (p = 1.73 × 10⁻⁹), whole blood −5.45 (p = 4.90 × 10⁻⁸) and uterus +1.55 (p = 0.122). The direction in the latter is imprecise. The earlier 2023 blood eQTL SMR result at rs550057 had p = 5.35 × 10⁻⁷, with HEIDI p = 0.0437 passing that study's multiplicity-adjusted criterion. The 2026 multi-ancestry index variant rs635634 had p = 7.92 × 10⁻¹⁴. In its clinical endometriosis analysis Table 11 lists five rs635634/methylation-probe colocalizations. Four had reported PPFC 0.999, and a fifth had 0.9923. Genotype, genetically predicted transcript, DNA methylation and assayed biopsy RNA represent separate quantities. Neither paper establishes the causal cell population. [2026 supplementary tables](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-026-02582-2/MediaObjects/41588_2026_2582_MOESM4_ESM.xlsx) · [2023 supplementary tables](https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-023-01323-z/MediaObjects/41588_2023_1323_MOESM4_ESM.xlsx)

## Blood group as a separate observable prediction

A Yale surgical case-control series reported a blood group A odds ratio of 2.89 (1.85–4.52). Geography, recruitment, outcome definitions, exposure reference and effect measure matter here. We extracted the endometriosis row from a Danish study with 482,914 blood-typed hospital patients, including 6,640 endometriosis events and 10,418,467 person-years in the female phecode 615 analysis. Supplementary file 7 gives the comparison to blood group O. Rates were adjusted for age, sex, age–sex interaction and birth year in the source modeling framework. The confidence intervals and p-values below are the authors' FDR-adjusted values. [Yale surgical series](https://pubmed.ncbi.nlm.nih.gov/19301022/) · [Danish cohort, data and methods](https://elifesciences.org/articles/83116)

| Blood group | Reference | Events in exposure group | Incidence rate ratio | Adjusted 95% interval | Adjusted p |
|---|---|---:|---:|---:|---:|
| A | O | 2,894 | 1.05 | 0.92–1.20 | 0.462 |
| B | O | 778 | 0.99 | 0.54–1.80 | 0.970 |
| AB | O | 284 | 0.94 | 0.56–1.59 | 0.836 |
| RhD positive | RhD negative | 5,618 | 0.99 | 0.63–1.55 | 0.956 |

The alternate A-versus-all-other-ABO contrast in the same people had IRR 1.059 (nominal p = 0.0119, study FDR q = 0.196). Both contrasts give a small point estimate. The Danish A-versus-O interval excludes large increases in this hospital population and accommodates a modest increase. Its retrospective design required survival to in-hospital blood typing, which limits inference about community incidence. Adding the Danish estimates into a meta-analysis with Yale would mix effect measures, comparator definitions and clinical ascertainment.

A UK Biobank lookup provided a distinct genotype-defined measurement. Li and Schooling's Table S2 used an additive linear-probability model for 1,516 main-diagnosis ICD-10 N80 cases and 359,678 controls. The combined-sex denominator limits interpretation as risk per woman. The A1 tag rs507666 had nominal p = 0.0154 (four-tag Holm p = 0.0617); A2, B and O tags all had p > 0.26. These UK participants may contribute to the GWAS meta-analyses. A serologic A group combines subtypes, so its estimate is not an A1 allele test. The Swedish SCANDAT3 paper listed N80 among its potential outcomes, while the figures aggregate 3,862 nonsignificant points into unlabeled traces and its validation table only lists selected associations. We retained it as a source-availability check and did not assign an endometriosis effect estimate by plot position. [UK Biobank antigen-tag study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7669452/) · [Swedish disease-wide study](https://elifesciences.org/articles/65658)

## Reproducing the aggregate lookup

Use Python 3.11+ and openpyxl (executed with 3.1.5). Fetch the four files specified in `genetics_sources.json` to a private input directory using the filenames shown. Run from the repository root. The [source manifest](../../evidence/locked-genetics-2026-11-15/genetics_sources.json) is beside the pinned crosswalk. The script verifies input digests before producing the aggregate JSON. The ten nominated sets were frozen before genetic lookup. The present comparisons were exploratory, and all signs and non-hits are retained.

```bash
python scripts/reconcile_genetic_evidence.py --source-dir data/private/genetics-lookup \
  --lock evidence/cross-cohort-2026-11-14/nominations.json \
  --crosswalk evidence/locked-genetics-2026-11-15/locked_gene_crosswalk.tsv --out results/genetics-check.json
```

Rendering these supplemental data does not fit a new GWAS, TWAS, eQTL or registry model. Raw participant records were unnecessary. The analytic contribution is the cross-source comparison of the originally frozen expression candidates, genetic evidence tiers, and a separate directly measurable prediction.
