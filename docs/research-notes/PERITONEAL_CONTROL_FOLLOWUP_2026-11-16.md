# Does a public cohort resolve the residual lesion effect after controlling for anatomy?

**Literature, metadata and physical-assay audit, November 16, 2026.** The normal organ contrast is established in the eight surgically negative Finnish donors studied. What remains open is whether there is an additional change inside a lesion when both tissues come from comparable peritoneal sites. We rechecked the case and control biopsies, evaluated additional primary cohorts and independently traced the assays in the two strongest whole-transcriptome collections. The full pathology-controlled question remains **NO-GO** within the audited sources. [Previous result](NORMAL_PERITONEUM_UTERUS_RESULTS_2026-11-16.md) · [Earlier ranked search](PERITONEAL_CONTROL_AUDIT_2026-11-16.md)

## Why the thirteen German pairs cannot adjudicate lesion specificity

In Erlangen (E-MTAB-694), the original specimen manifest identifies 13 donor-matched lesion–peritoneum pairs within 27 case donors. The primary study reports histological review. A secondary clinical table supplies five proliferative and eight secretory labels for the 13 under the documented study-number convention. All five panel genes have qualifying probes among the 54,675 features. These promising attributes led us to inspect the physical acquisition metadata. [Primary study](https://doi.org/10.1177/1933719112451147) · [Clinical audit, supplement recovery and limitations](SOHLER_CLINICAL_SOURCE_2026-11-16.md) · [Official public study](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-694)

| Physical scanner | Lesion arrays | Case peritoneum arrays | Uterine arrays |
| --- | ---: | ---: | ---: |
| Instrument used for lesions | 18 | 0 | 17 |
| Instrument used for comparison peritoneum | 0 | 22 | 0 |

**All 13 paired lesions were recorded 373–386 calendar days before their comparison peritoneum.** Those dates come from acquisition headers. Surgical dates were unavailable for this check. There were two instruments, five scanner-days and zero scanner-days containing both lesion and case peritoneum. Uterine arrays cannot connect the instruments because every uterine array was on the lesion instrument. Donor pairing removes stable donor effects but leaves tissue completely confounded with scanner and acquisition period. Within-array ranks could remove a uniform monotone scanner transformation. There is no shared tissue measurement from which to verify that assumption or estimate probe-specific differences. [Header census and probe mapping](../../evidence/peritoneal-followup/sohler-scans.json)

To check whether pairs crossing days could collectively rescue the design, we enumerated all 8,192 ways to swap labels within the 13 donors. Only the observed assignment preserves the recorded lesion and peritoneum counts in every scanner-day. Zero pairs can move, even when donors are swapped together. A conditional permutation test has no reference distribution for a tissue effect here. An additive scanner term is collinear with tissue in the two-arm design, and a batch-corrected score would depend on untestable assumptions. **NO-GO for an independent Erlangen lesion increment from these arrays.** Obtaining another table of pathology or cycle phase for these same arrays would leave the scanner problem intact. A cross-instrument bridge using the same biological material or a genuinely crossed independent series would address it. [Conditional enumeration](../../scripts/sohler_conditional_design.py) · [Design counts](../../evidence/peritoneal-followup/sohler-scans.json)

## Finnish run overlap and treatment

Seventeen Finnish case donors have both tissue labels. Eight have one common assay run and one specimen of each tissue in it. Seven pass the physical address gate in both arms. All three assays from the remaining donor contain repeated addresses. The surviving donors contribute exactly 14 WG-6 arrays. The previous finding in eight noncase donors did not test these case pairs. The new census read only assay IDs and address keys in the already cached archive. [Case and control capacity](../../evidence/peritoneal-followup/finnish-design.json) · [Clinical release](https://doi.org/10.6084/m9.figshare.12195975.v1)

| Native run | Case donors with both tissues | Distinct noncase peritoneum donors | Cases with recorded medication |
| --- | ---: | ---: | ---: |
| 9 | 1 | 5 | 1 |
| 10 | 5 | 2 | 4 |
| 12 | 1 | 1 | 1 |
| Total | 7 | 8 | 6 |

All eight noncase donors have stored false, the source database default. Only one case has a reference donor in the same run with the same recorded flag. None of the six recorded treated cases does. Four of the eight reference donors also contributed to the previous organ experiment. Substituting a three-group regression or residualizing treatment cannot create within-run treated controls, verify actual exposure, supply a subsite or certify negative pathology. The seven within-case pairs support a more limited comparison that keeps the donor and recorded run constant. Five pairs share one run. If all donor differences within a run can share a technical sign, only three independent run signs remain. Flipping those three signs together within run gives eight assignments and a smallest possible conventional two-sided p of 2/8 = 0.25. This is a design bound, independent of the outcomes. Only one proliferative and one secretory uterine classification are available among the seven, without a same-operation key. [Aggregate counts](../../evidence/peritoneal-followup/finnish-design.json)

The cohort registry, both investigators’ theses, the collection paper, the released application queries and the Figshare version inventory were checked for a public pathology crosswalk. Patient-level diagnostic confirmation and retained histology aliquots are documented. No public key recovered lesion-negative microscopy and location for the assayed comparison biopsies. All 163 released specimen-linked histology records are uterine. A separate primary pathology study found microscopic disease in visually normal peritoneum in 23 of 151 women with visible disease and four of 62 without visible disease. Those observations explain why specimen-level pathology matters. Transferring those proportions to Finland would require an unsupported assumption about sampling and population. [Finnish source audit](FINNISH_PERITONEAL_PATHOLOGY_2026-11-16.md) · [Khan et al.](https://doi.org/10.1093/humrep/det438)

## Matched biopsies beyond the original search

| Cohort or archive | What the source adds | Decision for the five-gene question |
| --- | --- | --- |
| Ibrahim et al., Charité Berlin | Twenty-three donors supplied histology-confirmed superficial lesions, five of those supplied additional unaffected peritoneum away from lesions, and ten surgical noncases supplied peritoneum. Methods specify lateral pelvic wall sampling and histological exclusion in noncases. Case biopsies span cycle phases, with three-month preoperative hormone exclusion. Two controls had intrauterine devices and exposure was unknown for three. | Strongest new anatomical benchmark. Tissue was fixed and immunostained for smooth-muscle, matrix and TGF-β receptor proteins. No matched public whole-transcriptome deposit was located. Availability of residual tissue or RNA is unknown. [Primary paper](https://doi.org/10.1007/s00404-018-4995-y) |
| Tran et al., Sydney | A 56-woman immunohistochemical study included lesions, normal peritoneum at least 5 cm from lesions in cases, and peritoneum from women without endometriosis. | An explicit distance criterion is feasible and valuable. Macrophage and nerve stains cannot reconstruct the mRNA score. [Methods](https://doi.org/10.1093/humrep/den483) |
| Barcena de Arellano et al., Berlin | Sixty lesion specimens, 60 distant case peritoneal biopsies and ten noncase specimens were studied by immunohistochemistry. | Specimen totals do not establish pair counts. Possible overlap with the other Berlin collection is unverified. No eligible RNA release was found. [Primary paper](https://pubmed.ncbi.nlm.nih.gov/21831994/) |
| Gunther et al., independent discovery catalogue | Its January 2025 GEO/ArrayExpress search yielded 122 included series. We inspected the control field of every row and found eight hits for adjacent tissue or peritoneum. All six peritoneum-labelled accessions appear in our prior report. | The other two hits were ovarian stromal chromatin (GSE197928) and ovarian lesion/adjacent ovary mRNA (GSE86534). Neither adds a pelvic peritoneum baseline. Each accession is adjudicated against its primary paper or official deposit. The catalogue is a discovery index with a date cutoff and a free-text field. [Catalogue and paper](https://doi.org/10.7554/eLife.103825) · [Eight record-level decisions](../../evidence/peritoneal-followup/catalogue.json) |
| Henze et al., 2026 preprint | The single-cell peritoneal reference reuses Tan's adjacent tissue. The bulk analysis reuses Finnish GSE141549. | Another analysis of known material. It cannot add independently recruited noncase peritoneum. [Preprint](https://doi.org/10.64898/2026.07.10.737856) |

Selection was outcome-blind for this task. We searched primary papers and deposits for peritoneal wall, peritoneal biopsy, normal-appearing, normal, unaffected, distant and adjacent tissue with superficial endometriosis, microarray, RNA-seq or spatial assays. We followed the Finnish registry and original investigators' descriptions, checked assay designs in pathology studies and traced the newer integration back to its contributors. The independent catalogue supplements the original 34-accession GDS screen. These searches can miss poorly indexed or newly released archives. The specific evidence that would reopen a full test is a public link from the profiled specimen to histological positivity or negativity and pelvic subsite, plus sufficient eligible donors across both case and noncase comparisons. An operative distance or separate distant site is needed to distinguish a local field from a wider case-peritoneum difference. [Original search and release criteria](PERITONEAL_CONTROL_AUDIT_2026-11-16.md#search-and-reproduction)

## Decision and execution boundaries

| Scientific question | Decision | Next action |
| --- | --- | --- |
| Is there an additional lesion signal above anatomically matched, histologically negative case and noncase peritoneum? | **NO-GO** for the audited public transcriptomes. | Reopen on a specimen-linked pathology and same-site whole-transcriptome release. Charité and Edinburgh document sampling designs that would discriminate the explanations. Remaining aliquots or unpublished RNA have not been established. |
| Can the independently recruited German arrays identify lesion versus case-peritoneum biology? | **NO-GO.** Scanner and acquisition period fully track tissue. | A crossed biological bridge or new cohort would be needed. The original 13 pairs are unsuitable as an independent specificity test. |
| Is the score higher in surgically labelled superficial lesions than in the same Finnish cases’ peritoneal samples? | **GO for a separate, prespecified seven-donor surgical-label comparison.** | Execute the [frozen within-case protocol](FINNISH_CASE_PERITONEUM_PLAN_2026-11-16.md) only under separate outcome authorization. Report the interval, direction and run sensitivity even if the result is inconclusive. |

A positive within-case result would quantify an increment beyond one peritoneal comparator. Location, local reaction, cellular mixture and microscopic disease could still explain it. A reversed effect would challenge the intended direction in these cases. A small estimate with wide uncertainty would leave the remaining lesion effect unresolved. The seven pairs, six prior core–margin pairs and eight prior noncase organ pairs represent different contrasts. Pooling their p-values or subtracting cross-cohort means would lose those distinctions.

## Reproducible checks

Task #55 started from main commit `817df03a5af765987892fbefe761cdf743016e81`. The queue, active issues, branches and PRs were checked against #49 and #53. The remaining open 3p21.31 locus PR covers a different question. Exact feature membership, checked input digests, cohort counts and scanner assignment support are in the [aggregate evidence directory](../../evidence/peritoneal-followup/). The [P0 task](../../research/tasks/peritoneal-followup-qualification.json) caps new source transfers at 40 MB, local RAM at 2 GB and execution at 120 minutes. All analyses below are deterministic.

Use the existing native inputs from the [normal-organ protocol](NORMAL_PERITONEUM_UTERUS_RESULTS_2026-11-16.md#reproduction-and-verification) and its inflate64 1.0.4 dependency. Cache the [BioStudies study JSON](https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-694) as `biostudies.json`, the [SDRF](https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/694/E-MTAB-694/Files/E-MTAB-694.sdrf.txt) under its published name, the [GPL570 feature manifest](https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz) as `GPL570.annot.gz`, and the [discovery workbook](https://cdn.elifesciences.org/articles/103825/elife-103825-supp1-v1.xlsx) as `elife-103825-supp1-v1.xlsx` under the ignored directory below. Use the [Jiang clinical workbook](https://pmc-oa-opendata.s3.amazonaws.com/PMC7060988.1/peerj-08-8730-s001.xlsx) at the indicated private path. EMBL-EBI provides public file access. A specific Erlangen reuse license remains unverified. The Finnish v1 deposit specifies CC BY 4.0. The repository includes newly computed aggregates and code.

The 57 archived CELs total 772,835,154 bytes. Their selected 512-byte prefixes total 29,184 bytes. The header reader requires an exact partial response and checks that those bytes lie in the declared text header, before any intensity cells. An initial single 1,024-byte format probe was also retrieved, from which only the declared header was inspected. The production census uses bytes 0–511 only. An initial transfer timed out and resumed using cached prefixes before the complete ordered digest check. The downloaded GPL570 feature annotation was 8,471,521 compressed bytes. The Finnish case audit reused the existing pinned native archive prefix without another large transfer.

```sh
d=data/private/peritoneal-followup
e=evidence/peritoneal-followup
python scripts/audit_sohler_scans.py "$d" --fetch-headers | diff - "$e/sohler-scans.json"
python scripts/audit_peritoneal_catalogue.py "$d/elife-103825-supp1-v1.xlsx" | diff - "$e/catalogue.json"
python scripts/audit_sohler_clinical.py "$d/E-MTAB-694.sdrf.txt" data/private/jiang-8730-s001.xlsx | diff - "$e/sohler-clinical.json"
python scripts/audit_finnish_case_peritoneum_design.py data/private/normal-peritoneum-uterus | diff - "$e/finnish-design.json"
python -m pytest -q
python scripts/research_queue.py validate
python scripts/verify_catalog.py
python scripts/repo_safety.py
```

After the prefixes have been saved, omit `--fetch-headers` for offline reproduction. The native case-pair lock command ran in 22.6 seconds and verified seven donors and 14 assays. A separately authorized follow-up can run the checked scorer offline. Its primary and independent arithmetic paths were exercised on synthetic records only. [Protocol frozen before outcomes](https://github.com/dovvnloading/Open-Disease/commit/1a1bbaa540b4159069a6dcd16f3cf7db9e31f78f) · [Execution command](../../scripts/finnish_case_peritoneum_score.py)

For that follow-up, the two commands use distinct modes. The first was executed on the native inputs in this audit. The second has been verified on synthetic records and is reserved for authorized outcome access.

```sh
python scripts/finnish_case_peritoneum_score.py lock data/private/normal-peritoneum-uterus
python scripts/finnish_case_peritoneum_score.py execute data/private/normal-peritoneum-uterus
```

The primary physical census was independently checked by parsing the native acquisition field, matching the embedded basename to each archive file and recomputing the donor and scanner-day tables plus all 8,192 permutations. Both implementations obtained 57 physical files, 13 cross-instrument pairs, one permissible assignment and the same ordered SHA-256. The clinical and catalogue tests also check source version, uniqueness and the interpretation of metadata fields. No case expression outcome was computed during qualification.
