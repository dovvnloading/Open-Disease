# Proliferative glandular donor coverage in the Human Endometrial Cell Atlas

**Literature and metadata assessment, November 14, 2026.** Two contributing studies in the Human Endometrial Cell Atlas (HECA) contain both endometriosis cases and controls with at least 20 cells from the proliferative functionalis glands in an untreated superficial biopsy. Under the atlas phase assignments, those studies contribute **seven cases and six controls** and **10,044 eligible gland cells**. The analysis unit is the donor. Only two independent Oxford pairs shared a capture library with any gland cells from both donors. Requiring 20 cells per donor within a library leaves one pair. Huang had no eligible mixed library. The Swiss capture-conditioned test cannot reach its 0.05 threshold here. [Primary study](https://doi.org/10.1038/s41588-024-01873-w) · [Workbook census](../../evidence/heca-glandular-2026-11-14/metadata_census.json)

## Anatomical definition

The Swiss study's retained proliferative epithelial population was labelled `glandular`. In HECA, the authors describe the cell state named `Glandular` as mid-secretory and `preGlandular` as early secretory. They anatomically assign two other states, `SOX9 functionalis I (CDH2+)` and `SOX9 functionalis II`, to proliferative functionalis glands. This census sums those two Table 3 columns within each donor. It excludes basalis, luminal and cycling states with distinct or mixed anatomical localization. The common anatomical compartment supports a prospective comparison across the two studies. Equivalence of the Swiss and HECA cell-state definitions remains unverified. [Authors’ supplementary notes 1 and 3](https://api.repository.cam.ac.uk/server/api/core/bitstreams/375f9e95-64ec-4b79-8139-20bf3b8ec10d/content) · [Swiss plan](SWISS_GLANDULAR_PLAN_2026-11-14.md)

## Which donors survive each restriction?

| Published source and eligibility requirement | Cases | Controls |
|---|---:|---:|
| Unique donors across single-cell and single-nucleus assays | 72 | 49 |
| Donors in the single-cell count table | 31 | 32 |
| Atlas-assigned proliferative phase | 11 | 14 |
| Superficial endometrium and no hormonal treatment for ≥3 months | 11 | 12 |
| At least 20 cells across both eligible gland states per donor | 9 | 10 |
| Component study has at least one eligible donor from both groups | **7** | **6** |

The clinical table records 121 distinct people. Five supplied both cells and nuclei. The 63 single-cell donor rows contain 313,527 cells summed over the 36 published states. The 25 proliferative donors agree with the biopsy count in Figure 1f. The 58 donors who supplied nuclei alone have no cell-state count rows in Table 3. [Study, Figure 1](https://doi.org/10.1038/s41588-024-01873-w) · [Supplementary Tables 1–3](https://pmc-oa-opendata.s3.amazonaws.com/PMC11387200.1/41588_2024_1873_MOESM3_ESM.xlsx)

| Contributing single-cell study | Eligible cases | Eligible controls |
|---|---:|---:|
| Marečková (Oxford) | 4 | 3 |
| Huang | 3 | 3 |
| Fonseca | 2 | 0 |
| Garcia-Alonso | 0 | 1 |
| Lai | 0 | 0 |
| Tan | 0 | 1 |
| Wang | 0 | 2 |

Within the two mixed studies, the workbook records surgical ascertainment for 12 donors and ascertainment through medical history for one Oxford control. None of these 13 has a recorded additional endometrial pathology in Table 1. The authors used the specimen's hormone and phase annotations as separate fields. A Table 1 footnote explains that the Huang phases were based on raw-file metadata amid inconsistencies between source records. The atlas paper also documents transcriptomic phase assignment when histology was unavailable. [Table 1, footnotes b, g–i](https://pmc-oa-opendata.s3.amazonaws.com/PMC11387200.1/41588_2024_1873_MOESM3_ESM.xlsx) · [Study, Methods](https://doi.org/10.1038/s41588-024-01873-w)

## Supplementary-table joins and provenance

Table 1 labels clinical donors, Table 3 prefixes each donor's count-row label with its contributing study, and Table 2 labels libraries. The analysis validated all 63 clinical-to-cell-count matches. Study-name forms were resolved explicitly. The Fonseca label already contains its study prefix. One Lai key carries the published phase-footnote marker; the code checks the specific footnote before reconciling that key. Disease stage zero denotes a control. Nonzero stages and one explicitly documented unstaged case agree with the pathology field. Only the workbook's literal no-treatment code qualifies. The group counts for all filters and all seven components are available in the census. [Analysis code](../../scripts/audit_heca_glandular.py) · [Aggregate results](../../evidence/heca-glandular-2026-11-14/metadata_census.json)

The 71 Table 2 QC rows cover 22 single-cell libraries, 41 nucleus libraries and eight spatial libraries. Table 2 contains no donor field or shared donor identifier connecting the physical libraries to Table 3. The source manuscript describes donors multiplexed in some libraries. Counting donors in mixed component studies therefore precedes the separate assessment of shared captures. The study lists six reused single-cell sources plus its newly generated anchor. That list includes none of GSE141549 (Finland), GSE25628 (Rome), or GSE266265 (Switzerland). The comparison established release-level separation. Individual donor independence across releases was not verified. [Methods and accessions](https://doi.org/10.1038/s41588-024-01873-w)

The workbook is 400,034 bytes, SHA-256 `f70030404af3f6070c10f22a064e20188b0fb3251c1042605563c0aea00d2ed3`. The reproducible runner reads the three named donor, library-QC and cell-count sheets from this pinned download. Run `python scripts/audit_heca_glandular.py data/private/heca/41588_2024_1873_MOESM3_ESM.xlsx` for aggregate JSON. No simulation or random seed enters the census. [Source workbook](https://pmc-oa-opendata.s3.amazonaws.com/PMC11387200.1/41588_2024_1873_MOESM3_ESM.xlsx) · [Issue 40](https://github.com/dovvnloading/Open-Disease/issues/40)

## What survives a comparison within capture libraries?

**The public annotations resolve the missing library link.** We decoded the dataset, donor, library and cell-state fields from the authors' 16,012,377,883-byte whole-cell atlas. These four annotation columns reproduced **all 2,268 donor-by-state counts** in Table 3, including zeros, across 313,527 cells. In Oxford the cell annotations reproduced all 35 donor–library links in the original sample-to-library manifest. The cell-state counts then established which eligible donors actually contributed gland cells to the same library. [Authors’ public dataset](https://www.reproductivecellatlas.org/endometrium_reference.html) · [Oxford sample-to-library manifest](https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/039/E-MTAB-14039/Files/E-MTAB-14039.sdrf.txt) · [Verified cell and library counts](../../evidence/heca-glandular-2026-11-14/cell-library-support.json)

Each donor in the table below first met the 20-cell **total per donor** criterion. The two additional cutoffs describe how many of those cells came from each individual library. The one-cell calculation preserves every observed pair. The 20-cell calculation tests how much overlap survives a stronger within-library measurement requirement. These are design diagnostics.

| Component | Minimum gland cells per donor within a library | Libraries with eligible cases and controls | Distinct cases | Distinct controls | Gland cells in these libraries |
| --- | ---: | ---: | ---: | ---: | ---: |
| Oxford | 1 | 3 | 2 | 2 | 416 |
| Oxford | 20 | 1 | 1 | 1 | 370 |
| Huang | 1 | 0 | 0 | 0 | 0 |
| Huang | 20 | 0 | 0 | 0 | 0 |

The seven Oxford candidates contributed 4,858 qualifying gland cells across eight libraries. Three mixed libraries represent two distinct pairs, because one pair was captured twice. Assigning each donor one disease label and preserving one case per mixed library leaves only four assignments. Complement symmetry gives a minimum attainable two-sided absolute-contrast p-value of 0.50 even under the permissive one-cell criterion. Conditioning additionally on the case counts of all eight libraries leaves two assignments across the seven donors. At the 20-cell within-library threshold only one pair remains. The authors document multiplexing donor suspensions before capture. [Oxford study, donor demultiplexing](https://doi.org/10.1038/s41588-024-01873-w) · [Manifest checks and assignment space](../../evidence/heca-glandular-2026-11-14/design-manifests.json)

The six Huang candidates contribute 5,186 gland cells in six separate libraries. The original GEO records agreed with all ten HECA Huang diagnoses and nine of ten coarse cycle-phase assignments. One control annotated as proliferative in HECA is annotated as secretory in its GEO record and linked ENA BioSample. The GEO and ENA entries reflect the same submitter. Requiring phase agreement leaves three cases and two controls. Both the source authors and GEO distinguish samples from unique participants. These records offer no shared-library case–control contrast. [Huang study and sample-level processing](https://doi.org/10.1093/cei/uxad029) · [GEO release](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214411) · [Source-based design audit](../../scripts/audit_heca_design_manifests.py)

**Feasibility decision.** A capture-conditioned external expression test at the Swiss study's 0.05 threshold cannot reach that threshold in this atlas under either cell cutoff. A comparison across separately processed libraries would require a different technical assumption and a separately specified analysis. This audit did not calculate disease-associated expression estimates. It establishes the independent units and attainable resolution before any further use of these records.

## Reproduce the capture and cell-state checks

The original workbook, Oxford SDRF and Huang SOFT file have pinned SHA-256 values and source URLs in the [workbook census](../../evidence/heca-glandular-2026-11-14/metadata_census.json) and [design-manifest results](../../evidence/heca-glandular-2026-11-14/design-manifests.json). The [HDF5 range receipt](../../evidence/heca-glandular-2026-11-14/annotation-ranges.json) identifies the source object, its ETag and the SHA-256 of each required range. A compact replay uses **1,802,240 bytes** in 15 byte ranges. Only four named observation fields are decoded. The program rejects an unexpected server version, full-file response or range digest. A separately available local copy can be used with `--h5ad`.

With the project Python environment, openpyxl and h5py installed, place the three source metadata files in an ignored directory. Run from the repository root. The last command explicitly opts into a small source transfer. Once the range cache exists, omit `--fetch` for an offline rerun.

```bash
d=data/private/heca
w="$d/41588_2024_1873_MOESM3_ESM.xlsx"
s="$d/E-MTAB-14039.sdrf.txt"
python scripts/audit_heca_glandular.py "$w"
python scripts/verify_heca_glandular.py "$w" --sdrf "$s"
python scripts/audit_heca_design_manifests.py "$w" "$s" "$d/GSE214411_family.soft.gz"
python scripts/audit_heca_cell_libraries.py "$w" "$s" --parts "$d/annotation-parts" --fetch
```

The two census implementations join donors and cell counts by separate methods and agree on every component. Independent reconstruction from the original annotation columns also reproduced the four library coverage rows above. The workbook pass took 0.34 seconds with peak RSS 37,080 KiB. The annotation check used Python 3.12.14, h5py 3.16.0 and openpyxl 3.1.5. Calendar-stage discordance and cohort-specific recruitment limit biological comparability even before an expression comparison. [Cell annotation code](../../scripts/audit_heca_cell_libraries.py) · [Independent workbook and state-count verification](../../scripts/verify_heca_glandular.py)
