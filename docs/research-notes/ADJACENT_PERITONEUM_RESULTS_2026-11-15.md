# Adjacent peritoneum carries the five-gene lesion-versus-uterus signal

**Exploratory analysis of a separately recruited human cohort, November 15, 2026.** In six women with all three tissue sites sampled, the frozen score was 11.55 percentile points higher in peritoneal lesions than in paired uterine tissue. Nearby peritoneum showed an equally prominent shift against the same uterine biopsies (+12.00 points). The direct lesion-minus-adjacent difference was **−0.46 points**, positive in two of six women. Its **95% constant-shift compatibility interval was −6.38 to +6.40 points**. These estimates weaken the interpretation that the whole-tissue panel identifies the lesion core. Shared anatomy, cell composition and a local disease field are compatible with the observations. [Aggregate results](../../evidence/adjacent-peritoneum/results.json) · [Frozen five-gene panel](../../evidence/lesion-transfer/panel.json) · [Local plan](ADJACENT_PERITONEUM_PLAN_2026-11-15.md)

| Same-patient comparison | Pairs | Mean score change, pp | Positive pairs | Leave-one-patient-out range, pp |
| --- | ---: | ---: | ---: | ---: |
| Peritoneal lesion minus nearby peritoneum, primary | 6 | −0.46 | 2 | −1.83 to +0.73 |
| Peritoneal lesion minus uterus | 6 | +11.55 | 5 | +6.12 to +14.87 |
| Nearby peritoneum minus uterus | 6 | +12.00 | 6 | +7.31 to +14.14 |

The prespecified increase toward lesion cores gave a one-sided paired sign-flip p-value of **37/64 = 0.578**. The two-sided p-value was **56/64 = 0.875**. The interval inverts this exact two-sided test over candidate shifts, assuming independent pairs and sign exchangeability of differences under a common additive shift. With six women it admits differences that could be biologically relevant. We specified no equivalence margin. The two comparisons against the same uterus describe the source of the primary contrast and have no separate inferential claim.

## What was actually matched

Tan et al. recruited patients undergoing surgery at the University of Connecticut Health Center. The study involved stage II–IV endometriosis; most participants received similar hormonal treatment. The GEO protocol describes oral contraceptive use. The team removed each visible peritoneal lesion with adjacent tissue and separated the lesion and adjacent regions for dissociation. They document H&E histology for the endometrium. The deposited metadata and study methods do not establish negative histology for every nearby aliquot or a uniform distance from the lesion. The comparator represents the lesion surroundings in someone with diagnosed endometriosis. It cannot establish that disease-free peritoneum would show the same pattern. [Tan et al., Methods and Extended Data Figure 1](https://doi.org/10.1038/s41556-022-00961-5) · [Official GEO series](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179640)

GSE179640 includes both bulk and single-cell sequencing. The bulk arm has six peritoneal lesion libraries and **zero adjacent peritoneal libraries**. We checked the clinical patient numbers in the deposited sample titles alongside tissue, method and location fields. Eight patients had peritoneal lesion single-cell libraries; six had adjacent tissue and all six also had a eutopic biopsy. We therefore used exactly the **19 filtered 10x libraries** from those six people: six lesions, seven adjacent libraries and six uterine biopsies. One patient's adjacent tissue had two libraries. They receive equal weight within that patient in the primary analysis. The processed matrices contain 24,798, 20,873, and 30,348 retained cell barcodes across the three sites, respectively. The replication units are the six women.

We collapsed raw filtered-matrix UMI counts across recovered cell barcodes, ranked all 33,538 gene-expression features within each library, and calculated each panel gene's rank percentile. Exact gene symbol and Ensembl pairs identify all five features. Higher TMEM255A, DIRAS3 and JCHAIN and lower GABRP and SALL1 each contribute equally to the signed score. Summing across the duplicate adjacent libraries before ranking gave −0.61 points for lesion minus adjacent tissue. Keeping just one of those libraries gave −0.51 and −0.40 points. These sensitivities address duplicate weighting. Viable-cell sorting, differential survival during dissociation and regional sampling all affect the recovered composition. [Scoring and checks](../../scripts/adjacent_peritoneum.py)

| Gene | Original predicted lesion direction | Lesion minus adjacent, pp | Lesion minus uterus, pp | Adjacent minus uterus, pp |
| --- | --- | ---: | ---: | ---: |
| TMEM255A | Higher | +1.30 | +7.21 | +5.90 |
| DIRAS3 | Higher | −2.59 | +2.60 | +5.19 |
| GABRP | Lower | −4.17 | −25.81 | −21.64 |
| SALL1 | Lower | +0.58 | −21.75 | −22.33 |
| JCHAIN | Higher | −4.59 | +0.37 | +4.96 |

The two genes measurable in the independent spatial gene tables, GABRP and SALL1, both had substantially lower relative ranks in nearby peritoneum than in uterus. The earlier epithelial result compares separate anatomical tissues. This same-site control therefore changes how far that result can support a lesion-specific expression program. Both analyses remain relative-expression measurements among treated or surgically selected participants. [Five-pair spatial compartment analysis](SPATIAL_PANEL_RESULTS_2026-11-15.md)

## Reproduction

The primary design was fixed locally before inspecting expression in this cohort. It follows knowledge of the bulk and spatial results. The interval was added after the estimate and p-value to quantify precision. The program verifies the SHA-256 digests of the pinned panel, native GEO family SOFT and file list. It derives 19 exact file URLs and byte lengths from these sources, checks patient pairing and per-file alignment, and verifies the combined source digest. H5 downloads total 315,533,135 bytes. The digest of the 19 hexadecimal SHA-256 strings concatenated in filename order is `2a2d62c1531355122be2de5bcfe1f748fef86a6e59e6e9e566f10880b0bc1bcb`. Sparse CSC shape, count integrity and identity of gene order are checked for every biopsy. No barcodes, patient-level scores or expression matrices are redistributed. The exact tests use no random seed. Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0 and h5py 3.15.1 were used.

```bash
python -m pip install h5py==3.15.1
export PYTHONPATH=src:scripts OPENBLAS_NUM_THREADS=1
mkdir -p results/adjacent-peritoneum
python scripts/adjacent_peritoneum.py --fetch > results/adjacent-peritoneum/results.json
cmp results/adjacent-peritoneum/results.json evidence/adjacent-peritoneum/results.json
python -m pytest -q tests/test_adjacent_peritoneum.py
```

Omit `--fetch` to replay with cached inputs. [Native GEO source manifest](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179640/suppl/filelist.txt)
