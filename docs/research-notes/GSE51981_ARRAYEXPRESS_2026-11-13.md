# GSE51981 ArrayExpress source check

**November 13, 2026 (US Eastern).** The complete EMBL-EBI sample and data relationship file (SDRF) contains 148 rows and 42 columns. Every `Assay Name` matches a distinct GSE51981 accession. Every `Source Name` equals that same assay accession followed by a space and the digit 1. The SDRF therefore supplies no independently documented participant key. Its columns have no field for medication or recent hormonal treatment. This resolves the file-access uncertainty in the [cohort audit](GSE51981_COHORT_AUDIT_2026-11-13.md). The donor-level model in [issue #18](https://github.com/dovvnloading/Open-Disease/issues/18) remains blocked. [EMBL-EBI study record](https://www.ebi.ac.uk/biostudies/api/v1/studies/E-GEOD-51981)

## Cross-check against GEO

Using the assay accession for the join, the sample titles, descriptions, source names and organism labels matched the corresponding GEO fields in all 148 records. The following counts use only characteristics explicitly present in GEO. Comparisons of characteristic values ignore letter case. [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51981)

| Characteristic | GEO records with field | Matching SDRF values |
| --- | ---: | ---: |
| Endometriosis status | 148 | 148 |
| Endometrial tissue and phase label | 148 | 148 |
| Endometriosis severity | 77 | 77 |
| Presence or absence of uterine/pelvic pathology | 71 | 71 |

ArrayExpress imported the GEO submission. This agreement verifies the transfer of the recorded labels. The investigation description file (IDF) reports 148 women at the study level. Its experimental factors cover disease status, severity, pelvic pathology and tissue. A mention of treatment in the extraction protocol refers to DNase treatment of an RNA preparation. The paper reports recent hormonal therapy as a cohort exclusion. Neither retrieved file supplies a participant-level exposure observation. The original clinical workbook could not be obtained from the publisher in this environment. The publication’s [BioStudies record](https://www.ebi.ac.uk/biostudies/api/v1/studies/S-EPMC4239429) lists no deposited file. Its contents remain unverified. [Original study, specimens and processing](https://pmc.ncbi.nlm.nih.gov/articles/PMC4239429/)

## Inputs and reproduction

The BioStudies archive returned the full SDRF and IDF on November 14, 2026 UTC. A second retrieval of the complete SDRF via the legacy ArrayExpress route returned the same digest. Individual source rows and clinical attributes stayed in ignored local storage. No expression measurements were inspected. [Study storage metadata](https://www.ebi.ac.uk/biostudies/api/v1/studies/E-GEOD-51981/info) · [Legacy SDRF](https://ftp.ebi.ac.uk/pub/databases/arrayexpress/data/experiment/GEOD/E-GEOD-51981/E-GEOD-51981.sdrf.txt)

| Input | Exact bytes | SHA-256 |
| --- | ---: | --- |
| [SDRF](https://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/981/E-GEOD-51981/Files/E-GEOD-51981.sdrf.txt) | 122,872 | `29cba159a95ed825643f01459a5992e9a687d19dd9c3c2807c8af2b18c844c79` |
| [IDF](https://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/981/E-GEOD-51981/Files/E-GEOD-51981.idf.txt) | 5,983 | `a55417e11c22be98b8364f05ae86eaf85417fef42ab0687da88e7ee3d8078dfa` |

Run from the repository root using CPython 3.11 or later. The [cohort audit](GSE51981_COHORT_AUDIT_2026-11-13.md) gives the opt-in acquisition command and digest for `gsm-quick.soft`. In a fresh working directory:

```bash
mkdir -p data/private/gse51981-audit
base=https://ftp.ebi.ac.uk/biostudies/fire/E-GEOD-/981/E-GEOD-51981/Files
for kind in sdrf idf; do
  name="E-GEOD-51981.${kind}.txt"
  curl --fail --silent --show-error --max-time 60 --max-filesize 150000 \
    --output "data/private/gse51981-audit/$name" "$base/$name"
done
(cd data/private/gse51981-audit && sha256sum E-GEOD-51981.*.txt)
```

```python
import csv
from pathlib import Path
import sys

sys.path.insert(0, "scripts")
from reconcile_gse51981_sources import CHAR, read_quick_view

directory = Path("data/private/gse51981-audit")
geo = read_quick_view(directory / "gsm-quick.soft")
with (directory / "E-GEOD-51981.sdrf.txt").open(encoding="utf-8", newline="") as stream:
    table = csv.reader(stream, delimiter="\t")
    columns = next(table)
    records = list(table)

assert len(columns) == 42 and len(records) == len(geo) == 148
assert all(len(row) == len(columns) for row in records)
assert {row[columns.index("Assay Name")] for row in records} == set(geo)
fields = {
    "endometriosis/no endometriosis": "Characteristics [endometriosis/no endometriosis]",
    "tissue": "Characteristics [organism part]",
    "endometriosis severity": "Characteristics [endometriosis severity]",
    "presence or absence of uterine/pelvic pathology": "Characteristics [presence or absence of uterine/pelvic pathology]",
}
counts = dict.fromkeys(fields, 0)
for row in records:
    accession = row[columns.index("Assay Name")]
    assert row[columns.index("Source Name")] == accession + " 1"
    attributes = dict(item.split(": ", 1) for item in geo[accession][CHAR])
    for key, column in fields.items():
        if key in attributes:
            assert attributes[key].casefold() == row[columns.index(column)].casefold()
            counts[key] += 1
print(counts)
```

The checks used CPython 3.12.14, deterministic parsing, and no random seed.
