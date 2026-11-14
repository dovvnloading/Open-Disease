# FOSB annotation in the combined Finnish workbook

**The combined workbook has exactly one FOSB annotation row.** On `Sheet1`, Excel row **21,565** pairs `Gene_symbol = FOSB` with `Probe_Id = ILMN_1751607`, the HumanHT-12 v4 identifier. Across all 28,247 data rows there is no `ILMN_1730706` label, the identifier used for FOSB by HumanWG-6 v2. The study's pooled table therefore identifies a single reporting row. Its two annotation columns do not document how the authors populated the older-generation measurements in that row. [Combined GEO workbook](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141549/suppl/GSE141549_batchCorrectednormalizedArrayscombined.xlsx) · [Source audit and manufacturer probe comparison](GSE141549_ELIGIBILITY_2026-11-13.md)

The author methods describe remapping probe sequences, aligning gene symbols across chip generations and adjusting batches. The publisher's versioned data release lists an SQL archive, a qPCR sample-count workbook and three plots. The study's GEO supplement list includes per-platform and combined data. We found no separately published two-platform probe-pair table or worksheet in the combined workbook documenting this particular pairing. A row labeled with the newer chip's identifier is insufficient evidence that the older chip physically measured the same 50-base oligonucleotide. The platform manifests already establish that their FOSB oligonucleotides differ. [Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC7455745/) · [Versioned author release file list](https://api.figshare.com/v2/articles/12195975/versions/1) · [GEO supplementary file index](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141549/suppl/)

## Independent check of the pooled annotation

The GEO combined workbook was downloaded locally (153,632,401 bytes, SHA-256 `9919fc5d6d6c09dbdae5a815da798b8250edf17a781fdb86dc25c2f9c9b21d50`). The audit decoded the shared-string dictionary and only the A and B cells in each worksheet row. The sheet declares `A1:OT28248` and all 28,248 A/B pairs, including the header, were present. The worksheet XML was decompressed as bytes; no sample measurement cells were converted to numbers, compared or reported.

Repeat with Python 3.11+ on the pinned workbook at an ignored local path. This streaming implementation checks that the annotated columns are complete and aligned by worksheet row.

```python
from pathlib import Path
import hashlib
import re
from xml.etree import ElementTree
from zipfile import ZipFile

book = Path("data/private/gse141549/combined.xlsx")
with book.open("rb") as source:
    assert hashlib.file_digest(source, "sha256").hexdigest() == (
        "9919fc5d6d6c09dbdae5a815da798b8250edf17a781fdb86dc25c2f9c9b21d50"
    )
pair = re.compile(
    rb'<row\s+[^>]+>\s*<c\s+r="A(\d+)"([^>]*)><v>(\d+)</v></c>'
    rb'\s*<c\s+r="B(\d+)"([^>]*)><v>(\d+)</v></c>'
)
with ZipFile(book) as archive:
    dictionary = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
    strings = ["".join(element.itertext()) for element in dictionary]
    carry = b""
    previous = 0
    matches = []
    with archive.open("xl/worksheets/sheet1.xml") as stream:
        for chunk in iter(lambda: stream.read(800_000), b""):
            buf = carry + chunk
            end = buf.rfind(b"</row>")
            if end < 0:
                carry = buf
                continue
            complete, carry = buf[:end + 6], buf[end + 6:]
            for a_row, a_type, a_index, b_row, b_type, b_index in pair.findall(complete):
                row = int(a_row)
                assert a_row == b_row and row == previous + 1
                assert b't="s"' in a_type and b't="s"' in b_type
                previous = row
                annotation = (strings[int(a_index)], strings[int(b_index)])
                if row == 1:
                    assert annotation == ("Gene_symbol", "Probe_Id")
                if annotation[0] == "FOSB" or annotation[1] in {"ILMN_1730706", "ILMN_1751607"}:
                    matches.append((row, *annotation))
    assert previous == 28_248
    print(matches)
```

## Separate HumanWG-6 v2 source

GEO directly serves the normalized WG-6 workbook (216,488,862 bytes, HTTP 200). Bounded byte ranges sufficed to check its first worksheet header and complete string dictionary. It has two annotation headers, `ID_REF` and `ILMN_Gene`, and 336 GEO sample columns including 62 primary case and 38 primary control endometrium labels. The string dictionary includes `FOSB` and `ILMN_1730706` and omits `ILMN_1751607`. We did not scan that entire worksheet or derive an expression result. The WG-6 product provides a chip-specific source independently of a join to the native patient tables. The prior audit's donor, treatment and phase eligibility gates still apply before fitting a model. [WG-6 normalized GEO workbook](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141549/suppl/GSE141549_Illumina_HumanWG-6_v2.0_expression_beadchip_normalized.xlsx)
