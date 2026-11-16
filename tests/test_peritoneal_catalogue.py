"""Synthetic controls distinguish catalogue mentions from eligible tissue groups."""

from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_peritoneal_catalogue import HEADERS, PRIOR_INVENTORY, audit, screen


def row(accession, tissue, control):
    return (accession, None, None, None, None, None, None, tissue, control, None, None)


def test_screen_uses_the_control_field_across_all_rows():
    rows = [
        HEADERS,
        row("Series GSE1", "peritoneal lesion", "uterus"),
        row("Series GSE2", "ovarian lesion", "Adjacent ovary"),
        row("E-MTAB-3", None, "ADJACENT PERITONEUM"),
        row("Series GSE4", None, None),
    ]
    assert screen(rows) == (4, 2, {"GSE2": ["adjacen"], "E-MTAB-3": ["adjacen", "periton"]})


def test_duplicate_study_and_changed_header_are_rejected():
    with pytest.raises(ValueError, match="duplicated"):
        screen([HEADERS, row("Series GSE1", None, None), row("GSE1", None, None)])
    with pytest.raises(ValueError, match="headers"):
        screen([tuple(reversed(HEADERS))])


def test_unversioned_file_rejected_before_workbook_parse(tmp_path):
    path = tmp_path / "changed.xlsx"
    path.write_bytes(b"not the published file")
    with pytest.raises(ValueError, match="checksum"):
        audit(path, PRIOR_INVENTORY)
