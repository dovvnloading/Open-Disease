"""Fail-closed clinical workbook projection and file integrity on invented data."""
from pathlib import Path
from types import SimpleNamespace
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import audit_gse266265_epithelial as atlas


def test_rds_requires_both_length_and_checksum(tmp_path, monkeypatch):
    rds = tmp_path / "toy.rds"
    rds.write_bytes(b"toy")
    with pytest.raises(ValueError, match="SHA-256"):
        atlas.verify(rds)
    monkeypatch.setattr(atlas, "SIZE", 3)
    with pytest.raises(ValueError, match="SHA-256"):
        atlas.verify(rds)
    monkeypatch.setattr(atlas, "DIGEST", atlas.sha256_file(rds))
    atlas.verify(rds)


def test_projection_reads_only_design_fields_and_validates_header(tmp_path, monkeypatch):
    bookpath = tmp_path / "toy"
    bookpath.write_bytes(b"synthetic workbook stand-in")
    monkeypatch.setattr(atlas, "WORKBOOK_DIGEST", atlas.sha256_file(bookpath))
    columns = [*atlas.SOURCE_COLUMNS, "unused"]
    sample = ["simulated", "non-ENDO", "proliferative", "simA", "simS", True, "ignored"]
    stub = SimpleNamespace(active=SimpleNamespace(values=[columns, sample, [None] * 7]), close=lambda: None)
    monkeypatch.setitem(sys.modules, "openpyxl", SimpleNamespace(load_workbook=lambda *a, **kw: stub))
    out = tmp_path / "study-design.tsv"
    atlas.write_sample_map(bookpath, out)
    text = out.read_text()
    assert len(text.splitlines()) == 2
    assert "simulated\tnon-ENDO\tproliferative\tsimA\tsimS\tTrue" in text
    assert "unused" not in text and "ignored" not in text
    columns[0] = "unknown"
    with pytest.raises(ValueError, match="headers"):
        atlas.write_sample_map(bookpath, out)
