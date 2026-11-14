"""Synthetic regressions for source selection and two-level multiplicity."""
import json
from pathlib import Path
import sys

from openpyxl import Workbook
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import analyze_melbourne_genes as melbourne


def row(probe, symbol, coef=-0.2, p=.001):
    return (probe, coef, 1., 2., p, 1., 2., 5., symbol, "1")


def as_rows(*contents):
    return iter([(melbourne.TITLE,), melbourne.HEADER, *contents])


def test_exact_symbol_mapping_and_two_level_correction():
    probes = melbourne.candidate_probes(as_rows(
        row("ILMN_4", " GLTP "), row("ILMN_5", "GLTP", p=.04),
        row("ILMN_7", "GLTPP1"), row("ILMN_10", "ABO", p=.04)), expected_count=4)
    genes = {v["symbol"]: v for v in melbourne.summarize(probes)}
    assert [x["probe_id"] for x in probes] == ["ILMN_4", "ILMN_5", "ILMN_10"]
    assert genes["GLTP"]["gene_p"] == .002 and genes["GLTP"]["holm_p"] == .01
    assert genes["GLTP"]["direction_agreement"] == "source_orientation_unresolved"
    assert genes["ABO"]["holm_p"] == .16
    assert genes["F11"]["holm_input_p"] == 1 and genes["F11"]["holm_p"] is None
    assert genes["F11"]["two_sided_result"] == "untested"


@pytest.mark.parametrize("records,error", [
    ([row("ILMN_5", "other"), row("ILMN_5", "other")], "repeated"),
    ([row("unmapped", "other")], "Invalid or repeated"),
    ([row("ILMN_1", "GLTP", p=float("nan"))], "Invalid coefficient"),
    ([row("ILMN_1", "KIF3C", coef=float("inf"))], "Invalid coefficient"),
    ([row("ILMN_1", "ABO", p=1.1)], "Invalid coefficient"),
    ([row("ILMN_1", "F11")], "universe"),
])
def test_malformed_rows_fail_closed(records, error):
    with pytest.raises(ValueError, match=error):
        melbourne.candidate_probes(as_rows(*records), expected_count=2)


def test_changed_nomination_or_sheet_rejected():
    with pytest.raises(ValueError, match="directions"):
        melbourne.check_nominations({"candidates": [{"id": "28", "direction": "lower", "family": "gene"}]})
    with pytest.raises(ValueError, match="title"):
        melbourne.candidate_probes(iter([("TableS21",), melbourne.HEADER]))


def test_end_to_end_streamed_workbook_and_pinned_gates(tmp_path, monkeypatch):
    path, lock, out = tmp_path / "synthetic.xlsx", tmp_path / "lock.json", tmp_path / "out"
    book = Workbook(write_only=True)
    sheet = book.create_sheet("TableS20")
    # Table S20 has formatting out to P. Its source dimensions must survive serialization.
    sheet.append([melbourne.TITLE])
    sheet.append(melbourne.HEADER)
    for i in range(1, 15263):
        values = list(row(f"ILMN_{i}", "GLTP" if i == 2 else "OTHER"))
        if i == 15262:
            values += [None] * 5 + ["synthetic"]
        sheet.append(values)
    book.save(path)
    # openpyxl write-only creates no dimensions header; add the pinned published dimension.
    from zipfile import ZipFile
    patched = tmp_path / "with-dimensions.xlsx"
    with ZipFile(path) as source, ZipFile(patched, "w") as dest:
        for info in source.infolist():
            data = source.read(info)
            if info.filename == "xl/worksheets/sheet1.xml":
                data = data.replace(b"<sheetViews>", b'<dimension ref="A1:P15264"/><sheetViews>', 1)
            dest.writestr(info, data)
    lock.write_text(json.dumps({"candidates": [
        {"id": str(gene), "family": "gene", "direction": direction}
        for gene, direction in melbourne.TARGETS.values()]}))
    with pytest.raises(ValueError, match="checksum"):
        melbourne.run(patched, lock, out)
    monkeypatch.setattr(melbourne, "SOURCE_HASH", melbourne.sha256_file(patched))
    genes, probes = melbourne.run(patched, lock, out)
    assert len(genes) == 5 and len(probes) == 1
    assert (out / "melbourne_probes.tsv").read_text().splitlines()[1].split("\t")[2] == "ILMN_2"
    assert len((out / "melbourne_genes.tsv").read_text().splitlines()) == 6
