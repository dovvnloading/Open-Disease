"""The frozen single-probe contrast, exercised on invented arrays and intensities."""

import csv
import gzip
import hashlib
from itertools import zip_longest
import subprocess
import sys
from pathlib import Path

import pytest

from scripts import analyze_gse25628_fosb as fosb


def matrix_table():
    cases = [("patient", "Normal (eutopic)", n) for n in (3, 4, 5, 5, 6, 6, 7, 8)]
    controls = [("healthy donor", "Normal (control)", n) for n in (2, 3, 3, 3, 3, 4)]
    lesions = [("patient", "Pathological (ectopic)", n) for n in (1e8, -1e8) * 4]
    # Interleave the three groups to ensure inclusion follows fields, not positions.
    samples = [s for trio in zip_longest(lesions, cases, controls) for s in trio if s]
    ids = [f"GSM9900{i:03d}" for i in range(len(samples))]
    return [
        ["!Series_geo_accession", "GSE25628"], ["!Series_platform_id", "GPL571"],
        ["!Sample_geo_accession", *ids], ["!Sample_platform_id", *(["GPL571"] * 22)],
        ["!Sample_data_processing", *(["RMA + quantile normalization"] * 22)],
        ["!Sample_characteristics_ch1", *["donor: " + s[0] for s in samples]],
        ["!Sample_characteristics_ch1", *(["tissue: endometrium"] * 22)],
        ["!Sample_characteristics_ch1", *["disease state: " + s[1] for s in samples]],
        ["!series_matrix_table_begin"], ["ID_REF", *ids],
        ["SYN_OTHER", *([1234.5] * 22)], [fosb.PROBE, *[s[2] for s in samples]],
        ["!series_matrix_table_end"],
    ]


def write_table(path, rows):
    with gzip.open(path, "wt", newline="", encoding="utf-8") as stream:
        csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)
    return path


def annotation_table():
    return [
        ["!platform_table_begin"], ["ID", "Gene symbol", "Gene ID"],
        [fosb.PROBE, "FOSB", "2354"], ["!platform_table_end"],
    ]


def test_frozen_single_probe_contrast(tmp_path, monkeypatch):
    matrix = write_table(tmp_path / "matrix.gz", matrix_table())
    annotation = write_table(tmp_path / "annotation.gz", annotation_table())
    monkeypatch.setattr(fosb, "MATRIX_HASH", hashlib.sha256(matrix.read_bytes()).hexdigest())
    monkeypatch.setattr(fosb, "ANNOTATION_HASH", hashlib.sha256(annotation.read_bytes()).hexdigest())
    report = fosb.analyze(matrix, annotation)

    assert report["counts"] == {"case": 8, "control": 6, "lesion": 8}
    assert (report["features_in_matrix"], report["features_tested"]) == (2, 1)
    # Base R t.test(c(3,4,5,5,6,6,7,8), c(2,3,3,3,3,4)), two-sided/unpaired.
    assert report["difference_log2"] == pytest.approx(2.5)
    assert report["ratio_geometric_means"] == pytest.approx(2**2.5)
    assert report["ci95_log2"] == pytest.approx([1.104567433941, 3.895432566059])
    assert report["welch_df"] == pytest.approx(9.625153964787)
    assert report["two_sided_p"] == pytest.approx(0.002662723638)
    assert report["leave_one_out_range_log2"] == pytest.approx([15/7, 20/7])
    assert report["decision"] == "supports_expected_direction"


@pytest.mark.parametrize(("case", "decision"), [
    ([0, 1, 0, 1, 0, 1, 0, 1], "supports_opposite_direction"),
    ([3, 4, 5, 5, 6, 6, 7, 8], "inconclusive"),
])
def test_decision_uses_entire_interval(case, decision):
    assert fosb.summarize(case, [4, 5, 5, 5, 5, 6])["decision"] == decision


@pytest.mark.parametrize(("key", "before", "after", "error"), [
    ("ID_REF", "GSM9900000", "GSM9999999", "sample order"),
    ("!Sample_geo_accession", "GSM9900000", "GSM9900001", "distinct"),
    ("!Sample_platform_id", "GPL571", "GPL570", "platform"),
    ("!Sample_data_processing", "RMA + quantile normalization", "unverified", "processing"),
    ("!Sample_characteristics_ch1", "donor: healthy donor", "donor: patient", "group"),
    ("!Sample_characteristics_ch1", "disease state: Normal (eutopic)",
     "disease state: Pathological (ectopic)", "count"),
    (fosb.PROBE, 3, "NA", "intensity"),
    (fosb.PROBE, 3, "inf", "intensity"),
])
def test_invalid_matrix_cell_rejected(tmp_path, key, before, after, error):
    rows = matrix_table()
    row = next(r for r in rows if r[0] == key and before in r[1:])
    row[row.index(before, 1)] = after
    with pytest.raises(ValueError, match=error):
        fosb.read_target(write_table(tmp_path / "invalid.gz", rows))


@pytest.mark.parametrize("bad_table", ["duplicate_target", "missing_target", "no_end"])
def test_target_must_exist_once_in_complete_matrix(tmp_path, bad_table):
    rows = matrix_table()
    if bad_table == "duplicate_target":
        rows.insert(-1, rows[-2].copy())
    elif bad_table == "missing_target":
        rows.pop(-2)
    else:
        rows.pop()
    with pytest.raises(ValueError, match="duplicate probe|probe is absent|end of matrix"):
        fosb.read_target(write_table(tmp_path / "bad.gz", rows))


def test_ambiguous_probe_assignment_is_rejected(tmp_path):
    rows = annotation_table()
    rows.insert(-1, ["SYN_OTHER", "FOSB", "2354"])
    with pytest.raises(ValueError, match="probe mapping"):
        fosb.check_annotation(write_table(tmp_path / "platform.gz", rows))


def test_cli_refuses_unpinned_inputs(tmp_path):
    matrix = write_table(tmp_path / "matrix.gz", matrix_table())
    annotation = write_table(tmp_path / "annotation.gz", annotation_table())
    result = subprocess.run(
        [sys.executable, str(Path(fosb.__file__)),
         "--matrix", str(matrix), "--annotation", str(annotation)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0 and not result.stdout
    assert "Input SHA-256 differs" in result.stderr
