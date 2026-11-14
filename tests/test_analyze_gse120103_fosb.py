"""Synthetic checks for the frozen four-cell percentile analysis."""
import csv
import gzip
from hashlib import sha256
from math import sqrt

import numpy as np
import pytest
from scipy.stats import t

from scripts import analyze_gse120103_fosb as fosb


# The interleaved order deliberately differs from the submitter's block order.
CELLS = [
    ("infertile_case", "Endometrium  from Stage IV Ovarian Endometriosis of Infertile women",
     "Fertile Secretory with Stage IV Ovarian Endometriosis (Group 2B)", 2),
    ("fertile_control", "Disease free Endometrium of fertile women",
     "Fertile Secretory without Endometriosis (Group 1A)", -1),
    ("fertile_case", "Endometrium  from Stage IV Ovarian Endometriosis of fertile women",
     "Fertile Secretory with Stage IV Ovarian Endometriosis (Group 2A)", 1),
    ("infertile_control", "Disease free Endometrium of Infertile women",
     "Infertile Secretory without Endometriosis (Group 1B)", 0),
]


def write_table(path, rows):
    with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
        csv.writer(f, delimiter="\t", lineterminator="\n").writerows(rows)
    return path


def matrix_rows():
    samples = [(rep, cell) for rep in range(1, 10) for cell in CELLS]
    ids = [f"GSM9900{i:03d}" for i in range(36)]
    return [
        ["!Series_geo_accession", "GSE120103"], ["!Series_platform_id", "GPL6480"],
        ["!Sample_geo_accession", *ids], ["!Sample_platform_id", *(["GPL6480"] * 36)],
        ["!Sample_channel_count", *(["1"] * 36)],
        ["!Sample_source_name_ch1", *["human endometrium, " + c[2] for _, c in samples]],
        ["!Sample_characteristics_ch1", *(["tissue: Endometrium"] * 36)],
        ["!Sample_characteristics_ch1", *["sample group: " + c[1] for _, c in samples]],
        ["!Sample_title", *[f"{c[1]}, replicate {rep}" for rep, c in samples]],
        ["!series_matrix_table_begin"], ["ID_REF", *ids],
        ["REF0", *[0] * 36], ["REF1", *[1] * 36], ["REF2", *[2] * 36],
        ["INCOMPLETE", "NA", *[-100] * 35],
        [fosb.PROBE, *[c[3] + int(rep % 3 == 0) for rep, c in samples]],
        ["!series_matrix_table_end"],
    ]


def annotation_rows():
    return [["!platform_table_begin"], ["ID", "Gene symbol", "Gene ID"],
            [fosb.PROBE, "FOSB", "2354"], ["!platform_table_end"]]


def row_for(rows, key):
    return next(row for row in rows if row[0] == key)


def test_full_analysis_and_common_reference(tmp_path, monkeypatch):
    matrix = write_table(tmp_path / "matrix.gz", matrix_rows())
    annot = write_table(tmp_path / "platform.gz", annotation_rows())
    with pytest.raises(ValueError, match="SHA-256"):
        fosb.analyze(matrix, annot)
    monkeypatch.setattr(fosb, "MATRIX_HASH", sha256(matrix.read_bytes()).hexdigest())
    monkeypatch.setattr(fosb, "ANNOTATION_HASH", sha256(annot.read_bytes()).hexdigest())
    monkeypatch.setattr(fosb, "MIN_COMMON_FEATURES", 4)
    cells, total, complete = fosb.read_percentiles(matrix)
    # For X in [-1, 0, 1, 2, 3], relative to [0, 1, 2, X], midrank percentiles
    # are 12.5, 25, 50, 75 and 87.5. Missing in one array excludes the whole row.
    values = {-1: 12.5, 0: 25, 1: 50, 2: 75, 3: 87.5}
    assert (total, complete) == (5, 4)
    for group, _, _, value in CELLS:
        np.testing.assert_equal(cells[group], [values[value], values[value], values[value + 1]] * 3)
    result = fosb.analyze(matrix, annot)
    assert result["features_tested"] == 1 and result["excluded_features"] == 1
    assert result["counts"] == {g: 9 for g, *_ in CELLS}
    assert result["stratum_differences_percentile_points"] == pytest.approx(
        {"fertile": 125 / 3, "infertile": 275 / 6})
    assert result["difference_percentile_points"] == pytest.approx(43.75)
    assert result["decision"] == "supports_expected_direction"


@pytest.mark.parametrize(("effect", "decision"), [(5, "supports_expected_direction"),
                                                  (-5, "supports_opposite_direction"),
                                                  (0, "inconclusive")])
def test_stratified_contrast_against_closed_form(effect, decision):
    baseline = np.arange(9, dtype=float) - 4
    cells = {"fertile_case": baseline + effect - 2, "fertile_control": baseline,
             "infertile_case": baseline + effect + 4, "infertile_control": baseline + 2}
    result = fosb.summarize(cells)
    assert result["standard_error_percentile_points"] == pytest.approx(sqrt(5 / 6))
    assert result["welch_df"] == pytest.approx(32)
    margin = t.ppf(0.975, 32) * sqrt(5 / 6)
    assert result["ci95_percentile_points"] == pytest.approx([effect - margin, effect + margin])
    assert result["two_sided_p"] == pytest.approx(2 * t.sf(abs(effect) / sqrt(5 / 6), 32))
    assert result["leave_one_out_range_percentile_points"] == pytest.approx([effect - .25, effect + .25])
    assert result["decision"] == decision


def test_stop_on_zero_contrast_variance():
    with pytest.raises(ValueError, match="variance"):
        fosb.summarize({g: np.ones(9) for g, *_ in CELLS})


@pytest.mark.parametrize(("key", "replacement", "error"), [
    ("!Series_geo_accession", "GSE000", "series"),
    ("!Sample_platform_id", "GPL0", "platform"),
    ("!Sample_channel_count", "2", "single-channel"),
    ("!Sample_geo_accession", "GSM9900001", "distinct"),
    ("!Sample_title", CELLS[0][1] + ", replicate 2", "replicates"),
    ("!Sample_title", CELLS[2][1] + ", replicate 1", "title"),
    ("!Sample_characteristics_ch1", "tissue: Myometrium", "tissue"),
    ("!Sample_source_name_ch1", "human endometrium, Infertile Secretory with Stage IV Ovarian Endometriosis (Group 2B)", "source"),
    ("ID_REF", "GSM9900001", "sample order"),
    (fosb.PROBE, "null", "incomplete"),
    ("INCOMPLETE", "junk", "Malformed submitted"),
])
def test_reject_invalid_fields(tmp_path, key, replacement, error):
    rows = matrix_rows()
    row_for(rows, key)[1] = replacement
    with pytest.raises(ValueError, match=error):
        fosb.read_percentiles(write_table(tmp_path / "input.gz", rows))


@pytest.mark.parametrize("missing", ["", "NA", "NaN", "null", "inf", "-Infinity"])
def test_missing_and_nonfinite_features_excluded_from_all_arrays(tmp_path, monkeypatch, missing):
    rows = matrix_rows()
    row_for(rows, "INCOMPLETE")[1] = missing
    monkeypatch.setattr(fosb, "MIN_COMMON_FEATURES", 4)
    assert fosb.read_percentiles(write_table(tmp_path / "matrix.gz", rows))[1:] == (5, 4)


@pytest.mark.parametrize("problem", ["duplicate", "missing_target", "wrong_width", "blank_key",
                                     "no_end", "extra_row", "gzip", "few_features"])
def test_reject_invalid_matrix_structure(tmp_path, problem):
    rows = matrix_rows()
    if problem == "duplicate":
        rows.insert(-1, row_for(rows, "REF1"))
    elif problem == "missing_target":
        rows.remove(row_for(rows, fosb.PROBE))
    elif problem == "wrong_width":
        row_for(rows, "REF1").pop()
    elif problem == "blank_key":
        row_for(rows, "REF1")[0] = " "
    elif problem == "no_end":
        rows.pop()
    elif problem == "extra_row":
        rows.append(["trailing"])
    path = write_table(tmp_path / "bad.gz", rows)
    if problem == "gzip":
        path.write_bytes(path.read_bytes()[:-5])
    errors = {"duplicate": "duplicate", "missing_target": "missing or incomplete",
              "wrong_width": "Invalid feature", "blank_key": "Invalid feature",
              "no_end": "end of matrix", "extra_row": "trailing", "few_features": "Insufficient"}
    with pytest.raises((EOFError, gzip.BadGzipFile) if problem == "gzip" else ValueError,
                       match=None if problem == "gzip" else errors[problem]):
        fosb.read_percentiles(path)


@pytest.mark.parametrize("other", [["SECOND", "ANOTHER /// FOSB", "12 /// 2354"],
                                   ["SECOND", "OTHER", "2354"], [fosb.PROBE, "FOSB", "2354"]])
def test_no_ambiguous_annotation(tmp_path, other):
    rows = annotation_rows()
    rows.insert(-1, other)
    with pytest.raises(ValueError, match="probe mapping"):
        fosb.check_annotation(write_table(tmp_path / "platform.gz", rows))


def test_percentiles_invariant_to_column_specific_increasing_transforms(tmp_path, monkeypatch):
    rows = matrix_rows()
    monkeypatch.setattr(fosb, "MIN_COMMON_FEATURES", 4)
    original, _, _ = fosb.read_percentiles(write_table(tmp_path / "original.gz", rows))
    for row in rows[rows.index(row_for(rows, "ID_REF")) + 1:-1]:
        row[1:] = [entry if entry == "NA" else (i + 1) * float(entry) + 4 * i
                   for i, entry in enumerate(row[1:])]
    changed, _, _ = fosb.read_percentiles(write_table(tmp_path / "transformed.gz", rows))
    for group in original:
        np.testing.assert_array_equal(original[group], changed[group])


def test_unequal_variances_use_satterthwaite_weights():
    # Cell variances 0, 1, 4, 9, with three observations per cell.
    cells = {"fertile_control": np.array([1, 1, 1]), "fertile_case": np.array([1, 2, 3]),
             "infertile_control": np.array([1, 3, 5]), "infertile_case": np.array([3, 6, 9])}
    result = fosb.summarize(cells)
    assert result["difference_percentile_points"] == pytest.approx(2)
    assert result["standard_error_percentile_points"] == pytest.approx(sqrt(14 / 12))
    assert result["welch_df"] == pytest.approx(4)
    assert result["leave_one_out_range_percentile_points"] == pytest.approx([1.25, 2.75])


def test_malformed_token_rejected_even_in_an_incomplete_feature(tmp_path):
    rows = matrix_rows()
    row_for(rows, "INCOMPLETE")[-1] = "typo"
    with pytest.raises(ValueError, match="Malformed submitted"):
        fosb.read_percentiles(write_table(tmp_path / "bad.gz", rows))
