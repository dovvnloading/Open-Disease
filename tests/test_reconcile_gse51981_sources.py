"""Artificial records check the metadata crosswalk and its aggregate-only output."""
import csv
import io
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from reconcile_gse51981_sources import reconcile


def sources(tmp_path, edit=lambda s: s):
    ids = [f"GSM{900000 + n}" for n in range(148)]
    phases = [
        "Late Secretory Endometrial tissue", "Unknown", "unrecognized SYNTHETIC term",
        "Proliferative Endometrial tissue", "Early Secretory Endometrial tissue",
        "Mid-Secretory Endometrial tissue",
    ]
    groups = ["Non-Endometriosis", "Endometriosis", "Non-Endometriosis",
              "Non-Endometriosis", "Endometriosis", "Non-Endometriosis"]
    pathologies = [
        "presence or absence of uterine/pelvic pathology: Uterine Pelvic Pathology",
        "endometriosis severity: Minimal/Mild",
        "presence or absence of uterine/pelvic pathology: unknown SYNTHETIC pathology",
        "presence or absence of uterine/pelvic pathology: No Uterine Pelvic Pathology",
        "endometriosis severity: Moderate/Severe",
        "presence or absence of uterine/pelvic pathology: Uterine Pelvic Pathology",
    ]
    field_rows = [
        ("!Sample_geo_accession", ids),
        ("!Sample_platform_id", ["GPL570"] * 148),
        ("!Sample_source_name_ch1", ["Endometriosis_Severity Not Available"] * 148),
        ("!Sample_characteristics_ch1", [f"tissue: {phases[i % 6]}" for i in range(148)]),
        ("!Sample_characteristics_ch1", [f"endometriosis/no endometriosis: {groups[i % 6]}" for i in range(148)]),
        ("!Sample_characteristics_ch1", [pathologies[i % 6] for i in range(148)]),
    ]
    header = io.StringIO()
    w = csv.writer(header, delimiter="\t", lineterminator="\n")
    w.writerow(["!Series_geo_accession", "GSE51981"])
    for k, v in field_rows:
        w.writerow([k, *v])
    matrix = tmp_path / "matrix.txt"
    matrix.write_text(header.getvalue() + "!series_matrix_table_begin\n", encoding="utf-8")
    soft = tmp_path / "quick.soft"
    lines = []
    for i in reversed(range(148)):
        lines.append("^SAMPLE = " + ids[i])
        lines.extend(k + " = " + v[i] for k, v in reversed(field_rows))
        lines.extend(["!Sample_series_id = GSE51981", "#VALUE = base-2 logarithm GC-RMA",
                      "!sample_table_begin", "!Sample_series_id = IGNORE_TABLE_CONTENT", "!sample_table_end"])
    soft.write_text(edit("\n".join(lines) + "\n"), encoding="utf-8")
    return matrix, soft


def test_reconciles_by_accession_and_key_even_if_reordered(tmp_path):
    counts = reconcile(*sources(tmp_path))
    assert counts["reconciled_arrays"] == counts["log2_gcrma_value_descriptions"] == 148
    assert counts["group_by_source_phase"] == {
        "case": {"unknown": 25, "early_secretory": 24}, "control": {"proliferative": 25},
        "other": {"late_secretory": 25, "mid_secretory": 24},
        "unresolved": {"unrecognized_or_missing": 25},
    }
    assert counts["cases_with_severity_characteristic_but_source_name_severity_unavailable"] == 49
    output = json.dumps(counts)
    assert "GSM" not in output and "SYNTHETIC" not in output


@pytest.mark.parametrize("before,after,error", [
    ("^SAMPLE = GSM900147", "^SAMPLE = GSM900145", "duplicate"),
    ("^SAMPLE = GSM900147", "^SAMPLE = GSM999999", "sets differ"),
    ("!Sample_geo_accession = GSM900147", "!Sample_geo_accession = GSM999999", "identifier mismatch"),
    ("!Sample_series_id = GSE51981", "!Sample_series_id = GSE11111", "identifier mismatch"),
    ("tissue: Unknown", "tissue: Late Secretory Endometrial tissue", "values differ"),
    ("!sample_table_end\n", "", "sets differ|Unterminated"),
])
def test_rejects_divergent_sources(tmp_path, before, after, error):
    with pytest.raises(ValueError, match=error):
        reconcile(*sources(tmp_path, lambda s: s.replace(before, after, 1)))
