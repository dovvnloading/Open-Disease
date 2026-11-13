"""GSE51981 metadata preflight tests: entirely artificial accession/clinical records."""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_gse51981_metadata import audit, check_manifest, _read_metadata  # noqa: E402

SAMPLES = [f"GSM{900000 + i}" for i in range(148)]


def tsv(*fields):
    s = io.StringIO()
    csv.writer(s, delimiter="\t", lineterminator="\n").writerow(fields)
    return s.getvalue()


def header(*, samples=SAMPLES, platforms=None, series="GSE51981",
           missing_char=False, bad_row=False):
    platforms = platforms or ["GPL570"] * 148
    lines = [tsv("!Series_geo_accession", series),
             tsv("!Series_platform_id", "GPL570"),
             tsv("!Sample_geo_accession", *samples),
             tsv("!Sample_platform_id", *platforms),
             tsv("!Sample_title", *["case" for _ in range(148)]),
             tsv("!Sample_characteristics_ch1", *["tissue: Early secretory endometrial tissue" for _ in range(148)]),
             tsv("!Sample_characteristics_ch1", *["endometriosis/no endometriosis: case (UNVERIFIED)" for _ in range(148)]),
             tsv("!Sample_characteristics_ch1", *["presence or absence of uterine/pelvic pathology: no" for _ in range(148)]),
             tsv("!Sample_data_processing", *["GCRMA" for _ in range(148)])]
    if missing_char:
        lines[6] = tsv("!Sample_characteristics_ch1", *["no key value" for _ in range(148)])
    if bad_row:
        lines[6] = tsv("!Sample_characteristics_ch1", "only-one")
    return "".join(lines) + "!series_matrix_table_begin\n" + "DONT_PARSE_MATRIX_INTENSITIES\n"


def write_source(tmp_path, payload=None, *, compressed=False):
    path = tmp_path / ("matrix.gz" if compressed else "matrix.txt")
    text = (payload if payload is not None else header()).encode("utf-8")
    path.write_bytes(gzip.compress(text) if compressed else text)
    return path


def write_manifest(tmp_path, *, included=0):
    from audit_gse51981_metadata import REQUIRED_MANIFEST
    fields = sorted(REQUIRED_MANIFEST)
    rows = []
    for i, gsm in enumerate(SAMPLES):
        inc = i < included
        d = dict.fromkeys(fields, "")
        d.update({"geo_accession": gsm, "donor_id": f"pseudonymous-SYN{i}",
                  "include": "true" if inc else "false",
                  "reviewed": "true" if inc else "false",
                  "disease_group": "case" if i % 2 else "control",
                  "phase": "proliferative", "tissue": "eutopic_endometrium",
                  "pathology": "SYNTHETIC no actual pathology", "hormone_treatment": "none",
                  "exclusion_reason": "" if inc else "pending_source_review",
                  "review_source_url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}"})
        rows.append(d)
    path = tmp_path / "private_SYNTHETIC_manifest.csv"
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return path, rows, fields


def rewrite_manifest(path, rows, columns):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


@pytest.mark.parametrize("compressed", [False, True])
def test_metadata_only_inspects_headers_without_expression_or_clinical_inference(tmp_path, compressed):
    source = write_source(tmp_path, compressed=compressed)
    report = audit(source)
    assert report["n_geo_accessions"] == 148
    assert report["source_file_sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    assert report["source_field_coverage"]["source_field_presence_n"] == {
        "tissue": 148, "endometriosis": 148, "other_pathology": 148,
        "explicit_cycle_phase": 0, "explicit_person_id": 0,
        "hormone_or_treatment": 0, "data_processing": 148,
    }
    assert report["source_rights_verified"] is False
    assert report["donor_independence_verified_by_script"] is False
    assert report["expression_intensities_parsed_or_analyzed"] is False
    assert report["biological_findings"] is False
    assert "case (UNVERIFIED)" not in json.dumps(report)
    assert report["preanalysis_approval"].startswith("NO_GO_")


@pytest.mark.parametrize("payload,match", [
    (header(series="GSE1234"), "Series accession"),
    (header(samples=SAMPLES[:-1]), "148 unique"),
    (header(samples=SAMPLES[:-1] + [SAMPLES[0]]), "148 unique"),
    (header(platforms=["GPL570"] * 147 + ["GPL96"]), "GPL570"),
    (header(bad_row=True), "misaligned"),
    (header().split("!series_matrix_table_begin")[0], "Missing matrix table begin"),
])
def test_corrupt_or_incorrect_source_is_rejected(tmp_path, payload, match):
    source = write_source(tmp_path, payload)
    with pytest.raises(ValueError, match=match):
        _read_metadata(source)


def test_unreviewed_manifest_reports_no_go_without_inventing_donor_labels(tmp_path):
    matrix = write_source(tmp_path)
    manifest, _, _ = write_manifest(tmp_path)
    report = audit(matrix, manifest=manifest)
    assert report["manifest"]["pending_review_n"] == 148
    assert report["manifest"]["included_n"] == 0
    assert report["manifest"]["source_label_authenticity_verified_by_script"] is False
    assert report["preanalysis_approval"].startswith("NO_GO_")


def test_structural_manifest_keeps_review_and_evidence_gates_separate(tmp_path):
    matrix = write_source(tmp_path)
    manifest, _, _ = write_manifest(tmp_path, included=5)
    report = audit(matrix, manifest=manifest)
    assert report["manifest"]["included_groups_n"] == {"case": 2, "control": 3, "other": 0}
    assert report["manifest"]["reviewed_n"] == 5
    assert report["preanalysis_approval"].startswith("NO_GO_")


@pytest.mark.parametrize("change,match", [
    (lambda x: x.update({"include": "true", "reviewed": "false"}), "not been reviewed"),
    (lambda x: x.update({"phase": "phase unknown"}), "unapproved group or phase"),
    (lambda x: x.update({"disease_group": "healthy"}), "unapproved group or phase"),
    (lambda x: x.update({"hormone_treatment": "unknown"}), "incompatible tissue or treatment"),
    (lambda x: x.update({"review_source_url": "https://example.org/unsupported"}), "cite its own primary GSM"),
    (lambda x: x.update({"donor_id": ""}), "mandatory review field"),
    (lambda x: x.update({"geo_accession": "GSM12345"}), "same 148"),
])
def test_manifest_fail_closed_when_source_review_cannot_support_inclusion(tmp_path, change, match):
    manifest, rows, fields = write_manifest(tmp_path, included=5)
    change(rows[0])
    rewrite_manifest(manifest, rows, fields)
    with pytest.raises(ValueError, match=match):
        check_manifest(manifest, SAMPLES)


def test_duplicate_donor_rejected_and_all_exclusions_have_reasons(tmp_path):
    manifest, rows, fields = write_manifest(tmp_path, included=5)
    rows[1]["donor_id"] = rows[0]["donor_id"]
    rewrite_manifest(manifest, rows, fields)
    with pytest.raises(ValueError, match="Duplicate donor"):
        check_manifest(manifest, SAMPLES)
    rows[1]["donor_id"] = "SYN1"
    rows[12]["exclusion_reason"] = ""
    rewrite_manifest(manifest, rows, fields)
    with pytest.raises(ValueError, match="All excluded"):
        check_manifest(manifest, SAMPLES)


def test_cli_writes_aggregate_only_to_ignored_results_and_refuses_overwrite(tmp_path):
    source = write_source(tmp_path)
    script = ROOT / "scripts" / "audit_gse51981_metadata.py"
    outside = tmp_path / "report.json"
    bad = subprocess.run([sys.executable, str(script), "--matrix", str(source),
                          "--out", str(outside)], capture_output=True, text=True)
    assert bad.returncode != 0
    assert not outside.exists()
    inside = ROOT / "results" / "synthetic-GSE51981-metadata-test.json"
    inside.parent.mkdir(exist_ok=True)
    try:
        inside.unlink(missing_ok=True)
        cmd = [sys.executable, str(script), "--matrix", str(source), "--out", str(inside)]
        ok = subprocess.run(cmd, capture_output=True, text=True)
        assert ok.returncode == 0, ok.stderr
        contents = json.loads(inside.read_text())
        assert contents["biological_findings"] is False
        assert "GSM900" not in json.dumps(contents)
        retry = subprocess.run(cmd, capture_output=True, text=True)
        assert retry.returncode != 0
    finally:
        inside.unlink(missing_ok=True)
