"""Only synthetic GPL570-shaped probe annotations, no real human data/network."""
from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from audit_gpl570 import classify, inspect  # noqa: E402


def artificial_table(*, soft=False):
    data = [
        "ID\tGene Symbol\tENTREZ_GENE_ID\tOther\n",
        "PROBE_SYN_1\tFAKEA\t101\ttext\n",
        "PROBE_SYN_2\tFAKEA /// FAKEB\t101 /// 102\ttext\n",
        "PROBE_SYN_3\t---\t---\ttext\n",
        "PROBE_SYN_4\tFAKEC\t---\ttext\n",
        "PROBE_SYN_5\tFAKE C\t103\ttext\n",
        "PROBE_SYN_6\tFAKED\t201 202\ttext\n",
    ]
    return ("^PLATFORM = GPL570\n!platform_table_begin\n" if soft else "") + "".join(data) + ("!platform_table_end\n" if soft else "")


@pytest.mark.parametrize("gzip_input,soft", [(False, False), (True, True)])
def test_synthetic_counts_and_provenance(tmp_path, gzip_input, soft):
    path = tmp_path / ("platform.gz" if gzip_input else "platform.tsv")
    payload = artificial_table(soft=soft).encode()
    path.write_bytes(gzip.compress(payload) if gzip_input else payload)
    result = inspect(path)
    assert result["total_probes"] == 6
    assert sum(result["classification"].values()) == 6
    assert result["classification"] == {
        "single_symbol_single_entrez": 1,
        "multi_gene": 2,
        "unmapped": 1,
        "incomplete": 1,
        "invalid_annotation": 1,
    }
    assert len(result["input_sha256"]) == 64
    assert result["source_verified"] is False
    assert result["declared_platform_id"] == ("GPL570" if soft else None)
    assert result["soft_platform_declaration_checked"] is soft
    assert result["biological_findings"] is False
    assert "PROBE_SYN" not in json.dumps(result)


def test_no_guessing_for_ambiguous_or_missing_field():
    assert classify("DDR1 /// MIR4640", "780 100616237") == "multi_gene"
    assert classify("FAKEA", "7 /// 8") == "multi_gene"
    assert classify("---", "") == "unmapped"
    assert classify("FAKEA", "NA") == "incomplete"
    assert classify("FAKEA/FAKEB", "101") == "invalid_annotation"


@pytest.mark.parametrize("payload,error", [
    ("ID\tGene Symbol\nP\tA\n", "Expected unique"),
    ("ID\tGene Symbol\tENTREZ_GENE_ID\nP\tA\t1\nP\tB\t2\n", "duplicate"),
    ("ID\tGene Symbol\tENTREZ_GENE_ID\nP\tA\n", "Malformed"),
    ("!platform_table_begin\nID\tGene Symbol\tENTREZ_GENE_ID\nP\tA\t1\n", "incomplete"),
])
def test_fail_closed_corrupt_platform(tmp_path, payload, error):
    path = tmp_path / "bad.tsv"
    path.write_text(payload)
    with pytest.raises(ValueError, match=error):
        inspect(path)


def test_cli_no_overwrite_and_only_ignored_results(tmp_path):
    source = tmp_path / "annotation.tsv"
    source.write_text(artificial_table())
    script = ROOT / "scripts" / "audit_gpl570.py"
    outside = tmp_path / "report.json"
    fail = subprocess.run([sys.executable, str(script), "--input", str(source),
                           "--out", str(outside)], capture_output=True, text=True)
    assert fail.returncode != 0
    assert not outside.exists()
    output = ROOT / "results" / "audit-gpl570-synthetic-test.json"
    output.parent.mkdir(exist_ok=True)
    try:
        output.unlink(missing_ok=True)
        command = [sys.executable, str(script), "--input", str(source), "--out", str(output)]
        ok = subprocess.run(command, capture_output=True, text=True)
        assert ok.returncode == 0, ok.stderr
        assert output.exists()
        retry = subprocess.run(command, capture_output=True, text=True)
        assert retry.returncode != 0
    finally:
        output.unlink(missing_ok=True)


@pytest.mark.parametrize("prefix,error", [
    ("^PLATFORM = GPL96\n!platform_table_begin\n", "Unexpected SOFT platform"),
    ("^PLATFORM = GPL570\n^PLATFORM = GPL570\n!platform_table_begin\n", "Repeated SOFT platform"),
    ("^PLATFORM = GPL570\n", "Missing SOFT platform table begin"),
    ("^PLATFORM = bogus\n!platform_table_begin\n", "Malformed SOFT platform"),
])
def test_sof_platform_identity_fails_closed(tmp_path, prefix, error):
    path = tmp_path / "wrong-platform.soft"
    path.write_text(prefix + artificial_table())
    with pytest.raises(ValueError, match=error):
        inspect(path)


def test_cli_expected_sha256_rejects_wrong_or_invalid_pins(tmp_path):
    source = tmp_path / "synthetic.soft"
    source.write_text(artificial_table(soft=True))
    real_digest = hashlib.sha256(source.read_bytes()).hexdigest()
    script = ROOT / "scripts" / "audit_gpl570.py"
    output = ROOT / "results" / "audit-gpl570-pin-test.json"
    output.parent.mkdir(exist_ok=True)
    try:
        output.unlink(missing_ok=True)
        base = [sys.executable, str(script), "--input", str(source), "--out", str(output),
                "--expected-sha256"]
        for bad in ("xyz", "0" * 64):
            failed = subprocess.run(base + [bad], capture_output=True, text=True)
            assert failed.returncode != 0
            assert not output.exists()
        passed = subprocess.run(base + [real_digest.upper()], capture_output=True, text=True)
        assert passed.returncode == 0, passed.stderr
        report = json.loads(output.read_text())
        assert report["input_sha256"] == real_digest
        assert report["source_verified"] is False  # A matching hash cannot verify release/rights.
    finally:
        output.unlink(missing_ok=True)


def artificial_geo_annotation():
    return (
        "^Annotation\n"
        "!Annotation_date = Aug 09 2016\n"
        "!Annotation_platform = GPL570\n"
        "!Annotation_platform_title = Synthetic annotation\n"
        "#Gene ID = Entrez Gene identifier\n"
        "!platform_table_begin\n"
        "Gene ID\tID\tGene symbol\n"
        "33\tSYN1\tFAKEA\n"
        "33///44\tSYN2\tFAKEA///FAKEB\n"
        "\tSYN3\t\n"
        "44///\tSYN4\tFAKEB\n"
        "44///33\tSYN5\tFAKE@///FAKEA\n"
        "33\tSYN6\tFAKEA\n"
        "!platform_table_end\n"
    )


def test_geo_generated_annotation_uses_declared_format_and_release(tmp_path):
    path = tmp_path / "annotation.gz"
    path.write_bytes(gzip.compress(artificial_geo_annotation().encode()))
    report = inspect(path)
    assert report["table_schema"] == "geo_annotated_platform"
    assert report["declared_annotation_date"] == "2016-08-09"
    assert report["declared_platform_id"] == "GPL570"
    assert report["soft_platform_declaration_checked"]
    assert report["total_probes"] == 6
    assert report["distinct_entrez_ids_among_single_mapped_probes"] == 1
    assert report["probes_with_multiple_valid_entrez_ids"] == 2
    assert report["classification"] == {
        "single_symbol_single_entrez": 2, "multi_gene": 1, "unmapped": 1,
        "incomplete": 0, "invalid_annotation": 2,
    }


@pytest.mark.parametrize("original,replacement", [
    ("GPL570", "GPL96"),
    ("!Annotation_platform = GPL570\n", ""),
    ("!Annotation_platform = GPL570", "^PLATFORM = GPL570"),
    ("!Annotation_platform = GPL570", "!Annotation_platform = GPL570\n!Annotation_platform = GPL570"),
    ("!Annotation_date = Aug 09 2016\n", ""),
    ("Aug 09 2016", "unknown"),
    ("^Annotation\n", ""),
    ("Gene ID\tID\tGene symbol", "ENTREZ_GENE_ID\tID\tGene Symbol"),
    ("!platform_table_end\n", ""),
    ("!platform_table_end\n", "!platform_table_end\n^Annotation\n"),
])
def test_fail_closed_inconsistent_geo_annotation(tmp_path, original, replacement):
    path = tmp_path / "annotation.soft"
    path.write_text(artificial_geo_annotation().replace(original, replacement))
    with pytest.raises(ValueError):
        inspect(path)
