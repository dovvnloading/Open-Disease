"""Synthetic-only GEO format and provenance tests; never contact NCBI."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

import pytest

from open_disease.geo_microarray import (
    MANIFEST_COLUMNS,
    SERIES_MATRIX_URL,
    download_series_matrix,
    prepare_series_matrix,
)


SAMPLES = [f"GSM90000{i}" for i in range(6)]


def fixture_bytes(*, header=None, rows=None, include_end=True) -> bytes:
    def tsv(*columns):
        out = io.StringIO()
        csv.writer(out, delimiter="\t", lineterminator="\n").writerow(columns)
        return out.getvalue()

    lines = [tsv("!Series_geo_accession", "GSE51981"),
             tsv("!Series_last_update_date", "Nov 12 2026"),
             tsv("!Sample_geo_accession", *SAMPLES),
             tsv("!Sample_title", *["case in source title" for _ in SAMPLES]),
             tsv("!Sample_characteristics_ch1", *["phase: source text" for _ in SAMPLES]),
             tsv("!Sample_characteristics_ch1", *["pathology: unverified" for _ in SAMPLES]),
             tsv("!Sample_platform_id", *["GPL570"] * 6),
             "!series_matrix_table_begin\n",
             tsv("ID_REF", *(header or SAMPLES))]
    for row in (rows if rows is not None else [
        ("SYN_PROBE_01", "4.1", "4.2", "", "5.0", "3.7", "4.4"),
        ("SYN_PROBE_02", "8", "9", "7", "8", "9", "7"),
    ]):
        lines.append(tsv(*row))
    if include_end:
        lines.append("!series_matrix_table_end\n")
    return "".join(lines).encode("utf-8")


def prepare(tmp_path, *, payload=None, expected_samples=6):
    raw = tmp_path / "data" / "raw" / "matrix.txt.gz"
    raw.parent.mkdir(parents=True, exist_ok=True)
    raw.write_bytes(gzip.compress(payload if payload is not None else fixture_bytes()))
    out = tmp_path / "data" / "processed"
    paths = [out / "manifest.csv", out / "intensities.tsv", out / "audit.json"]
    def run():
        return prepare_series_matrix(raw, manifest=paths[0], expression=paths[1],
                                     audit=paths[2], expected_samples=expected_samples)
    return raw, paths, run


def test_prepares_source_only_and_never_infers_phenotypes(tmp_path):
    raw, paths, run = prepare(tmp_path)
    report = run()
    assert (report["n_samples"], report["n_probes"], report["n_missing_intensities"]) == (6, 2, 1)
    assert report["source_file_sha256"] == hashlib.sha256(raw.read_bytes()).hexdigest()
    assert report["status"] == "metadata_unreviewed_no_analysis"
    with paths[0].open(newline="") as stream:
        manifest = list(csv.DictReader(stream))
    assert len(manifest) == 6
    assert list(manifest[0]) == MANIFEST_COLUMNS
    assert [r["geo_accession"] for r in manifest] == SAMPLES
    assert all(r["include"] == "false" and r["reviewed"] == "false" for r in manifest)
    assert all(r["donor_id"] == r["disease_group"] == r["phase"] == "" for r in manifest)
    assert manifest[0]["source_title"] == "case in source title"  # NOT a diagnosis
    assert json.loads(manifest[0]["source_characteristics_json"]) == {
        "!Sample_characteristics_ch1": ["phase: source text", "pathology: unverified"]
    }
    assert paths[1].read_text().splitlines()[1] == "SYN_PROBE_01\t4.1\t4.2\tNA\t5.0\t3.7\t4.4"
    assert json.loads(paths[2].read_text())["manual_review_required"] is True


@pytest.mark.parametrize("payload,match", [
    (fixture_bytes(header=SAMPLES[::-1]), "header differs"),
    (fixture_bytes(header=[SAMPLES[0]] * 6), "duplicate GSM"),
    (fixture_bytes(rows=[("p", "1", "2", "3", "4", "5", "6"),
                         ("p", "1", "2", "3", "4", "5", "6")]), "duplicate ID_REF"),
    (fixture_bytes(rows=[("p", "1", "2", "Infinity", "4", "5", "6")]), "Invalid numeric"),
    (fixture_bytes(rows=[("p", "1", "2")]), "row length"),
    (fixture_bytes(include_end=False), "Incomplete"),
])
def test_rejects_corrupt_or_ambiguous_matrix_without_partial_outputs(tmp_path, payload, match):
    _, paths, run = prepare(tmp_path, payload=payload)
    with pytest.raises(ValueError, match=match):
        run()
    assert not any(path.exists() for path in paths)


def test_sample_count_gate_and_no_overwrites(tmp_path):
    _, paths, run = prepare(tmp_path, expected_samples=148)
    with pytest.raises(ValueError, match="Unexpected count"):
        run()
    assert not any(path.exists() for path in paths)
    _, _, run = prepare(tmp_path)
    run()
    with pytest.raises(ValueError, match="already exists"):
        run()


def test_download_explicitly_records_checksum_and_prevents_overwrite(tmp_path, monkeypatch):
    response = gzip.compress(fixture_bytes())
    requested = []
    def fake_open(request, *, timeout):
        requested.append((request.full_url, timeout))
        return io.BytesIO(response)
    monkeypatch.setattr("open_disease.geo_microarray.urlopen", fake_open)
    dest = tmp_path / "data/raw/GSE51981_series_matrix.txt.gz"
    record = download_series_matrix(dest)
    assert requested == [(SERIES_MATRIX_URL, 90)]
    assert dest.read_bytes() == response
    assert record["sha256"] == hashlib.sha256(response).hexdigest()
    assert json.loads(Path(str(dest) + ".provenance.json").read_text()) == record
    with pytest.raises(ValueError, match="already exists"):
        download_series_matrix(dest)
    assert len(requested) == 1


def test_download_size_guard_cleans_partial_file(tmp_path, monkeypatch):
    monkeypatch.setattr("open_disease.geo_microarray.urlopen",
                        lambda *_args, **_kwargs: io.BytesIO(gzip.compress(fixture_bytes())))
    dest = tmp_path / "data/raw/matrix.txt.gz"
    with pytest.raises(ValueError, match="max_bytes"):
        download_series_matrix(dest, max_bytes=10)
    assert not dest.exists()
    assert not Path(str(dest) + ".provenance.json").exists()
    assert list(dest.parent.iterdir()) == []



def test_racing_provenance_is_not_deleted(tmp_path, monkeypatch):
    """If another process creates provenance after preflight, preserve its file."""
    dest = tmp_path / "data/raw/matrix.txt.gz"
    sidecar = Path(str(dest) + ".provenance.json")

    def fake_open(*_args, **_kwargs):
        sidecar.write_text("owned by another process")
        return io.BytesIO(gzip.compress(fixture_bytes()))

    monkeypatch.setattr("open_disease.geo_microarray.urlopen", fake_open)
    with pytest.raises(FileExistsError):
        download_series_matrix(dest)
    assert sidecar.read_text() == "owned by another process"
    assert not dest.exists()
    assert list(dest.parent.iterdir()) == [sidecar]


def test_racing_download_destination_is_never_overwritten(tmp_path, monkeypatch):
    """A different process wins the destination race after the initial check."""
    dest = tmp_path / "data/raw/matrix.txt.gz"
    sidecar = Path(str(dest) + ".provenance.json")

    def fake_open(*_args, **_kwargs):
        dest.write_bytes(b"competing owner's file")
        return io.BytesIO(gzip.compress(fixture_bytes()))

    monkeypatch.setattr("open_disease.geo_microarray.urlopen", fake_open)
    with pytest.raises(ValueError, match="already exists"):
        download_series_matrix(dest)
    assert dest.read_bytes() == b"competing owner's file"
    assert not sidecar.exists()


def test_racing_prepare_output_is_never_overwritten(tmp_path, monkeypatch):
    """Partially published owned files roll back; competing files stay intact."""
    _, paths, run = prepare(tmp_path)
    import open_disease.geo_microarray as geo
    actual_link = geo.os.link

    def racing_link(src, dest):
        if Path(dest) == paths[0]:
            paths[0].write_text("competing owner's manifest")
            raise FileExistsError("concurrent output")
        return actual_link(src, dest)

    monkeypatch.setattr(geo.os, "link", racing_link)
    with pytest.raises(ValueError, match="already exists"):
        run()
    assert paths[0].read_text() == "competing owner's manifest"
    assert not paths[1].exists()
    assert not paths[2].exists()
