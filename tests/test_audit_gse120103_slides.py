"""Offline slide-confounding and bounded-header extraction checks."""
import gzip
import hashlib

import pytest

from scripts import audit_gse120103_slides as slides


def prefix(barcode, protocol, slide=None):
    slide = slide or barcode[:12]
    name = f"NA_{slide}_S01"
    # Pad the compressed prefix past the FEPARAMS record, as in an Agilent file.
    text = (f"TYPE\ttext\ttext\nFEPARAMS\tFeatureExtractor_Barcode\tProtocol_Name\tFeatureExtractor_ArrayName\tFeatureExtractor_ScanFileName\n"
            f"DATA\t{barcode}\t{protocol}\t{name}\t{name}.tif\n*\nFEATURES\tignored\n")
    return gzip.compress(text.encode()).ljust(slides.PREFIX_SIZE, b"\x00")


def test_header_identity_and_missing_fields():
    barcode = "000000000001_1_2"
    assert slides.original_slide(prefix(barcode, "proto")) == ("000000000001", barcode, "proto")
    with pytest.raises(ValueError, match="disagree"):
        slides.original_slide(prefix(barcode, "proto", slide="000000000002"))
    with pytest.raises(ValueError, match="Missing"):
        slides.original_slide(prefix("", "test-protocol"))
    with pytest.raises(ValueError, match="Incomplete"):
        slides.original_slide(prefix("slide", "protocol")[:-2])


def test_audit_requires_distinct_barcodes_to_establish_overlap(tmp_path, monkeypatch):
    source = []
    for group in slides.GROUPS:
        for i in range(9):
            url = f"synthetic/{group}/{i}"
            source.append((group, url))
            # Real overlap only between the first control and first case.
            k = slides.GROUPS.index(group) + 1
            if group in ("1A", "2A") and i == 0:
                barcode = f"000000000000_{k}_1"
            else:
                barcode = f"{k:012d}_1_{i}"
            protocol = "first" if i < 5 else "second"
            path = tmp_path / (hashlib.sha256(url.encode()).hexdigest() + ".prefix")
            path.write_bytes(prefix(barcode, protocol))
    monkeypatch.setattr(slides, "sample_headers", lambda _: source)
    result = slides.audit(None, tmp_path)
    assert result["unique_slides"] == 5
    assert result["slides_spanning_clinical_groups"] == 1
    assert result["slides_per_group"] == [2, 2, 1, 1]
    assert result["protocol_group_counts"] == {"first": [5] * 4, "second": [4] * 4}
    source.pop()
    with pytest.raises(ValueError, match="Incomplete"):
        slides.audit(None, tmp_path)
