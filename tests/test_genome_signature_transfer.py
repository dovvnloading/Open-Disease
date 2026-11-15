import gzip
import hashlib
from pathlib import Path
import struct
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from fetch_genome_transfer_sources import reconcile_ena
from genome_signature_transfer import score_profiles, transfer
from rome_hybridization_controls import CONTROL
from rome_scan_design import scan_fields


def test_score_is_aligned_and_normalized_on_measured_overlap():
    genes = np.array([30, 10, 100])
    values = np.array([[1., 10., 999.], [20., 4., -555.]])
    ids = np.array([10, 20, 30])
    weights = np.array([2., 100., -1.])
    n, retained, scores = score_profiles(genes, values, ids, weights, {10})
    assert (n, retained) == (2, 1)
    np.testing.assert_allclose(scores[:,0], [19/3, -4])
    np.testing.assert_allclose(scores[:,1], [-1, -20])
    np.testing.assert_allclose(scores[:,2], [4.5, -8])


def test_equal_day_score_excludes_single_group_day():
    y = np.array([0, 1, 2, 8, 12, 99], dtype=float)[:,None]
    cases = np.array([False, False, True, False, True, True])
    blocks = np.array(['A', 'A', 'A', 'B', 'B', 'C'])
    r = transfer(np.array([1]), y, cases, blocks, np.array([1]), np.ones(1), set())
    assert r['score_difference_pp'] == pytest.approx((1.5 + 4)/2)
    assert (r['cases'], r['controls'], r['strata'], r['assignments']) == (2, 3, 2, 6)
    assert r['admissible_deletions'] == 2


def test_native_scan_parser_checks_geometry_and_timestamp():
    hdr = b'Cols=732\nDatHeader=field HG-U133A_2.1sq 03/17/10 12:34:56 123  M10\n'
    payload = struct.pack('<6i', 64, 4, 732, 732, 535824, len(hdr)) + hdr
    scanner, dt = scan_fields(gzip.compress(payload).ljust(4096, b'\x00'))
    assert scanner == '123' and dt.isoformat() == '2010-03-17T12:34:56'
    with pytest.raises(ValueError, match='4 KiB'):
        scan_fields(gzip.compress(payload))


def test_ena_reorders_but_never_repairs_content():
    lines = [b'header\n'] + [f'row{i}\n'.encode() for i in range(7)]
    canonical = b''.join(lines)
    digest = hashlib.sha256(canonical).hexdigest()
    shuffled = lines[0] + b''.join(reversed(lines[1:]))
    assert reconcile_ena(shuffled, digest) == canonical
    with pytest.raises(ValueError, match='beyond a row rearrangement'):
        reconcile_ena(shuffled.replace(b'row3', b'row8'), digest)


def test_external_controls_exclude_human_reference_and_labeling_controls():
    accepted = ('AFFX-BioDn-5_at', 'AFFX-r2-Ec-bioC-3_at', 'AFFX-r2-P1-cre-3_at', 'AFFX-CreX-5_at')
    excluded = ('AFFX-HUMGAPDH/M33197_3_at', 'AFFX-r2-Bs-lys-5_at', 'AFFX-DapX-3_at')
    assert all(CONTROL.fullmatch(p) for p in accepted)
    assert not any(CONTROL.fullmatch(p) for p in excluded)
