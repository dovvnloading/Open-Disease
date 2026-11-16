from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from normal_peritoneum_design import Pair
from normal_peritoneum_preflight import blocks, pair_lock


def test_existing_lock_required_and_disagreement_never_replaced(tmp_path):
    pair = Pair("fictional", "4", "spec-P", "spec-E", ("p",), ("u",), None, "unknown")
    path = tmp_path / "synthetic.json"
    with pytest.raises(ValueError, match="Missing or changed"):
        pair_lock(path, [pair], create=False)
    pair_lock(path, [pair], create=True)
    pair_lock(path, [pair], create=True)
    pair_lock(path, [pair], create=False)
    with pytest.raises(ValueError, match="changed"):
        pair_lock(path, [], create=True)
    pair_lock(path, [pair], create=False)


def test_bad_prefix_rejected_before_first_block(tmp_path, monkeypatch):
    import normal_peritoneum_preflight as p
    (tmp_path / "raw-prefix.zip.part").write_bytes(b"spurious")
    monkeypatch.setattr(p, "ZIP_PREFIX_BYTES", 8)
    with pytest.raises(ValueError, match="SHA-256 before stream"):
        next(blocks(tmp_path))
