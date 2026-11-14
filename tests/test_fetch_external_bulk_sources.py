from hashlib import sha256
from io import BytesIO
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import fetch_external_bulk_sources as acquire


RECORD = {"url": "https://example.org/frozen", "bytes": 3, "sha256": sha256(b"abc").hexdigest()}


def test_offline_first_and_verified_cache(tmp_path, monkeypatch):
    path = tmp_path / "new/entry"
    monkeypatch.setattr(acquire, "urlopen", lambda *a, **kw: pytest.fail("Unexpected network access"))
    with pytest.raises(FileNotFoundError, match="--fetch"):
        acquire.verify_or_fetch(path, RECORD, False)
    monkeypatch.setattr(acquire, "urlopen", lambda *a, **kw: BytesIO(b"abc"))
    assert acquire.verify_or_fetch(path, RECORD, True) == 3
    monkeypatch.setattr(acquire, "urlopen", lambda *a, **kw: pytest.fail("Unexpected cache refetch"))
    assert acquire.verify_or_fetch(path, RECORD, True) == 0
    path.write_bytes(b"bad")
    with pytest.raises(ValueError, match="Existing file"):
        acquire.verify_or_fetch(path, RECORD, True)
    assert path.read_bytes() == b"bad"


@pytest.mark.parametrize("payload", [b"ab", b"abcdefg", b"bad"])
def test_changed_download_never_installed(tmp_path, monkeypatch, payload):
    path = tmp_path / "entry"
    monkeypatch.setattr(acquire, "urlopen", lambda *a, **kw: BytesIO(payload))
    with pytest.raises(ValueError, match="Upstream response"):
        acquire.verify_or_fetch(path, RECORD, True)
    assert not path.exists()
