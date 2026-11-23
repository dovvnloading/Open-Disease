"""Replay the published aggregate comparison and check the two tipping points."""

import json
from pathlib import Path
import subprocess
import sys


EVIDENCE = Path(__file__).resolve().parents[1] / "evidence" / "outflow-obstruction-endometriosis"


def test_reproduce_published_contrasts():
    output = subprocess.check_output([sys.executable, str(EVIDENCE / "calculate.py")])
    assert output == (EVIDENCE / "results.json").read_bytes()

    rows = json.loads(output)
    assert len(rows) == 5
    thresholds = {
        row["study_id"]: row["tipping_point"]["additional_reference_cases"]
        for row in rows if row["outcome"] == "endometriosis"
    }
    assert thresholds == {
        "beijing-peking-union-obstructed-hemivagina": 17,
        "london-unicornuate-uterus-ultrasound": 16,
    }
