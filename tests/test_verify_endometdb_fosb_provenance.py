"""Synthetic aggregate-only regressions; never load expression or donor records."""
import copy
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from verify_endometdb_fosb_provenance import evaluate


@pytest.fixture
def manifest():
    return json.loads((ROOT / "evidence" / "endometdb-fosb-measurement-audit" /
                       "measurement-provenance.json").read_text(encoding="utf-8"))


def test_combined_row_is_blocked_but_native_route_is_distinct(manifest):
    result = evaluate(manifest)
    assert result["pooled_cross_platform_decision"] == "NO_GO"
    assert result["ht12_case_control_decision"] == "NO_GO_ZERO_GEO_PRIMARY_CONTROLS"
    assert (result["native_wg6_shared_strata"], result["native_wg6_shared_cases"],
            result["native_wg6_shared_controls"]) == (7, 12, 13)
    assert result["new_expression_values_inspected"] is False
    assert len(result["blockers"]) == 3


@pytest.mark.parametrize("change, error", [
    (lambda d: d["physical_probes"]["GPL13376"].update(probe="ILMN_1751607"), "GPL13376"),
    (lambda d: d["physical_probes"]["GPL13376"].update(address=7160239), "GPL13376"),
    (lambda d: d.update(probe_oligos_distinct=False), "Distinct"),
    (lambda d: d["combined_workbook"].update(reported_fosb_probe="ILMN_1730706"), "Combined"),
    (lambda d: d["combined_workbook"].update(fosb_rows=2), "Combined"),
    (lambda d: d["native_cohort"].update(complete_wg6=31), "source-pinned"),
    (lambda d: d["native_cohort"]["shared_run_phase_cells"][0].update(controls=0), "shared cell"),
    (lambda d: d["native_cohort"]["shared_run_phase_cells"][0].update(cases=3), "totals"),
    (lambda d: d["native_cohort"]["shared_run_phase_cells"][0].update(phase="mid-secretory"), "histology"),
    (lambda d: d["native_cohort"]["shared_run_phase_cells"][0].update(run=11,phase="proliferative"), "Repeated"),
    (lambda d: d["native_cohort"].update(hormone_medication_lookback="3 months"), "washout"),
    (lambda d: d["native_cohort"].update(no_donor_cross_release_matching=False), "donor"),
    (lambda d: d.update(ht12_geo_primary_endometrium_controls=1), "HT-12"),
    (lambda d: d["sources"].update(geo="http://invalid"), "HTTPS"),
    (lambda d: d["native_cohort"]["shared_run_phase_cells"][0].update(cases=True), "integer"),
])
def test_source_provenance_or_cohort_mismatch_fails_closed(manifest, change, error):
    changed = copy.deepcopy(manifest)
    change(changed)
    with pytest.raises(ValueError, match=error):
        evaluate(changed)


def test_new_documentation_cannot_automatically_authorize_results(manifest):
    manifest["combined_workbook"]["documented_wg6_to_pooled_fosb_mapping"] = (
        "https://example.org/new-source")
    result = evaluate(manifest)
    assert result["pooled_cross_platform_decision"] == "RE_REVIEW_REQUIRED"
    assert result["new_expression_values_inspected"] is False


def test_provenance_evidence_requires_a_real_url(manifest):
    manifest["combined_workbook"]["documented_native_to_geo_assay_crosswalk"] = "inferred by age"
    with pytest.raises(ValueError, match="HTTPS"):
        evaluate(manifest)


def test_missing_wg6_evidence_fails_closed(manifest):
    del manifest["physical_probes"]["GPL13376"]
    with pytest.raises(ValueError, match="Both GEO platforms"):
        evaluate(manifest)
