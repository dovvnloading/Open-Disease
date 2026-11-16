"""Fail-closed, aggregate-only FOSB provenance decision for GSE141549.

This is a *source-annotation gate*, not a patient/assay data parser or an
expression test. It never authorizes a pooled analysis automatically.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


# Previously pinned manufacturer audits (see GSE141549_ELIGIBILITY_2026-11-13.md).
EXPECTED = {
    "GPL13376": {"probe": "ILMN_1730706", "address": 6040086,
                  "sha256": "c6f2c87ab7bc094b22270c907e1f9a6d74e9990b62726b53bf44660b0db700c1"},
    "GPL10558": {"probe": "ILMN_1751607", "address": 7160239,
                  "sha256": "770fb74bb90e548ee731924310d8f5c48589e63ee530dc11c7686d28a3ba5ba6"},
}
COMBINED_SHA256 = "9919fc5d6d6c09dbdae5a815da798b8250edf17a781fdb86dc25c2f9c9b21d50"
VALID_PHASES = {"proliferative", "secretory"}


def _count(value: object, name: str) -> int:
    if type(value) is not int or value < 0:
        raise ValueError(f"{name} must be a nonnegative integer")
    return value


def _source(value: object, name: str) -> None:
    if not isinstance(value, str) or not value.startswith("https://"):
        raise ValueError(f"{name} requires a public HTTPS provenance URL")


def evaluate(data: dict) -> dict:
    """Validate locked source observations and return an aggregate-only decision.

    Source URLs/fields here document the existing *published audit*. Adding a
    crosswalk URL does NOT authorize clinical or quantitative execution: it
    changes the outcome to RE_REVIEW_REQUIRED, never GO.
    """
    if not isinstance(data, dict) or data.get("schema_version") != 1:
        raise ValueError("Expected provenance schema_version 1")
    sources = data.get("sources")
    if not isinstance(sources, dict) or set(sources) != {"publication", "geo", "figshare"}:
        raise ValueError("Expected publication, GEO and versioned Figshare sources")
    for name, source in sources.items():
        _source(source, name)
    physical = data.get("physical_probes")
    if not isinstance(physical, dict) or set(physical) != set(EXPECTED):
        raise ValueError("Both GEO platforms must be present")
    for platform, expected in EXPECTED.items():
        if physical[platform] != expected:
            raise ValueError(f"{platform} physical probe/address or source pin changed")
    if data.get("probe_oligos_distinct") is not True:
        raise ValueError("Distinct manufacturer FOSB oligos must not be conflated")

    combined = data.get("combined_workbook")
    if not isinstance(combined, dict):
        raise ValueError("Missing combined workbook audit")
    if (combined.get("sha256") != COMBINED_SHA256 or
        combined.get("rows") != 28247 or
        combined.get("fosb_rows") != 1 or
        combined.get("fosb_excel_row") != 21565 or
        combined.get("reported_fosb_probe") != EXPECTED["GPL10558"]["probe"] or
        combined.get("wg6_probe_labels") != 0):
        raise ValueError("Combined workbook observations differ from the pinned audit")
    # A source string is not independently verified. It can ONLY reopen review.
    for key in ("documented_wg6_to_pooled_fosb_mapping", "documented_native_to_geo_assay_crosswalk"):
        value = combined.get(key)
        if value is not None:
            _source(value, key)

    native = data.get("native_cohort")
    if not isinstance(native, dict):
        raise ValueError("Missing native eligible cohort")
    if ([_count(native.get(key), key) for key in
         ("clinically_eligible", "complete_wg6", "duplicated_address_assays",
          "complete_wg6_cases", "complete_wg6_controls")] != [38, 32, 6, 15, 17]):
        raise ValueError("Native cohort counts conflict with the source-pinned census")
    if native["complete_wg6"] + native["duplicated_address_assays"] != native["clinically_eligible"]:
        raise ValueError("Clinically eligible cohort not reconciled")
    if native["complete_wg6_cases"] + native["complete_wg6_controls"] != native["complete_wg6"]:
        raise ValueError("WG-6 case/control counts not reconciled")
    cells = native.get("shared_run_phase_cells")
    if not isinstance(cells, list) or len(cells) != 7:
        raise ValueError("Expected seven shared native run-by-phase cells")
    visited = set()
    cases = controls = 0
    for cell in cells:
        if not isinstance(cell, dict) or set(cell) != {"run", "phase", "cases", "controls"}:
            raise ValueError("Unexpected run/phase cell keys")
        _count(cell["run"], "run")
        if cell["phase"] not in VALID_PHASES:
            raise ValueError("Unexpected histology phase")
        ncase, nctrl = _count(cell["cases"], "cases"), _count(cell["controls"], "controls")
        if ncase == 0 or nctrl == 0:
            raise ValueError("A shared cell needs cases and controls")
        key = (cell["run"], cell["phase"])
        if key in visited:
            raise ValueError("Repeated run-by-phase cell")
        visited.add(key)
        cases += ncase
        controls += nctrl
    if (cases, controls) != (12, 13):
        raise ValueError("Shared-cell totals must be 12 cases and 13 controls")
    if native.get("hormone_medication_lookback") != "unknown":
        raise ValueError("Cannot assume a known hormone washout period")
    if native.get("no_donor_cross_release_matching") is not True:
        raise ValueError("Cross-release donor matching must not be inferred")
    if _count(data.get("ht12_geo_primary_endometrium_controls"), "HT-12 GEO controls") != 0:
        raise ValueError("HT-12 case/control eligibility requires source re-audit")

    crosswalk_fields = ("documented_wg6_to_pooled_fosb_mapping", "documented_native_to_geo_assay_crosswalk")
    has_new_evidence = any(combined[key] is not None for key in crosswalk_fields)
    return {
        "task": "endometdb-fosb-measurement-audit",
        "pooled_cross_platform_decision": "RE_REVIEW_REQUIRED" if has_new_evidence else "NO_GO",
        "blockers": [] if has_new_evidence else [
            "No publicly documented WG-6-to-combined FOSB row derivation",
            "No source-backed native assay to combined GEO column crosswalk",
            "Gene symbol and HT-12 probe label do not establish WG-6 physical measurement",
        ],
        "ht12_case_control_decision": "NO_GO_ZERO_GEO_PRIMARY_CONTROLS",
        "native_wg6_distinct_probe": EXPECTED["GPL13376"]["probe"],
        "native_wg6_shared_strata": len(cells),
        "native_wg6_shared_cases": cases,
        "native_wg6_shared_controls": controls,
        "native_raw_route": "SEPARATELY_ANALYZED_NOT_A_POOLED_NORMALIZED_VALIDATION",
        "new_expression_values_inspected": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence", type=Path, help="Aggregate, source-pinned JSON; no patient fields")
    parser.add_argument("--require-no-go", action="store_true", help="Fail if a new source requires review")
    args = parser.parse_args()
    outcome = evaluate(json.loads(args.evidence.read_text(encoding="utf-8")))
    print(json.dumps(outcome, indent=2))
    if args.require_no_go and outcome["pooled_cross_platform_decision"] != "NO_GO":
        parser.exit(2, "New provenance needs independent review before any comparison.\n")


if __name__ == "__main__":
    main()
