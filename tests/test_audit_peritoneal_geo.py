from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_peritoneal_geo import aggregate, group_map


def test_repeated_and_discordant_assays_do_not_become_controls_or_extra_pairs():
    rows = {}
    for title, role, tissue, stage, phase in [
        ("a", "Control endometrium", "Endometrium", "Healthy", "secretory"),
        ("b", "Control peritoneum", "Peritoneum", "Healthy", "proliferative"),
        ("c", "Control peritoneum", "Peritoneum", "4", "secretory"),
        ("b Replicate", "Control peritoneum", "Peritoneum", "Healthy", "secretory"),
        ("d", "Patient peritoneum", "Peritoneum", None, "secretory"),
    ]:
        rows[title] = {"description": role, "tissue": tissue, "disease stage": stage,
                       "cycle phase": phase, "platform_id": "array"}
    result = aggregate(rows, {k: "same" for k in rows})
    assert result["by_platform"]["array"]["link_groups_with_both_control_tissues_by_phase"] == {
        "discordant_or_multiple": 1}
    assert result["by_platform"]["array"]["arrays"]["control_peritoneum"] == 1
    assert result["exclusions_within_candidate_descriptions"] == {
        "missing_or_discordant_stage": 2, "marked_replicate": 1}


def test_workbook_requires_bijection_to_sample_titles():
    header = ("GEO Sample link", "GEO Sample code")
    with pytest.raises(ValueError, match="more than once"):
        group_map([header, ("group1", "a"), ("group2", "a")], {"a", "b"})
    with pytest.raises(ValueError, match="Incomplete"):
        group_map([header, ("group1", "a")], {"a", "b"})
