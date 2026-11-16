"""Offline GSE141549 audit. Never interpret GEO link-group labels as native donor IDs."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path
import json

from openpyxl import load_workbook

from audit_endometdb import pinned
from genome_endometdb import load_wg6_annotation

BRIEF = "15b63510077781d8dca7ff9c08e3050f8101b29d463001c5968d4bea6359a394"
LINK = "44cad654a10df5def4b79356f5c0cb2845f66ee9b699372bd1dd0871bd215018"
PANEL = "4c634a65d648a83af8418e7ad1bf50e03d1a83982d0265d47d9ccbd06aba0114"
FIELDS = {"title", "description", "platform_id", "series_id", "characteristics_ch1"}
ROLES = {
    "Patient peritoneum": ("case_peritoneum", "Peritoneum"),
    "Control peritoneum": ("control_peritoneum", "Peritoneum"),
    "Control endometrium": ("control_endometrium", "Endometrium"),
    **{f"Peritoneal endometriosis {color} lesion": ("lesion", "Peritoneal endometriosis lesion")
       for color in ("red", "black", "white")},
}
ARMS = ("lesion", "case_peritoneum", "control_peritoneum", "control_endometrium")


def records(brief):
    result = {}
    for block in brief.split("^SAMPLE = ")[1:]:
        row = {}
        for line in block.splitlines():
            field, sep, value = line.partition(" = ")
            if sep and field.removeprefix("!Sample_") in FIELDS:
                key = field.removeprefix("!Sample_")
                if key == "characteristics_ch1":
                    key, value = value.split(": ", 1)
                if key in row:
                    raise ValueError("Duplicated SOFT metadata field")
                row[key] = value
        if row["series_id"] != "GSE141549" or row["title"] in result:
            raise ValueError("Unexpected series or duplicated sample title")
        result[row["title"]] = row
    return result


def group_map(rows, titles):
    iterator = iter(rows)
    if tuple(next(iterator)[:2]) != ("GEO Sample link", "GEO Sample code"):
        raise ValueError("Unexpected link workbook")
    groups, labels = {}, set()
    for label, *members in iterator:
        if not label or label in labels:
            raise ValueError("Duplicate or missing link-group label")
        labels.add(label)
        for member in (v for v in members if v is not None):
            if member in groups:
                raise ValueError("Sample occurs more than once in link workbook")
            groups[member] = label
    if set(groups) != set(titles):
        raise ValueError("Incomplete linkage between GEO workbook and brief")
    return groups


def aggregate(rows, groups):
    counts = defaultdict(Counter)
    by_arm = defaultdict(lambda: defaultdict(set))
    phase_by_group = defaultdict(lambda: defaultdict(set))
    phase_counts = defaultdict(lambda: defaultdict(Counter))
    exclusions = Counter()
    for title, row in rows.items():
        role = ROLES.get(row["description"])
        if role is None:
            continue
        arm, expected = role
        if title.endswith(" Replicate") or row["tissue"].endswith(" Replicate"):
            exclusions["marked_replicate"] += 1
            continue
        if row["tissue"] != expected:
            raise ValueError("Conflicting sample description and tissue")
        stage = row.get("disease stage")
        valid = stage == "Healthy" if arm.startswith("control_") else stage in {"1", "2", "3", "4"}
        if not valid:
            exclusions["missing_or_discordant_stage"] += 1
            continue
        platform = row["platform_id"]
        label = groups[title]
        phase = row.get("cycle phase", "missing")
        counts[platform][arm] += 1
        by_arm[platform][arm].add(label)
        phase_counts[platform][arm][phase] += 1
        phase_by_group[platform, label][arm].add(phase)
    result = {}
    for platform, arms in sorted(by_arm.items()):
        cp, ce = arms["control_peritoneum"], arms["control_endometrium"]
        concordance = Counter()
        for label in cp & ce:
            phases = phase_by_group[platform, label]
            p, e = phases["control_peritoneum"], phases["control_endometrium"]
            phase = next(iter(p)) if len(p) == 1 and p == e else "discordant_or_multiple"
            concordance[phase] += 1
        result[platform] = {
            "arrays": {arm: counts[platform][arm] for arm in ARMS},
            "link_groups": {arm: len(arms[arm]) for arm in ARMS},
            "link_groups_with_lesion_and_case_peritoneum": len(arms["lesion"] & arms["case_peritoneum"]),
            "link_groups_with_both_control_tissues_by_phase": dict(concordance),
            "array_cycle_category_by_arm": dict(phase_counts[platform]),
        }
    return {"geo_arrays_all_tissues": len(rows), "geo_link_groups_all_tissues": len(set(groups.values())),
            "exclusions_within_candidate_descriptions": dict(exclusions), "by_platform": result}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    directory = parser.parse_args().directory
    rows = records(pinned(directory / "geo-brief.soft", BRIEF).decode("utf-8"))
    pinned(directory / "geo-link.xlsx", LINK)
    book = load_workbook(directory / "geo-link.xlsx", read_only=True)
    try:
        groups = group_map(book.active.values, rows)
    finally:
        book.close()
    result = aggregate(rows, groups)
    mapping = load_wg6_annotation(directory)
    panel = Path(__file__).resolve().parents[1] / "evidence/lesion-transfer/panel.json"
    if sha256(panel.read_bytes()).hexdigest() != PANEL:
        raise ValueError("Unexpected frozen panel")
    result["wg6_human_refseq_addresses_for_frozen_entrez"] = {
        str(g["entrez"]): sorted(address for address, gene in mapping.items() if gene == g["entrez"])
        for g in json.loads(panel.read_text())
    }
    result["wg6_manifest_address_count"] = len(mapping)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
