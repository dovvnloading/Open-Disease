"""Audit shared-run EndometDB case pairs using clinical rows and probe addresses only."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
import json
from pathlib import Path

from audit_endometdb import HEAD_SHA256, PHASES, TAIL_SHA256, index, pinned, read_metadata
from audit_endometdb_chip_census import PREFIX_SHA256, RAW_START, classify_addresses, raw_keys
from audit_peritoneal_controls import ARMS, native_assays
from normal_peritoneum_design import candidate_pairs as noncase_pairs, verified_pairs as verified_noncase_pairs
from normal_peritoneum_preflight import (
    PRIVATE_ROOT, RAW_COPY, RAW_FK, blocks, platform_manifests,
)

LESION, PERITONEUM = ARMS[:2]


@dataclass(frozen=True)
class CasePair:
    donor: str
    run: str
    lesion: tuple[str, ...]
    peritoneum: tuple[str, ...]
    lesion_specimen: str
    peritoneum_specimen: str


def case_pairs(assays):
    """Select on the native clinical design before assigning a chip generation."""
    by_donor = defaultdict(lambda: defaultdict(list))
    for assay in assays:
        if assay.arm in (LESION, PERITONEUM):
            by_donor[assay.donor][assay.arm].append(assay)

    selected = []
    excluded = Counter()
    for donor, arms in sorted(by_donor.items()):
        if len(arms) == 1:
            excluded["only_" + next(iter(arms))] += 1
            continue
        common = {a.run for a in arms[LESION]} & {a.run for a in arms[PERITONEUM]}
        if len(common) != 1:
            excluded["multiple_shared_runs" if common else "no_shared_run"] += 1
            continue
        run = common.pop()
        lesion, peritoneum = ([a for a in arms[arm] if a.run == run] for arm in (LESION, PERITONEUM))
        if any(len({a.specimen for a in group}) != 1 for group in (lesion, peritoneum)):
            excluded["multiple_specimens_in_shared_run"] += 1
            continue
        selected.append(CasePair(donor, run, tuple(sorted(a.assay_id for a in lesion)),
                                 tuple(sorted(a.assay_id for a in peritoneum)),
                                 lesion[0].specimen, peritoneum[0].specimen))
    return selected, {
        "donors_with_either_tissue": len(by_donor),
        "donors_with_both_tissues": sum(len(arms) == 2 for arms in by_donor.values()),
        "donors_with_one_shared_run_and_one_specimen_per_arm": len(selected),
        "excluded_donors": dict(excluded),
    }


def donor_uterine_phase(tables):
    """Describe uterine histology by donor. No operation linkage is available."""
    samples = index(tables["patient_sample"], "id")
    codes = index(tables["code_value"], "id")
    levels = defaultdict(set)
    for sample_id, row in index(tables["pathology_report_endometrium"], "patient_sample_id").items():
        sample = samples[sample_id]
        if codes[sample["code_tissue_type_id"]]["value"] == "E":
            label = codes.get(row["code_histology_class_id"], {}).get("title_en")
            levels[sample["patient_id"]].add(PHASES.get(label, "other_or_unclassified"))
    return {donor: next(iter(values)) if len(values) == 1 else "mixed_or_inconsistent"
            for donor, values in levels.items()}


def chip_gate(pairs, statuses, medication, donor_phase):
    """Both tissues must contain a complete WG-6 assay in the selected run."""
    assigned = [a for p in pairs for a in (*p.lesion, *p.peritoneum)]
    if len(assigned) != len(set(assigned)) or set(assigned) != statuses.keys():
        raise ValueError("Expected one chip classification per candidate assay")
    if len({p.donor for p in pairs}) != len(pairs):
        raise ValueError("Expected one shared-run comparison per donor")

    chip_counts = {LESION: Counter(), PERITONEUM: Counter()}
    paired_statuses = Counter()
    retained = []
    retained_assays = Counter()
    for pair in pairs:
        tissue_statuses = []
        for tissue, ids in ((LESION, pair.lesion), (PERITONEUM, pair.peritoneum)):
            labels = Counter(statuses[a] for a in ids)
            chip_counts[tissue].update(labels)
            tissue_statuses.append(labels)
        left, right = tissue_statuses
        paired_statuses[",".join(sorted(left)) + " / " + ",".join(sorted(right))] += 1
        if left["WG6"] and right["WG6"]:
            retained.append(replace(pair, lesion=tuple(a for a in pair.lesion if statuses[a] == "WG6"),
                                    peritoneum=tuple(a for a in pair.peritoneum if statuses[a] == "WG6")))
            retained_assays[LESION] += left["WG6"]
            retained_assays[PERITONEUM] += right["WG6"]

    return retained, {
        "candidate_assay_status_by_tissue": {k: dict(v) for k, v in chip_counts.items()},
        "candidate_pairs_by_chip_status_lesion_then_peritoneum": dict(paired_statuses),
        "donors_passing_wg6_both_arms": len(retained),
        "donors_excluded_by_chip": len(pairs) - len(retained),
        "retained_wg6_assays_by_tissue": dict(retained_assays),
        "retained_donors_by_run": dict(Counter(p.run for p in retained)),
        "retained_donors_by_survey_hormone_flag": dict(Counter(medication[p.donor] for p in retained)),
        "retained_uterine_histology_by_donor_unverified_as_same_operation":
            dict(Counter(donor_phase.get(p.donor, "no_uterine_report") for p in retained)),
    }


def noncase_reference(retained_cases, noncase_assays, statuses, medication, prior_controls):
    """Count donors within each run, then union across runs before counting overlap."""
    cases = defaultdict(set)
    for p in retained_cases:
        cases[p.run].add(p.donor)
    controls = defaultdict(set)
    specimens = set()
    accepted_assays = 0
    candidates = Counter()
    for a in noncase_assays:
        if a.arm != ARMS[2]:
            raise ValueError("Expected only confirmed control peritoneum")
        if a.run in cases:
            status = statuses[a.assay_id]
            candidates[status] += 1
            if status == "WG6":
                controls[a.run].add(a.donor)
                specimens.add((a.donor, a.specimen))
                accepted_assays += 1

    unique = {donor for donors in controls.values() for donor in donors}
    overlap_on_flag = Counter()
    for run, donors in cases.items():
        flags_present = {medication[d] for d in controls[run]} - {"unknown"}
        overlap_on_flag.update(medication[d] for d in donors if medication[d] in flags_present)
    return {
        "run_summary": {
            run: {"case_pairs": len(cases[run]), "noncase_peritoneum_donors": len(controls[run]),
                  "case_survey_flags": dict(Counter(medication[d] for d in cases[run])),
                  "noncase_survey_flags": dict(Counter(medication[d] for d in controls[run]))}
            for run in sorted(cases)
        },
        "candidate_noncase_assays_in_retained_case_runs_by_chip_status": dict(candidates),
        "case_pairs_in_runs_with_wg6_noncase_peritoneum": sum(len(ds) for run, ds in cases.items() if controls[run]),
        "case_pairs_with_same_run_noncase_sharing_reported_nonmissing_flag": dict(overlap_on_flag),
        "stored_false_interpretation": "Schema defaults to false. A stored false does not establish absence of exposure.",
        "distinct_noncase_donors": len(unique),
        "distinct_noncase_specimens": len(specimens),
        "retained_noncase_wg6_assays": accepted_assays,
        "noncase_donors_in_multiple_case_runs": sum(sum(d in ds for ds in controls.values()) > 1 for d in unique),
        "distinct_noncase_donors_by_survey_flag": dict(Counter(medication[d] for d in unique)),
        "prior_paired_organ_control_donors": len(prior_controls),
        "overlap_with_prior_paired_organ_control_donors": len(unique & prior_controls),
        "noncase_donors_outside_prior_paired_organ_controls": len(unique - prior_controls),
    }


def audit_with_pairs(directory):
    from inflate64 import Inflater

    directory = directory.resolve()
    if directory == PRIVATE_ROOT or not directory.is_relative_to(PRIVATE_ROOT):
        raise ValueError("Use an ignored project data/private input directory")
    head = directory / "archive-start.bin"
    header = Inflater().inflate(pinned(head, HEAD_SHA256)[48:])[:RAW_START]
    if RAW_COPY not in header or RAW_FK not in header:
        raise ValueError("Native address column order or assay foreign key changed")
    tables = read_metadata(head, directory / "archive-tail.bin")
    assays, medication, _ = native_assays(tables)
    pairs, design = case_pairs(assays)
    case_ids = [a for p in pairs for a in (*p.lesion, *p.peritoneum)]
    controls = [a for a in assays if a.arm == ARMS[2]]
    prior_candidates, _ = noncase_pairs(tables)
    prior_ids = {a for p in prior_candidates for a in (*p.peritoneum_assays, *p.uterus_assays)}
    selected = set(case_ids) | {a.assay_id for a in controls} | prior_ids
    native_ids = {r["id"].encode() for r in tables["sample_analysis_result"]}
    classifications, census = classify_addresses(raw_keys(blocks(directory)), [a.encode() for a in selected],
                                                 native_ids, platform_manifests(directory))
    statuses = {a.decode(): status for a, status in classifications.items()}
    retained_cases, gate = chip_gate(pairs, {a: statuses[a] for a in case_ids},
                                     medication, donor_uterine_phase(tables))
    prior_pairs, _ = verified_noncase_pairs(prior_candidates, {a: statuses[a] for a in prior_ids})
    report = {
        "source": "https://doi.org/10.6084/m9.figshare.12195975.v1",
        "head_sha256": HEAD_SHA256, "tail_sha256": TAIL_SHA256, "prefix_sha256": PREFIX_SHA256,
        "clinical_design": design, "address_census_case_control_and_prior_organ_candidates": census,
        "chip_gate": gate,
        "same_run_noncase_reference": noncase_reference(retained_cases, controls, statuses, medication,
                                                        {p.donor for p in prior_pairs}),
    }
    return report, retained_cases


def audit(directory):
    report, _ = audit_with_pairs(directory)
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    print(json.dumps(audit(parser.parse_args().directory), indent=2, sort_keys=True))
