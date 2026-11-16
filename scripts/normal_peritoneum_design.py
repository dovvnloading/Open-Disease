"""Select native noncase pairs, then lock physical-chip eligibility without signal."""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass, replace
import json

from audit_endometdb import PHASES
from audit_peritoneal_controls import ARMS, UTERUS, native_assays

PERITONEUM = ARMS[2]
CHIP_STATUSES = {"WG6", "HT12_R1", "HT12_R2", "missing", "repeated_addresses",
                 "no_exact_manifest", "ambiguous_manifest"}


@dataclass(frozen=True)
class Pair:
    donor: str
    run: str
    peritoneum_specimen: str
    uterus_specimen: str
    peritoneum_assays: tuple[str, ...]
    uterus_assays: tuple[str, ...]
    uterine_phase: str | None
    medication: str

    @property
    def sensitivity_phase(self):
        if self.medication == "stored_false" and self.uterine_phase in PHASES.values():
            return self.uterine_phase
        return None


def frozen_pairs_json(pairs):
    """Private local freeze only; contains native donor, specimen and assay keys."""
    return json.dumps([asdict(p) for p in sorted(pairs, key=lambda p: p.donor)],
                      indent=2, sort_keys=True) + "\n"


def phase_counts(pairs):
    return dict(sorted(Counter(p.sensitivity_phase for p in pairs if p.sensitivity_phase).items()))


def candidate_pairs(tables):
    assays, medication, _ = native_assays(tables)
    by_donor = defaultdict(lambda: defaultdict(list))
    for a in assays:
        if a.arm in (PERITONEUM, UTERUS):
            by_donor[a.donor][a.arm].append(a)

    excluded = Counter()
    candidates = []
    for donor, arms in sorted(by_donor.items()):
        missing = {PERITONEUM, UTERUS} - arms.keys()
        if missing:
            excluded["missing_" + missing.pop()] += 1
            continue
        p, u = arms[PERITONEUM], arms[UTERUS]
        shared = {a.run for a in p} & {a.run for a in u}
        if len(shared) != 1:
            reason = "no_shared_run" if not shared else "multiple_shared_runs"
            excluded[reason] += 1
            continue
        run = shared.pop()
        p = [a for a in p if a.run == run]
        u = [a for a in u if a.run == run]
        if any(len({a.specimen for a in group}) != 1 for group in (p, u)):
            excluded["multiple_specimens_in_shared_run"] += 1
            continue
        phases = {a.endometrial_histology for a in u}
        if len(phases) != 1:
            raise ValueError("Conflicting linked histology for a single selected uterine specimen")
        histology = phases.pop()
        candidates.append(Pair(
            donor, run, p[0].specimen, u[0].specimen,
            tuple(sorted(a.assay_id for a in p)), tuple(sorted(a.assay_id for a in u)),
            histology.lower() if histology else None, medication[donor],
        ))

    report = {
        "healthy_donors_with_either_tissue": len(by_donor),
        "healthy_donors_with_both_tissues": sum(len(arms) == 2 for arms in by_donor.values()),
        "candidate_pairs": len(candidates),
        "excluded_donors_by_design": dict(sorted(excluded.items())),
        "candidate_assays_by_tissue": {
            "peritoneum": sum(len(p.peritoneum_assays) for p in candidates),
            "uterus": sum(len(p.uterus_assays) for p in candidates),
        },
        "candidate_sensitivity_phase_counts": phase_counts(candidates),
    }
    return candidates, report


def verified_pairs(candidates, statuses):
    """Keep every verified assay of each selected specimen; drop pairs missing an arm."""
    assigned = [a for p in candidates for a in p.peritoneum_assays + p.uterus_assays]
    if len(set(assigned)) != len(assigned) or len({p.donor for p in candidates}) != len(candidates):
        raise ValueError("Repeated donor or assay in private candidate design")
    if set(statuses) != set(assigned) or set(statuses.values()) - CHIP_STATUSES:
        raise ValueError("Chip census must classify exactly the candidate assays with known statuses")
    retained = []
    excluded = Counter()
    classification = {"peritoneum": Counter(), "uterus": Counter()}
    for p in candidates:
        classification["peritoneum"].update(statuses[a] for a in p.peritoneum_assays)
        classification["uterus"].update(statuses[a] for a in p.uterus_assays)
        peritoneum = tuple(a for a in p.peritoneum_assays if statuses[a] == "WG6")
        uterus = tuple(a for a in p.uterus_assays if statuses[a] == "WG6")
        if peritoneum and uterus:
            retained.append(replace(p, peritoneum_assays=peritoneum, uterus_assays=uterus))
            continue
        missing_arms = []
        if not peritoneum:
            missing_arms.append("peritoneum")
        if not uterus:
            missing_arms.append("uterus")
        excluded["no_verified_" + "_or_".join(missing_arms)] += 1

    report = {
        "candidate_pairs": len(candidates), "verified_pairs": len(retained),
        "candidate_assay_status_by_tissue": {k: dict(sorted(v.items())) for k, v in classification.items()},
        "excluded_donors_by_chip": dict(sorted(excluded.items())),
        "retained_assays_by_tissue": {"peritoneum": sum(len(p.peritoneum_assays) for p in retained),
                                     "uterus": sum(len(p.uterus_assays) for p in retained)},
        "retained_medication_codes": dict(sorted(Counter(p.medication for p in retained).items())),
        "retained_uterine_histology": dict(sorted(Counter(p.uterine_phase or "missing" for p in retained).items())),
        "retained_sensitivity_phase_counts": phase_counts(retained),
        "minimum_pairs": 6, "passes_pair_count_gate": len(retained) >= 6,
    }
    return retained, report
