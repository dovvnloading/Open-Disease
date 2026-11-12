"""Simulated toy data: this does not constitute human or biological evidence."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd


def generate_synthetic(out_dir: str | Path, *, seed: int = 1337) -> tuple[Path, Path]:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    genes = [f"SYN{i:03d}" for i in range(50)]
    planted = {"SYN000": 5.0, "SYN001": 3.5, "SYN002": 0.24}
    rows = []
    for study in ("synthetic-discovery", "synthetic-holdout"):
        for phase in ("proliferative", "secretory"):
            for diagnosis in ("case", "control"):
                for index in range(6):
                    donor = f"{study}-{phase}-{diagnosis}-{index:02d}"
                    donor_factor = float(rng.lognormal(0, 0.11))
                    n_cells = int(rng.integers(35, 60))
                    counts = []
                    for gene in genes:
                        factor = planted.get(gene, 1) if diagnosis == "case" else 1
                        offset = 1.3 if phase == "secretory" else 1
                        counts.append(int(rng.poisson(90 * donor_factor * offset * factor)))
                    library_size = sum(counts)
                    for gene, count in zip(genes, counts):
                        rows.append(dict(study_id=study, donor_id=donor, diagnosis=diagnosis,
                                         phase=phase, hormone_treatment="none",
                                         cell_type="epithelial", gene_id=gene, counts=count,
                                         library_size=library_size, n_cells=n_cells))
    pseudobulk = out_dir / "pseudobulk.csv"
    pd.DataFrame(rows).to_csv(pseudobulk, index=False)
    evidence = out_dir / "synthetic_evidence.csv"
    pd.DataFrame([
        dict(gene_id=gene, evidence_type=typ, support_score=score,
             source_id="SYNTHETIC_NOT_REAL_GENETICS", safety_review="not_assessed")
        for gene, typ, score in [
            ("SYN000", "fine_mapping", 0.9), ("SYN001", "colocalization", 0.8),
            ("SYN002", "integrated_multiomics", 0.7), ("SYN003", "locus_only", 0.95),
            ("SYN020", "fine_mapping", 0.7)
        ]
    ]).to_csv(evidence, index=False)
    (out_dir / "SYNTHETIC_ONLY.json").write_text(json.dumps({
        "seed": seed, "data": "100% SIMULATED, NO HUMAN PARTICIPANTS",
        "planted_case_multipliers": planted, "not_scientific_evidence": True,
        "note": "The two toy cohorts share a programmed effect; not true independent validation."
    }, indent=2) + "\n", encoding="utf-8")
    return pseudobulk, evidence
