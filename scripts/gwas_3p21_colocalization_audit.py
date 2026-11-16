"""Verify coordinate coverage and overlap of published 3p21.31 gene lists offline."""

import argparse
import hashlib
import json
from pathlib import Path


EVIDENCE = Path(__file__).resolve().parents[1] / "evidence"
OBSERVATIONS = EVIDENCE / "gwas-3p21-colocalization-audit" / "observations.json"
TWAS_ARCHIVE = EVIDENCE / "uterine-genetic-direction-2026-11-15" / "results.json"
GRCH37_CHR3_LENGTH = 198022430  # RefSeq NC_000003.11


def positive_integer(value, field):
    if type(value) is not int or value < 1:
        raise ValueError(f"{field} must be a positive integer")
    return value


def gene_symbols(values):
    if not isinstance(values, list) or not all(
        isinstance(symbol, str) and symbol and symbol.strip() == symbol for symbol in values
    ):
        raise ValueError("Gene lists must contain nonempty gene symbols")
    return set(values)


def verify(observations, archive_bytes):
    """Return descriptive comparisons using 1-based inclusive GRCh37 coordinates."""
    digest = hashlib.sha256(archive_bytes).hexdigest()
    if digest != observations["uterus_twas_source_sha256"]:
        raise ValueError("Koller uterus TWAS archive SHA-256 differs from the source pin")

    lead = observations["disease_lead"]
    candidate = observations["coloc_candidate"]
    for marker in (lead, candidate):
        if (marker["assembly"], marker["chromosome"]) != ("GRCh37", "3"):
            raise ValueError("Both markers must use GRCh37 chromosome 3")
        if positive_integer(marker["position"], marker["id"]) > GRCH37_CHR3_LENGTH:
            raise ValueError(f"{marker['id']} is outside GRCh37 chromosome 3")

    distance = abs(candidate["position"] - lead["position"])
    windows = {}
    for label in ("published", "audit"):
        radius = positive_integer(observations[f"{label}_radius_bp"], f"{label} radius")
        start = max(1, lead["position"] - radius)
        end = min(GRCH37_CHR3_LENGTH, lead["position"] + radius)
        windows[label] = {
            "radius_bp": radius,
            "start_inclusive": start,
            "end_inclusive": end,
            "contains_candidate": start <= candidate["position"] <= end,
            "distance_outside_bp": max(0, start - candidate["position"], candidate["position"] - end),
        }

    selected = json.loads(archive_bytes)["selected"]
    twas = gene_symbols([row["symbol"] for row in selected if row["band"] == "3p21.31"])
    eqtl = gene_symbols(observations["reported_endometrial_eqtl_genes"])
    return {
        "assembly": "GRCh37",
        "chromosome": "3",
        "window_center": lead["id"],
        "window_center_position": lead["position"],
        "coloc_candidate": candidate["id"],
        "coloc_candidate_position": candidate["position"],
        "physical_distance_bp": distance,
        "windows": windows,
        "published_gene_annotations": {
            "reported_mtag_endometrial_eqtl": sorted(eqtl),
            "selected_koller_uterus_twas": sorted(twas),
            "both": sorted(eqtl & twas),
            "mtag_endometrial_eqtl_only": sorted(eqtl - twas),
            "koller_uterus_twas_only": sorted(twas - eqtl),
        },
        "uterus_twas_source_sha256": digest,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", type=Path, default=OBSERVATIONS)
    parser.add_argument("--twas-archive", type=Path, default=TWAS_ARCHIVE)
    args = parser.parse_args()
    observations = json.loads(args.observations.read_text(encoding="utf-8"))
    print(json.dumps(verify(observations, args.twas_archive.read_bytes()), indent=2))


if __name__ == "__main__":
    main()
