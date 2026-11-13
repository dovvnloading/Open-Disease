"""Commands generate research artifacts only; NOT clinical software."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import pandas as pd
from .analysis import analyze_pseudobulk
from .geo_microarray import download_series_matrix, prepare_series_matrix
from .pseudobulk import aggregate_cells
from .synthetic import generate_synthetic
from .targets import rank_targets, validate_targets


def write_csv(frame: pd.DataFrame, path: str) -> None:
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(dest, index=False)
    print(f"Wrote {len(frame)} rows to {dest}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Open-Disease nonclinical research pipeline")
    commands = parser.add_subparsers(dest="command", required=True)
    s = commands.add_parser("simulate", help="Create synthetic donor-level demo")
    s.add_argument("--out-dir", required=True)
    s.add_argument("--seed", type=int, default=1337)
    a = commands.add_parser("aggregate", help="Small-table cell-to-donor aggregation")
    a.add_argument("--cells", required=True)
    a.add_argument("--out", required=True)
    a = commands.add_parser("analyze", help="Donor-level exploratory differential expression")
    a.add_argument("--pseudobulk", required=True)
    a.add_argument("--out", required=True)
    a.add_argument("--min-donors", type=int, default=3)
    a = commands.add_parser("rank", help="Rank pre-specified mapped candidates in discovery")
    a.add_argument("--de", required=True)
    a.add_argument("--evidence", required=True)
    a.add_argument("--discovery-study", required=True)
    a.add_argument("--phase", required=True)
    a.add_argument("--cell-type", required=True)
    a.add_argument("--fdr", type=float, default=0.10)
    a.add_argument("--top-n", type=int, default=25)
    a.add_argument("--out", required=True)
    a = commands.add_parser("validate", help="Test frozen candidate list in separate study")
    a.add_argument("--ranked", required=True)
    a.add_argument("--de", required=True)
    a.add_argument("--validation-study", required=True)
    a.add_argument("--fdr", type=float, default=0.10)
    a.add_argument("--out", required=True)
    a = commands.add_parser("sha256", help="Compute local input checksum")
    a.add_argument("files", nargs="+")
    g = commands.add_parser("geo-download", help="Opt-in public GSE51981 Series Matrix download (local only)")
    g.add_argument("--out", required=True, help="Ignored local path, e.g. data/raw/GSE51981.txt.gz")
    g = commands.add_parser("geo-prepare", help="Offline GEO probe table + UNREVIEWED metadata manifest")
    g.add_argument("--matrix", required=True)
    g.add_argument("--manifest", required=True)
    g.add_argument("--expression", required=True)
    g.add_argument("--audit", required=True)
    g.add_argument("--expected-samples", type=int, default=148)
    args = parser.parse_args(argv)
    try:
        if args.command == "simulate":
            paths = generate_synthetic(args.out_dir, seed=args.seed)
            print("SIMULATED ONLY: no patient data or disease discoveries")
            for path in paths:
                print(path)
        elif args.command == "aggregate":
            write_csv(aggregate_cells(pd.read_csv(args.cells)), args.out)
        elif args.command == "analyze":
            write_csv(analyze_pseudobulk(pd.read_csv(args.pseudobulk),
                                         min_donors_per_group=args.min_donors), args.out)
        elif args.command == "rank":
            write_csv(rank_targets(pd.read_csv(args.de), pd.read_csv(args.evidence),
                                   discovery_study=args.discovery_study, phase=args.phase,
                                   cell_type=args.cell_type, fdr=args.fdr, top_n=args.top_n), args.out)
        elif args.command == "validate":
            write_csv(validate_targets(pd.read_csv(args.ranked), pd.read_csv(args.de),
                                       validation_study=args.validation_study, fdr=args.fdr), args.out)
        elif args.command == "sha256":
            for filename in args.files:
                digest = hashlib.sha256()
                with Path(filename).open("rb") as stream:
                    for chunk in iter(lambda: stream.read(1024 * 1024), b""):
                        digest.update(chunk)
                print(digest.hexdigest(), filename)
        elif args.command == "geo-download":
            print(json.dumps(download_series_matrix(args.out), indent=2, sort_keys=True))
        elif args.command == "geo-prepare":
            print(json.dumps(prepare_series_matrix(
                args.matrix, manifest=args.manifest, expression=args.expression,
                audit=args.audit, expected_samples=args.expected_samples
            ), indent=2, sort_keys=True))
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
