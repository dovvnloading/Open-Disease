"""CI guard for accidentally tracked restricted data; not a privacy guarantee."""
from __future__ import annotations

import subprocess
from pathlib import PurePosixPath

FORBIDDEN_DIRS = ("data/raw/", "data/processed/", "data/private/", "results/", "outputs/")
FORBIDDEN_EXTENSIONS = (".rds", ".h5ad", ".bam", ".cram", ".fastq", ".fastq.gz",
                        ".vcf", ".vcf.gz", ".parquet", ".h5", ".pgen")


def check_paths(paths: list[str]) -> list[str]:
    bad = []
    for path in paths:
        p = path.lower().replace("\\", "/")
        if (p.startswith(FORBIDDEN_DIRS) or p.endswith(FORBIDDEN_EXTENSIONS) or
                PurePosixPath(p).name in {".env", "id_rsa", "id_ed25519"}):
            bad.append(path)
    return bad


def main() -> None:
    tracked = subprocess.check_output(["git", "ls-files", "-z"]).decode("utf-8").split("\0")
    violations = check_paths([p for p in tracked if p])
    if violations:
        raise SystemExit("Restricted paths/types tracked: " + ", ".join(violations))
    print(f"Tracked-file check passed ({len(tracked) - 1} files). Separate human-data audit required.")


if __name__ == "__main__":
    main()
