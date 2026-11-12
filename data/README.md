# Local research data directory

**No clinical or participant-level data belong in this Git repository.**

Subdirectories `data/raw/`, `data/processed/` and `data/private/` are ignored by Git. Use a suitably secured and lawful analysis environment for permitted data. The source catalog `sources.json` is a publicly shareable metadata inventory, not an automated download manifest or a grant of access.

Synthetic examples are generated into the ignored `results/` directory by `bash scripts/run_demo.sh`. No real genetic, clinical or endometriosis data are included in the demo.

The donor-level pseudobulk schema and provenance requirements are documented in [DATA_SOURCES.md](../DATA_SOURCES.md).
