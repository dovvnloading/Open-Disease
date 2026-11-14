# Canonical source registry moved

The machine-readable public source inventory lives at [data/sources.json](../data/sources.json).
The canonical narrative review is **DATA_SOURCES.md** at the repository root.

    python scripts/verify_catalog.py

For a source used in an expression analysis, the catalog records `human_data_used_by_open_disease: true` and links to an aggregate research report. Metadata audits alone leave this field false. The input arrays and individual records stay in ignored local paths. See the [GSE25628 FOSB result](../docs/research-notes/GSE25628_FOSB_RESULTS_2026-11-13.md) and the [GSE51981 eligibility audit](../docs/research-notes/GSE51981_COHORT_AUDIT_2026-11-13.md).
