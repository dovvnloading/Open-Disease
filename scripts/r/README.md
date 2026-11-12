# Offline R ingestion utilities

- pseudobulk_core.R — independently authored, sparse-aware donor × cell-type count export. Requires Matrix; uses strictly synthetic smoke tests in CI.
- export_seurat_pseudobulk.R — Seurat-object wrapper requiring an **explicit, locally reviewed sample-to-donor map** and local raw count object.

Read docs/ATLAS_EXPORT.md before any real data processing. NO dataset is downloaded by these scripts. They are not scientific validation of endometriosis mechanisms and must not be used on patient data outside an approved setting.

The repository also provides a different **cell-level mapping** exporter (scripts/seurat_to_pseudobulk.R) and a synthetic cross-check of both approaches in tests/r/test_cross_exporters.R. See docs/ATLAS_EXPORT.md. Avoid inventing 466,371 patient mappings: every source participant/sample assignment must be independently verified.
