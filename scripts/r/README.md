# Offline R ingestion utilities

- pseudobulk_core.R — independently authored, sparse-aware donor × cell-type count export. Requires Matrix; uses strictly synthetic smoke tests in CI.
- export_seurat_pseudobulk.R — Seurat-object wrapper requiring an **explicit, locally reviewed sample-to-donor map** and local raw count object.

Read docs/ATLAS_EXPORT.md before any real data processing. NO dataset is downloaded by these scripts. They are not scientific validation of endometriosis mechanisms and must not be used on patient data outside an approved setting.
