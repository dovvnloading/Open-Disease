# Metadata-only source registry

catalog/sources.json records public evidence, access status and licensing questions.
No patient-level, pseudonymized individual counts, genomes, restricted files, tokens or fabricated checksums belong here.

    python scripts/verify_catalog.py

The validator performs **local structural checks only**. It does not download or verify any human data. For consent, reuse, conflict and statistical constraints see docs/DATA_SOURCES.md and docs/PROTOCOL.md.
