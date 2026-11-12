# Research data policy

This directory contains **documentation only**. It must never store real patient-derived counts, raw genetic data, individual-level metadata, confidential study files or consent materials in the public repository.

Use local data/raw, data/processed or approved controlled-computing storage. These paths are gitignored, but permissions and access agreements remain mandatory. Synthetic examples are generated under results/ by scripts/run_demo.sh and are also gitignored.

See docs/DATA_SOURCES.md and docs/SCHEMA.md. Run the local checksum utility on approved files rather than sharing any sensitive bytes in GitHub.
