#!/usr/bin/env bash
set -euo pipefail
# Toy results are SIMULATED and are NOT endometriosis findings.
OUT="results/synthetic-demo"
if [ "$#" -ge 1 ]; then OUT="$1"; fi
mkdir -p "$OUT"
python -m open_disease simulate --out-dir "$OUT"
python -m open_disease analyze --pseudobulk "$OUT/pseudobulk.csv" --out "$OUT/de.csv"
python -m open_disease rank --de "$OUT/de.csv" --evidence "$OUT/synthetic_evidence.csv" \
  --discovery-study synthetic-discovery --phase proliferative \
  --cell-type epithelial --out "$OUT/ranked.csv"
python -m open_disease validate --ranked "$OUT/ranked.csv" --de "$OUT/de.csv" \
  --validation-study synthetic-holdout --out "$OUT/validation.csv"
echo "SIMULATED ONLY: $OUT/validation.csv"
