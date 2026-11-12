# Candidate evidence input

The target-ranking example intentionally uses **fake SYN### genes and fake evidence**. It must not be confused with genetic or experimental findings.

A real assessment must have exactly one evidence row per gene with columns:

```
gene_id,evidence_type,support_score,source_id,safety_review
```

- `gene_id`: stable mapped gene identifier, with genome/transcript release documented in the study protocol.
- `evidence_type`: fine_mapping, colocalization, integrated_multiomics, or locus_only. The software excludes locus_only entries; merely nearby a significant SNP is not fine mapping.
- `support_score`: a transparent **heuristic 0–1 assessment**, with a predeclared rubric and sensitivity analysis. This is not P(causal), therapeutic response probability, or an objective genetic likelihood.
- `source_id`: full DOI/accession/version/auditable reference, not a claim invented by a model.
- `safety_review`: not_assessed, needs_review or known_concern. All values require human domain review before laboratory selection.

No verified real-gene evidence file exists in the project yet. The ranking formula is a software prototype, **not registered biological inference**.
