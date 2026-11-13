# Contributing

Open-Disease is an **independent, agent-first computational research effort**. Contributions can be made by an AI agent or a person on ordinary hardware. No grant, institution, lab, specialist recruitment or outreach is expected.

## Choose work

1. Read `AGENTS.md`, `docs/AGENT_WORKFLOW.md`, and the relevant scientific protocol.
2. Run `python scripts/research_queue.py ready` and pick an eligible JSON task from `research/tasks/`. Check its linked issue for concurrent work; use a separate branch/PR.
3. Prefer a source-linked audit, reproducible public-data result (including null/negative), or a targeted engineering fix over speculative architecture.
4. Attach a public-safe handoff: question and task ID, exact commands, code and package versions, source URLs, provenance/checksums **only when actually obtained**, sample/exclusion assumptions, test status, negatives, uncertainty and remaining blockers.
5. Run `python scripts/research_queue.py validate` and `python -m pytest -q`. Add synthetic regression tests for code changes. Do not claim unrun tests passed.

## Hard limits

Do not invent phenotypes, donor independence, licensing permission, gene mappings or scientific results. Do not commit raw/sample-level biomedical data, private clinical records, secrets, restricted assets, or unlicensed upstream code. Keep raw and derived human-data files in ignored local directories. Do not publish patient-level outputs in PRs or CI. No diagnosis, treatment advice, therapeutic claim or autonomous patient contact.

No expert review or lab partnership is required to submit **clearly labeled exploratory computation**. Stronger causal, independent-validation, therapeutic or clinical claims require the relevant scientific evidence and, where appropriate, qualified and ethical review. Do not contact anyone without explicit operator authorization. A blocker is a legitimate output; do not quietly weaken acceptance checks to get a green result.
