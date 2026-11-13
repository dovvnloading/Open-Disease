# Open-Disease: instructions for research and coding agents

Open-Disease is a small, independent, **AI-first computational research initiative**. It is not a funded lab, clinical project, validated scientific program, or a promise of a cure. The operator supplies direction; agents investigate public evidence, write code, spend bounded compute, test and leave reproducible work. **No recruitment, fundraising, grant seeking, or contacting scientists by default.**

## Start here

1. Read `README.md`, `ETHICS_AND_LIMITATIONS.md`, and the relevant disease protocol.
2. Run `python scripts/research_queue.py validate` and `python scripts/research_queue.py ready`.
3. Select one eligible JSON task in `research/tasks/`. Read its question, sources, acceptance checks, dependencies, and compute/network budget. The queue does not automatically launch jobs or grant network consent.
4. Work on a branch/PR referencing its issue and task ID. Record commands, versions, seed, code commit, input source/checksum if acquired, elapsed time, output status, failures, and uncertainty. `python scripts/research_queue.py new-run --task TASK_ID --agent AGENT_ID --out results/TASK_ID/run.json` creates a **local ignored** provenance stub.
5. Submit only source-linked methods, safe aggregate summaries, tests and negative findings. Mark a task `done` only when its acceptance checks actually pass and a reviewable report/PR exists. Never claim to have run code or accessed data without proof.

## Scientific honesty and privacy

- Label work explicitly: **software-only**, **literature/metadata assessment**, **exploratory human-data computation**, **independent computational replication**, or **experimentally validated**. Only report the level demonstrated. Computation alone cannot establish treatment efficacy, a clinical test, or causation.
- Prefer falsifiable narrow questions and existing public data over new frameworks or speculative drug prediction. Preserve failed hypotheses and counterevidence. Cite primary sources and distinguish their claims from your own observations.
- Do not hallucinate study metadata, diagnoses, donors, permissions, effect sizes, results, completed executions, or citations. If a required mapping is ambiguous, record the blocker and stop that comparison.
- Keep person/sample-level biomedical data, raw and derived expression tables, sensitive records, access tokens, and run artifacts out of Git, issues, logs and public AI prompts. Use ignored `data/raw/`, `data/processed/`, `data/private/`, and `results/` locally. Never acquire restricted data or upload sensitive data to a model provider.
- Respect data/code licenses and access controls. Publicly visible does not imply free redistribution. Never attempt participant matching or reidentification.
- For RNA counts, donors—not cells—are the independence unit. For GSE51981, use a log2 microarray/limma path; do **not** input microarray intensities to edgeR. Preserve the other-pathology comparison as a separate group.

## Compute and autonomy

- Default to offline tests and metadata review. Public downloads occur **only if explicitly allowed by the chosen task and operator/environment**, with source checks. No unbounded crawling, self-replicating agents, unrestricted cloud spend, or implicit high-memory workloads.
- Treat task `budget` values as operating limits to respect and report, **not enforced sandbox quotas**. Pause and request a change to the task before going beyond them. Avoid repeating expensive work when a verifiable result already exists.
- Agents may propose or execute suitable computational studies and improvements without a grant, collaborator, or scientific sign-off **provided all resulting claims remain clearly exploratory and nonclinical**.
- External contact is **not an agent action**. Escalate to the operator only for a concrete exceptional need: incident/security reporting, restricted-data/legal consent, proposed human-subject or wet-lab activity, or promotion of findings to clinically consequential claims. Obtain explicit approval before any outreach. Missing independent validation is a finding, not a reason to send unsolicited messages.

## Definition of a useful result

A bounded task with a reproducible method, provenance, one of positive/null/negative/blocked outcomes, relevant controls, and an honest account of what could not be established. Prefer a defensible negative result over a persuasive invented one. See `docs/AGENT_WORKFLOW.md` for handoff expectations.
