# Consolidation plan, question B ("Roots or beaches?")

Goal: every number measured, robust and error-free BEFORE choosing the report structure.
Spec: pipeline/config.py (frozen). Design: reports/diseno_B_una_pagina.md.

## Shared rules for every agent
- Python: `.venv/bin/python -I`. Do not commit. Do not edit or delete files in explore/ (legacy, read-only).
- Write only to your own paths (listed in your brief). Decisions go to decisions.md (append, one line each).
- Coastal flag: column coastal_true of explore/out_p86_coastal_flags.csv. Never the buffer method of p73/p92.
- Every number: script + output file that produced it. Seeds from config.

## Wave 1 (parallel)
| ID | Agent | Task | Writes to |
|---|---|---|---|
| W1-core | Kiro Opus 5.5 medium | pipeline core: networks, basic properties, rewiring, representativeness | pipeline/b01_*, b02_*, tests/, results/b01*, results/b02* |
| W1-audit | Kiro Opus 5.5 medium | audit of legacy explore/ scripts and reported numbers | research_notes/audit_B_existing_code.md |
| W1-lit | Antigravity Neon high | literature, novelty, exodus figure | research_notes/literature_B.md, research_notes/references_B.bib |
| W1-me | orchestrator | read lectures + past examples; write guide | _local/notes_course_methods.md (private) |

## Wave 2 (after W1-core)
R1 k-core promotion + gravity null; R2 roots alignment across sister nights; R3 summer territories robustness; pattern mining.

## Wave 3
Adversarial review as Prof. Shi Zhou.
