# Helix.AI — Differentiation Strategy

## Problem Statement

Bioinformaticians can already use general-purpose AI tools (Claude Code, Cursor, ChatGPT) to
generate data analysis code. Helix.AI must offer something fundamentally different — not better
code generation, but **autonomous biological data analysis** that requires no coding at all.

---

## Core Value Proposition

> Upload your data. Describe your goal. Helix proposes a plan, you approve it, Helix executes it
> and delivers interpreted biological findings — no code, no environment setup, no debugging.

---

## Competitive Differentiation

| Dimension | Claude Code / Cursor | **Helix.AI (target)** |
|---|---|---|
| **Output** | Code | Results |
| **Who runs it** | The user | Helix (autonomously) |
| **Domain knowledge** | Generic | Bioinformatics-native |
| **Data understanding** | None | Upload-time profiling, format awareness |
| **Requires coding?** | Yes | No |
| **Reproducibility** | Manual | Built-in provenance & audit trail |
| **Iteration** | User rewrites prompt | Helix revises its own plan |
| **Result delivery** | Raw stdout / code diff | Interpreted biological findings |

---

## The Autonomous Analysis Loop

The differentiating user experience is a **Plan → Approve → Execute** loop:

```
1. User uploads data      →  Helix profiles it (format, schema, quality)
2. User states goal       →  Helix proposes a biological analysis plan
3. User approves plan     →  Helix executes autonomously
4. Helix delivers results →  Interpreted findings, not raw output
5. User gives feedback    →  Helix revises and re-runs
```

No coding required at any step. The user is a scientist, not an engineer.

---

## Current State vs Target

### Done (foundation)
- File upload + profiling (tabular, FASTA, VCF, BED, single-cell, alignment)
- Code interpreter for tabular QA
- Session & history management
- Tool routing infrastructure

### Critical Gaps

| # | Gap | Why it matters |
|---|---|---|
| 1 | **Plan → Approve → Execute loop** | The defining UX that separates Helix from code assistants |
| 2 | **Domain-aware planner** | Agent must know what analysis is appropriate per data type (VCF ≠ h5ad ≠ clinical table) |
| 3 | **Result interpretation** | Raw plots and CSVs are not findings; Helix must summarize what results mean biologically |
| 4 | **Reproducibility** | Every analysis must be replayable on new data without re-prompting |
| 5 | **Integrated tool suite** | BLAST, pathway enrichment, differential expression, variant annotation as first-class tools |

---

## Anchor Use Case (MVP Demo Target)

> **"Find genes that correlate with treatment response"**
>
> 1. User uploads `clinical_trial_data.xlsx`
> 2. User types: *"find which genes correlate with treatment response"*
> 3. Helix proposes plan: load data → filter by response column → compute correlations →
>    plot top hits → summarize findings
> 4. User approves
> 5. Helix executes and returns: a ranked gene list, a correlation heatmap, and a one-paragraph
>    biological summary

This single loop, executed reliably, is the demo that differentiates Helix from any
general-purpose AI tool.

---

## Next Steps

1. Spec and implement the **Plan → Approve → Execute loop** for the anchor use case
2. Build the **domain-aware planner** that maps file type + user intent → analysis plan
3. Add **result interpretation** layer (LLM-generated biological summary of outputs)
4. Define reproducibility contract (workflow ID, replayable execution, audit log)
