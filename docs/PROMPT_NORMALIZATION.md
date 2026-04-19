# Prompt normalization / prompt engineering

## Limitation of keyword-based routing

The command router keeps a **small** set of high-confidence keyword checks (e.g. toolbox, consensus phrases). Relying on keyword matching alone would limit us to our imagination of what scientists could ask. Prefer **LLM-based routing** (`HELIX_USE_LLM_ROUTER=1`) and extending `ROUTER_TOOLS` and the router prompt over adding more hard-coded phrases. The keyword path is a fallback for mock/CI and when the LLM is unavailable.

---

## Problem

The same user intent can produce very different behavior depending on **formatting**:

| User prompt | What happens |
|-------------|----------------|
| “Perform MSA and calculate consensus of the following sequences: \\n>seq1\\nACGT\\n…” | `_split_workflow_command` splits on **every newline** → 7 steps (intro, “>seq1”, first sequence line, “>seq2”, …). All steps route to `handle_natural_command`. Poor UX. |
| “Perform MSA of the following: >seq1 ACGT >seq2 TGCA **then** calculate the consensus” | Split only on **“ then ”** → 2 steps: `sequence_alignment` + `handle_natural_command`. Correct. |

So formatting (newlines in FASTA vs. explicit “then”) shouldn’t dictate whether we get a sensible multi-step plan.

## Approaches

### 1. Structure-aware workflow splitting (implemented)

**Idea:** When splitting the command into workflow steps, don’t treat **every** newline as a step separator. Treat newlines inside “data blocks” (e.g. FASTA) as part of the payload.

- **Option A (current change):** Only split on **explicit workflow delimiters**: `" then "`, `" and then "`, `";"`, `"->"`, `"→"`. Do **not** split on bare `\n`. So a single sentence with FASTA (with newlines) stays **one** chunk and routes once (e.g. `handle_natural_command`); the agent or a later step can still plan sub-steps internally.
- **Option B (future):** Detect FASTA (and optionally code blocks): mask internal newlines (e.g. replace with space or placeholder) before splitting, then split on newlines as well, then unmask in each part. That would allow “step one\\nstep two” to become 2 steps while “intro\\n>seq1\\nACGT” stays 1 step.

The codebase uses **Option A** in `command_router._split_workflow_command`: the regex no longer includes `\n`, so FASTA newlines no longer create extra steps.

### 2. Optional prompt-rewriting step (future)

**Idea:** After the user submits the prompt, run an analysis step (heuristic or LLM) that may **rewrite** the prompt so the rest of the pipeline interprets it better.

- **Detect composite intent:** e.g. “align and calculate consensus” (no “then”). Optionally rewrite to “Align the following sequences: &lt;FASTA&gt; then calculate the consensus” so `route_plan` produces 2 steps.
- **Use router’s `suggested_steps`:** When the LLM router returns `handle_natural_command` with `suggested_steps = ["Align sequences", "Compute consensus"]`, we could either:
  - Build a **synthetic plan** from those steps (no rewrite), or
  - **Rewrite** the user prompt to a “then”-separated form and re-run routing.

This would live in a small “prompt normalization” or “prompt engineering” layer that runs **before** (or alongside) the router, and could be gated by an env flag.

## Where it lives

- **Splitting behavior:** `backend/command_router.py` → `_split_workflow_command()`. Only explicit delimiters (`then`, `and then`, `;`, `->`, `→`) are used; newlines are not step separators.
- **Future rewrite/normalization:** Could be a function in `command_router` or a separate module called from `main` before `route_command` / `route_plan`.
