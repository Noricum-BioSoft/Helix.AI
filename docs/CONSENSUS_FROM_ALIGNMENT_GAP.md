# Consensus-from-alignment: current gap and target behaviour

## What’s wrong today

For a plan “align sequences **then** calculate the consensus sequence”:

1. **Step 2 is routed to the wrong tool**  
   When `handle_command("calculate the consensus sequence")` runs, the agent doesn’t get a tool call from the LLM, falls back to the router, and the **LLM router** (or keyword router in some paths) returns **sequence_alignment**. So step 2 runs **alignment again** (often with placeholder/default sequences), not a consensus computation. You see: *“Deterministic router fallback identified tool: sequence_alignment”* for step 2.

2. **No data flow from step 1 to step 2**  
   Step 2’s arguments are only `{ "command": "calculate the consensus sequence", "session_id": "..." }`. The broker does **not** pass step 1’s result (aligned sequences) into step 2. So even if we had a consensus tool or generated code, it wouldn’t receive the alignment to derive consensus from. The plan IR has no `$ref` from step 2 to `steps.step1.result`, and we don’t inject previous step outputs into `session_context` for `handle_natural_command`.

3. **No consensus in the toolbox**  
   There is no dedicated “consensus from alignment” tool. A smart system would either call a small consensus routine (e.g. majority-rule per column, with thresholding) or use an alignment tool that can output consensus; today we only re-run alignment for step 2.

So: step 2 neither **receives** the alignment nor **performs** a real consensus step; it just runs alignment again with wrong/placeholder data.

---

## What a very smart system would do

- **Semantics:** Treat “calculate the consensus sequence” as: *derive a single sequence from the **already aligned** sequences* (e.g. majority base per column, with thresholding or ambiguity codes), not as “run alignment again”.
- **Routing:** Do **not** map “consensus” to `sequence_alignment`. Map it to something that can **compute consensus**:
  - `handle_natural_command` → tool generator (generates code that takes aligned sequences and computes consensus), or
  - A dedicated **consensus** tool that takes `aligned_sequences` (and optional threshold) and returns one consensus sequence.
- **Data flow:** Step 2 must receive step 1’s output:
  - Either the plan has step 2 arguments like `{ "command": "...", "aligned_sequences": { "$ref": "steps.step1.result.alignment" } }` and the broker resolves `$ref` before calling the tool, or
  - The broker injects previous step results into `session_context` (e.g. `session_context["previous_step_results"] = context["steps"]`) so `handle_natural_command` / tool generator can see the alignment.
- **Execution:** Either run a small consensus implementation (e.g. Biopython or custom majority-rule) or generate code that does it, using the **actual** aligned sequences from step 1.

---

## Proposed fixes (short term)

1. **Never route “consensus” to sequence_alignment**  
   - In the **router**: when the user command clearly asks for consensus (e.g. “calculate the consensus”, “compute consensus”), always return `handle_natural_command`, not `sequence_alignment` (keyword path already does this; LLM path can override LLM’s choice when the command contains consensus phrases).  
   - In the **agent** router fallback: when the command contains consensus phrases and the router returns `sequence_alignment`, **ignore** that and do not return “tool_mapped”; let the agent continue to the **tool generator** so it can generate consensus-from-alignment code.

2. **Pass step 1 result into step 2**  
   - For plan steps with `tool_name == "handle_natural_command"`, the broker (or plan builder) should add to the step’s arguments a reference to previous step outputs, e.g. `aligned_sequences: { "$ref": "steps.step1.result.alignment" }`, or the broker should inject `context["steps"]` into `session_context` when calling the tool executor so the agent/tool generator can see step 1’s result and use it for “calculate the consensus sequence”.

3. **Optional: dedicated consensus tool**  
   - Add a small tool, e.g. `consensus_from_alignment(aligned_sequences, threshold=0.5)`, that takes aligned FASTA (or list of aligned sequences) and returns one consensus sequence (e.g. majority-rule with IUPAC ambiguity when below threshold). Then the router can map “calculate the consensus” to this tool when alignment is already in context; the plan would pass step 1’s output via `$ref` into this tool’s arguments.

---

## Summary

| Issue | Current | Target |
|-------|--------|--------|
| Step 2 routing | Router/agent return `sequence_alignment` → run alignment again | Return `handle_natural_command` or a consensus tool → run consensus (or generated code) |
| Step 2 input | Only `command` + `session_id`; no alignment data | Step 2 receives step 1’s alignment (via `$ref` or injected context) |
| Consensus logic | None (alignment run twice) | Majority-rule (and optional threshold) on aligned columns, or alignment tool that outputs consensus |

Implementing (1) in the agent so that when the command clearly asks for consensus we do not accept `sequence_alignment` from the router fallback and instead proceed to the tool generator is the minimal change to stop step 2 from being “alignment again” and to allow generated consensus code (once (2) is in place so that code can receive the alignment).

---

## Fix: data flow from step 1 → step 2 (implemented)

- **ExecutionBroker._discover_inputs** now treats `session_context["previous_plan_steps"]` as an input source: if a previous step result has `alignment` (list of `{name, sequence}`) or `aligned_sequences` (string), that is written to a temporary FASTA file and added as a discovered input. The tool generator receives it via `INPUT_FILES`, so "calculate the consensus sequence" runs on the actual alignment from step 1 instead of falling back to the script demo.
- Temp FASTA files are removed after execution (success or error) in the tool generator.
