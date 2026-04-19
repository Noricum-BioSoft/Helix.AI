# How `handle_natural_command` Treats Multi-Step Pipelines

## What you’re asking

When the router returns `tool = "handle_natural_command"`, the request is **not** mapped to a single predefined tool. You want to know:

1. How multi-step pipelines are handled (some steps = toolbox tools, others = not in toolbox).
2. When and how code is generated for steps that don’t have a tool.
3. Whether the system can mix “call tool X” and “generate code for this step” in one workflow.

Yes: we support that mix. Below is where it comes from and where the gaps are.

---

## What `handle_natural_command` means

- **Router contract:** The router must return exactly one `tool` from a fixed list. For composite or ambiguous requests it returns `"handle_natural_command"` to mean: “no single tool; let the **agent/planner/tool-generator** handle this.”
- **Downstream:** That one “tool” is not in `safe_tools`, so the agent doesn’t short-circuit to a single tool. Execution goes through the **full agent path** (and/or the **plan execution** path when a multi-step plan is built). So “handle_natural_command” is the entry point for multi-step and mixed toolbox/code-gen behavior.

---

## Two ways multi-step can happen

### 1. Explicit workflow (user says “A then B”)

- **Detection:** In `main.py`, `_looks_like_workflow(command)` is true when the command contains delimiters like `" and then "`, `" then "`, `"->"`, `"\n"`, `";"` (and length > 20).
- **Plan building:** `command_router.route_plan(command, session_context)` splits the command on those delimiters and calls `route_command(part)` for each part. So each step gets a **tool name** (or `handle_natural_command`).
- **Execution:** The broker runs `__plan__` with that plan. For each step it calls `_tool_executor(step.tool_name, step.arguments)`:
  - If `step.tool_name` is a **toolbox tool** (e.g. `sequence_alignment`, `fastqc_quality_analysis`), that tool is invoked directly.
  - If `step.tool_name` is **`handle_natural_command`**, the executor calls the **agent** for that step (e.g. `handle_command(step.arguments["command"], ...)`). The agent can then map to a tool or trigger the **tool generator** for that sub-command (e.g. “calculate consensus from the alignment”).
- So in one plan you can have:
  - Step 1: `sequence_alignment` (toolbox).
  - Step 2: `handle_natural_command` with `command = "calculate consensus from the alignment"` → agent/tool-generator (possibly generated code).

### 2. Single composite sentence (no “then”)

- Example: “perform a multi sequence alignment and calculate the consensus sequence of the following sequences: …”
- **No workflow split:** There are no “then”/newline/semicolon delimiters, so `_looks_like_workflow` is false. We do **not** build a plan from the router; we only call `route_command(whole_command)` once and get `handle_natural_command`.
- **Single agent invocation:** The whole command is sent to the agent via `handle_command(full_command, ...)`. So **one** agent call must handle both “align” and “consensus.”
- **Inside the agent:** The execute path can call `_handle_multi_step_workflow`, which uses:
  - **workflow_planner_agent** → builds a `WorkflowPlan` (list of operations with `tool_name`, etc.),
  - **workflow_executor** → for each step:
    - **If the step’s tool is in the toolbox** → `_execute_tool(tool_name, params)`.
    - **If the step’s tool is not available** → `_generate_and_execute_tool(...)` (generate code and run it).
- So even for a single sentence, the agent can **internally** plan a multi-step workflow and mix toolbox tools and generated code; the split is driven by the **workflow planner**, not by the router’s “then” split.

---

## Where “tool in toolbox” vs “generate code” is decided

- **Plan-based path** (explicit “A then B”):
  - Each step’s `tool_name` comes from **router** (`route_plan` → `route_command(part)`). So a step can be `sequence_alignment` (toolbox) or `handle_natural_command` (agent/code-gen for that step).
- **Agent workflow path** (single sentence → `_handle_multi_step_workflow`):
  - **WorkflowExecutor** (e.g. in `workflow_executor.py`) does:
    - `_check_tool_availability(operation.tool_name)` (e.g. via `agent_tools` / tool registry).
    - If **available** → `_execute_tool(tool_name, params)`.
    - If **not available** → `_generate_and_execute_tool(operation, parameters, ...)` (tool generator produces and runs code for that step).
  - So the same workflow can have some steps as toolbox tools and some as generated code.

---

## Summary

| Scenario | Who builds steps? | Toolbox vs code-gen |
|----------|-------------------|----------------------|
| User says “A then B” | Router (`route_plan` splits; each part → `route_command` → one tool or `handle_natural_command`) | Steps with a toolbox tool → run that tool; steps with `handle_natural_command` → agent (and possibly tool generator) for that step. |
| User says “A and B” (one sentence) | No router plan. Agent gets full command. | Agent’s workflow planner + executor: steps with an available tool → `_execute_tool`; steps with no tool → `_generate_and_execute_tool`. |

So:

- **`handle_natural_command`** is the label that means “don’t use a single router tool; use the agent (and possibly multi-step planning and code generation).”
- **Multi-step** can mix toolbox tools and non-toolbox steps: either via an explicit plan (steps from router, some of which are `handle_natural_command`) or via the agent’s internal workflow (WorkflowExecutor chooses execute vs generate per step).
- **Code generation** is used when a step’s tool is not in the toolbox: WorkflowExecutor calls `_generate_and_execute_tool` for that step.

---

## Possible improvement: use router’s `suggested_steps` for planning

Right now, when the router returns `handle_natural_command` it can also return **`router_reasoning.suggested_steps`** (e.g. `["Align sequences", "Compute consensus"]`). That array is **not** yet used to build a Plan when the user didn’t say “then.”

A possible enhancement: when `tool === "handle_natural_command"` and `router_reasoning.suggested_steps` has multiple items, the backend could **synthesize a Plan** from those steps (e.g. one step per suggested_steps entry, with descriptions) and run it through the same plan executor. That would give you explicit multi-step execution and mixed toolbox/code-gen even for single-sentence composite requests, without requiring “then” in the user text.
