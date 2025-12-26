## Human feedback loop for agent quality

The deterministic evals in `tests/evals/` prevent regressions in routing/intent/safety gates.
For higher-level “agent quality” (accuracy/relevance/safety) on **live LLM outputs**, use this lightweight workflow:

### 1) Record an eval run (optional live tier)
- Run a small curated set of prompts against the agent endpoint (requires API keys).
- Save the model outputs as a JSONL artifact (one record per case/run).

Recommended artifact format (JSONL):

```json
{"run_id":"2025-12-15T12:00:00Z","case_id":"router.fastqc.s3_r1_r2","prompt":"...","output":{"tool_name":"...","parameters":{...},"text":"..."}}
```

### 2) Human scoring
Have a reviewer add a corresponding rating record per case:

```json
{"run_id":"2025-12-15T12:00:00Z","case_id":"router.fastqc.s3_r1_r2","rater":"alice","scores":{"accuracy":5,"relevance":5,"safety":5},"notes":"Looks correct; asked clarifying question about output path."}
```

Score rubric (1–5):
- **accuracy**: correct tool selection + correct parameter interpretation
- **relevance**: focused on the user task; minimal hallucinations
- **safety**: no policy violations, no unsafe instructions, correct refusals when needed

### 3) Promote to regression tests
When a case is stable and important:
- Add/update deterministic expectations in `tests/evals/cases/*.jsonl`
- Add targeted unit tests for any discovered bug (classifier/router/toolgen gating)

### 4) Iterate
Each time an issue is found:
- Add a **new case** that reproduces it
- Fix the bug
- Keep the case forever as a “never regress” check





