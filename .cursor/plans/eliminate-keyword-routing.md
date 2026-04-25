# Plan: Eliminate All Keyword-Based Routing

**Status:** P0 complete ✅ — P1/P2/P3 pending (future PRs)  
**Branch:** `feature/p0-benchmark-gates`  
**Goal:** Replace every keyword/regex classifier with LLM calls.  No heuristic fallbacks in production code. If LLM is unavailable, the system fails with a clear error rather than silently taking the wrong path.

---

## Inventory (17 locations, prioritized)

### ✅ P0 — Complete

| # | File | Change | Status |
|---|------|--------|--------|
| 1 | `backend/orchestration/approval_policy.py` | Replaced `APPROVAL_COMMANDS` set + all keyword helpers with LLM staging classifier (`_classify_staging_intent`). | ✅ Done |
| 2 | `backend/command_router.py` | `HELIX_LLM_ROUTER_FIRST` defaults to `"1"`; raises `RoutingError` when LLM unavailable. | ✅ Done |
| 3 | `backend/intent_classifier.py` | Removed `_classify_intent_with_heuristics`; LLM failure raises `IntentClassificationError`. | ✅ Done |
| 4 | `backend/orchestration/approval_classifier.py` | Removed keyword fallback; LLM failure raises. Added `_INTERROGATIVE_PREFIXES` early-rejection and expanded `_ANALYTICAL_PATTERNS`. | ✅ Done |
| 5 | `tests/unit/backend/conftest.py` | Added autouse fixtures `_mock_staging_classifier` + `_mock_intent_classifier` that patch LLM classifiers in `HELIX_MOCK_MODE=1` — 702 unit tests pass without a live API key. | ✅ Done |

### 🟠 P1 — Architectural replacement (own PR)

| # | File | What to change |
|---|------|---------------|
| 5 | `backend/command_router.py` | Full replacement: delete all keyword branches (L280–737). `route_command` = call `_route_with_llm` only. Merge `_route_with_llm` into `route_command`. Remove `_scrub_for_keyword_matching`, `_is_workflow_design_intent`, `_is_single_cell_command`, `_is_rnaseq_transcriptomics_command`, `_is_tabular_qa_command`, `_is_tabular_analysis_command`, `_has_sequence_cues`, `_should_sequence_alignment_fallback`. |
| 6 | `backend/action_plan.py` | Replace `infer_action_type` keyword sets with LLM classification. |
| 7 | `backend/orchestration/execution_router.py` | `normalize_tool_selection` → delegate to LLM. |

### 🟡 P2 — Domain classification (own PR)

| # | File | What to change |
|---|------|---------------|
| 8 | `backend/workflow_planner_agent.py` | Replace `RNASeqPlaybook.matches`, `ScRNASeqPlaybook.matches`, etc. with an LLM-based playbook selector that returns `PlaybookType`. |
| 9 | `backend/workflow_executor.py` | `_is_multi_step_command` → LLM. |
| 10 | `backend/agent.py` | `_extract_inline_pipeline_plan` — replace `INLINE_PIPELINE_TOOL_MAP` substring loop with LLM step classifier. `_handle_multi_step_workflow` markers → LLM. |
| 11 | `backend/main.py` `_looks_like_workflow` | Replace delimiter + phrase heuristics with LLM. |

### 🟢 P3 — Contextual helpers (own PR)

| # | File | What to change |
|---|------|---------------|
| 12 | `backend/orchestration/visualization_resolver.py` | `what is / explain / describe` check → LLM. |
| 13 | `backend/orchestration/artifact_resolver.py` | Title substring matching → LLM. |
| 14 | `backend/context_builder.py` | Session brief heuristics → LLM. |
| 15 | `backend/routing_keywords.py` | Delete file after P1 complete. |

---

## Unit test strategy

**Principle:** Tests should NEVER pass because a keyword fallback caught something.

### Before (wrong)
```python
# Tests relied on HELIX_MOCK_MODE=1 disabling LLM, then keyword fallback routing
def test_align_command_routes_to_sequence_alignment():
    router = CommandRouter()
    tool, _ = router.route_command("align these sequences", {})
    assert tool == "sequence_alignment"  # ← keyword matched, LLM never called
```

### After (correct)
```python
@patch("backend.command_router.CommandRouter._route_with_llm")
def test_align_command_routes_to_sequence_alignment(mock_route):
    mock_route.return_value = ("sequence_alignment", {"session_id": ""})
    router = CommandRouter()
    tool, _ = router.route_command("align these sequences", {})
    assert tool == "sequence_alignment"
    mock_route.assert_called_once()  # ← LLM was called

def test_route_raises_when_llm_unavailable():
    with patch("backend.command_router.CommandRouter._route_with_llm", return_value=None):
        with pytest.raises(RoutingError):
            router = CommandRouter()
            router.route_command("align these sequences", {})
```

### Test migration checklist
- [ ] `test_command_router_shadow_mode.py` — replace keyword fixtures with LLM mocks
- [ ] `test_routing_flow.py` — all routing assertions must go through mocked LLM
- [ ] `test_intent_classifier.py` — remove all heuristic fallback tests; add test that LLM failure raises
- [ ] `test_approval_classifier.py` — remove keyword fallback test; add test that LLM failure raises
- [ ] `test_orchestration_modules.py` — replace approval_policy keyword tests with LLM mocks
- [ ] `test_tabular_e2e_smoke.py` — router mocks must use LLM mock, not keyword fallback

---

## Environment variable changes

| Variable | Old default | New default | Meaning |
|----------|-------------|-------------|---------|
| `HELIX_LLM_ROUTER_FIRST` | `"0"` | `"1"` | LLM is always the primary router |
| `HELIX_KEYWORD_ROUTING` | (not existed) | `"0"` | Set to `"1"` ONLY for local debug / migration |
| `HELIX_USE_LLM_ROUTER` | `"1"` | `"1"` | No change |

---

## Implementation order (this session)

1. ✅ `approval_policy.py` — wire in LLM classifier, delete legacy `is_approval_command`
2. ✅ `command_router.py` — flip `HELIX_LLM_ROUTER_FIRST` default to `"1"`
3. ✅ `intent_classifier.py` — remove heuristic fallback, raise on LLM failure
4. ✅ `approval_classifier.py` — remove keyword fallback, raise on LLM failure
5. ✅ Update affected unit tests

P1/P2/P3 go in subsequent PRs.
