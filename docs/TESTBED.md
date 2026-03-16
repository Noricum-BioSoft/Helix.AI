# Helix.AI Testbed — Scenario-Based Test Definition

This document defines the **testbed** for Helix.AI: every test derived from [USER_SCENARIOS.md](USER_SCENARIOS.md) to ensure the system works correctly on both **backend** and **frontend**. Each scenario ID (A1–E4, F1–F6, G1–G3, H1–H2, I1–I2) is mapped to concrete test cases with layer, endpoint/component, assertions, and dependencies.

**Test layers**

| Layer | Scope | When run |
|-------|--------|----------|
| **Backend unit** | Single module or function; TestClient or mocks; no live LLM/S3 | Every commit; fast |
| **Backend integration** | Full request/response via TestClient(app); may use mock mode | Every commit; medium |
| **Backend e2e** | Live backend (optional live LLM/S3); real tools | CI or nightly; slow |
| **Frontend unit** | React component or hook in isolation | Every commit; fast |
| **Frontend integration** | UI + mock API or test backend | PR / nightly |
| **E2E (browser)** | Real browser vs running app (Playwright/Cypress) | Nightly / release |

---

## 1. Cross-cutting tests (session, state, correctness)

These apply across multiple scenarios.

| Test ID | Layer | Description | Endpoint / code | Key assertions |
|---------|--------|-------------|------------------|----------------|
| **CC1** | Backend integration | New session creates session_id; no runs yet | POST /session/create or POST /create_session | 200; body has session_id |
| **CC2** | Backend integration | Execute with session_id attaches run to that session | POST /execute, GET /session/{sid}/runs | run in runs list; run_id present |
| **CC3** | Backend integration | needs_inputs response does not create a run | POST /execute with prompt that yields needs_inputs; GET /session/{sid}/runs | status=needs_inputs; run count unchanged (or only router/agent run, no tool run) |
| **CC4** | Backend unit | build_standard_response sets error when result.status=error | Call build_standard_response with success=False or status=error | response["error"] is set; response["success"] is False |
| **CC5** | Frontend | New Session button clears history and gets new session_id | UI click "New Session" | history empty; session_id updated; no stale run links |

---

## 2. A. Exploration & tool discovery

| Scenario | Test ID | Layer | Trigger (command or UI) | Key assertions |
|----------|---------|--------|--------------------------|----------------|
| **A1** | A1_be | Backend integration | POST /execute "I have bulk RNA-seq count data and sample metadata. What can you do with it?" | status or result.status = needs_inputs; text contains "Required inputs", "count_matrix", "sample_metadata", table pipes \| |
| **A1** | A1_fe | Frontend | Same prompt via UI; inspect response block | Rendered markdown shows tables (GFM); "Provide the required inputs above" visible |
| **A2** | A2_be | Backend integration | POST /execute "What do I need to run single-cell RNA-seq analysis?" | needs_inputs; text contains data_file, data_format or resolution |
| **A2** | A2_fe | Frontend | Same via UI | Tables and required/optional parameters visible |
| **A3** | A3_be | Backend integration | POST /execute "How do I run FastQC on my FASTQ files?" | needs_inputs; text contains input_r1, input_r2 |
| **A3** | A3_fe | Frontend | Same via UI | S3 path examples visible |
| **A4** | A4_be | Backend integration | POST /execute "What's the difference between bulk and single-cell RNA-seq?" | success; no tool execution (intent=ask); text is markdown; no run created for a tool like bulk_rnaseq |
| **A4** | A4_fe | Frontend | Same via UI | Answer shown as markdown; no plot/table from execution |

**Correctness (all A)**  
- A1_be, A2_be, A3_be: response text must come from _TOOL_INPUT_REQUIREMENTS (no hallucinated params).  
- A4_be: must not call bulk_rnaseq_analysis or single_cell_analysis.

---

## 3. B. First run (supply data)

| Scenario | Test ID | Layer | Trigger | Key assertions |
|----------|---------|--------|---------|----------------|
| **B1** | B1_be | Backend integration | POST /execute with bulk RNA-seq demo prompt + follow-up with count_matrix, sample_metadata, design_formula (or use followUpPrompt from demo) | success; run_id in response; data.results or data.visuals; data.links or top-level links contain "analysis.py" and "bundle.zip"; session has run |
| **B1** | B1_fe | Frontend | Load Demo 1, click "Load & Run" | Response shows plots/tables; download links for analysis.py and bundle.zip visible |
| **B2** | B2_be | Backend integration | POST /execute with single-cell demo prompt + follow-up (data_file S3 H5, data_format 10x, resolution) | success; run_id; visuals or results; script saved under session/run; links present |
| **B2** | B2_fe | Frontend | Load Demo 2, "Load & Run" | UMAP/markers or placeholder; download links visible |
| **B3** | B3_be | Backend integration | POST /execute "Run FastQC on s3://bucket/R1.fastq.gz and s3://bucket/R2.fastq.gz" (or mock paths) | success or job_id; result references FastQC/report; no fabricated metrics |
| **B3** | B3_fe | Frontend | Same command or demo 3 | Report or link shown; no raw "No plots" only when no artifacts |
| **B4** | B4_be | Backend integration | POST /execute with FASTA sequences + "Build a phylogenetic tree from these sequences." | success; tree_newick or ete_visualization in response; tree from real tool |
| **B4** | B4_fe | Frontend | Paste FASTA, send tree request | Tree visualization rendered; Newick or iframe present |
| **B5** | B5_be | Backend integration | POST /execute "Trim my reads: R1 and R2 on S3, adapter AGATCGGAAGAGC, quality 20" (with paths or inline reads) | success; trimmed output or summary stats in response |
| **B5** | B5_fe | Frontend | Same | Trimming result or summary visible |

**Correctness (all B)**  
- All artifacts (plots, tables, Newick, HTML) must be produced by the actual tool (ExecutionBroker), not LLM prose.  
- If S3 unreachable for B1/B2: synthetic fallback + notice in text (e.g. "Synthetic demo data used").

---

## 4. C. Iterative refinement

| Scenario | Test ID | Layer | Trigger | Key assertions |
|----------|---------|--------|---------|----------------|
| **C1** | C1_be | Backend integration | After B1: POST /execute "Change the significance threshold to 0.01." | success; new run_id; tool or result indicates patch_and_rerun; parent_run_id on new run; parameter_diff or output_diff present; alpha=0.01 in script or result |
| **C1** | C1_fe | Frontend | After B1: type "Change the significance threshold to 0.01", submit | New plot/table shown; run updated; diff or message about change |
| **C2** | C2_be | Backend integration | After B1: POST /execute "Switch the volcano plot x-axis to linear fold change." | success; patch_and_rerun; new run; X_SCALE or x_scale in script/params |
| **C2** | C2_fe | Frontend | Same | Volcano or results view updated |
| **C3** | C3_be | Backend integration | After B1: POST /execute "Use design formula ~infection_status + time_point only, no interaction." | success; design/formula in result matches; contrasts consistent |
| **C4** | C4_be | Backend integration | After B2: POST /execute "Increase clustering resolution to 0.8." | Either patch_and_rerun with resolution=0.8 or clear message to re-run with new param; no silent ignore |

**Correctness (all C)**  
- Patches must apply to the actual analysis.py; re-execution uses patched file.  
- GET /session/{sid}/runs returns parent_run_id chain.

---

## 5. D. Validation & interpretation

| Scenario | Test ID | Layer | Trigger | Key assertions |
|----------|---------|--------|---------|----------------|
| **D1** | D1_be | Backend integration | After B1: POST /execute "Which genes are most significant for the infection effect?" | success; answer in text; no new tool run for bulk_rnaseq; answer consistent with DE (e.g. gene names from result) |
| **D2** | D2_be | Backend integration | POST /execute "What does padj < 0.05 mean?" | success; markdown explanation; no execution of analysis tool |
| **D3** | D3_be | Backend integration | After B1: POST /execute "Show me the same volcano but with log2 scale on the x-axis again." | success; patch_and_rerun or equivalent; x_scale log2 |

**Correctness (all D)**  
- Q&A must not mutate run state or overwrite artifacts.  
- D3 must use last scriptable run for patch.

---

## 6. E. Reproducibility & handoff

| Scenario | Test ID | Layer | Trigger | Key assertions |
|----------|---------|--------|---------|----------------|
| **E1** | E1_be | Backend integration | After B1: GET /download/script?path=... (path from run artifact or response) | 200; Content-Type text/plain or application/x-python; body contains "# ── Parameters ──" and run_id or RUN_DIR |
| **E1** | E1_fe | Frontend | After B1: click "analysis.py" link | Download starts; file is Python script |
| **E2** | E2_be | Backend integration | GET /download/bundle?session_id=...&run_id=... (after B1→C1→C2) | 200; Content-Type application/zip; ZIP contains analysis.py, README.md, run_manifest.json, iteration_history or equivalent; manifest lists inputs/outputs |
| **E2** | E2_fe | Frontend | Click "bundle.zip" | Download starts; ZIP structure as above |
| **E3** | E3_be | Backend unit / integration | (Documentation test) Downloaded analysis.py: parameter block has COUNT_MATRIX, SAMPLE_METADATA, RUN_DIR; script runs with python when data available | Run script in test with mock data; exits 0 or produces expected outputs |
| **E4** | E4_be | Backend integration | Unzip bundle; README describes how to run; manifest has correct run_id and artifact paths | README not placeholder; manifest.run_id matches requested run |

---

## 7. F. Sequence, alignment & plasmid

| Scenario | Test ID | Layer | Trigger | Key assertions |
|----------|---------|--------|---------|----------------|
| **F1** | F1_be | Backend integration | POST /execute "Align these FASTA sequences" + FASTA body | success; alignment in result; usable for tree |
| **F2** | F2_be | Backend integration | POST /execute "Build a tree from these sequences" + unaligned FASTA | success; tree_newick or ete_visualization |
| **F3** | F3_be | Backend integration | POST /execute "Generate 50 random variants of this sequence." + sequence | success; variants or mutated_sequences in result; session has mutation_results |
| **F4** | F4_be | Backend integration | After F3: POST /execute "Select 10 variants with best conservation." (session_id same) | success; 10 sequences in output; from variant_selection |
| **F5** | F5_be | Backend integration | POST /execute "Show my gene in a plasmid map — vector pUC19, EcoRI/BamHI sites." | success; plasmid image or SVG in result |
| **F6** | F6_be | Backend integration | POST /execute "Put these representative sequences into pUC19 and show the plasmid maps." + sequences | success; multiple plasmid visuals |

---

## 8. G. Literature & knowledge

| Scenario | Test ID | Layer | Trigger | Key assertions |
|----------|---------|--------|---------|----------------|
| **G1** | G1_be | Backend integration | POST /execute "Fetch the BRCA1 human sequence from NCBI." | success; FASTA or accession in result; from NCBI API or cached |
| **G2** | G2_be | Backend integration | POST /execute "Look up protein P53 in UniProt." | success; structured summary from UniProt |
| **G3** | G3_be | Backend integration | POST /execute "What is GO:0006915?" | success; term name/definition from go_tools |

---

## 9. H. DNA synthesis

| Scenario | Test ID | Layer | Trigger | Key assertions |
|----------|---------|--------|---------|----------------|
| **H1** | H1_be | Backend integration | POST /execute "Research vendors for synthesizing a 500 bp sequence." | success; vendor options or structured response |
| **H2** | H2_be | Backend integration | POST /execute "Submit this sequence for synthesis to Twist." + sequence | success; submission confirmation or mock indicator; UI must not imply real order if mock |

---

## 10. I. Multi-step pipelines

| Scenario | Test ID | Layer | Trigger | Key assertions |
|----------|---------|--------|---------|----------------|
| **I1** | I1_be | Backend integration | POST /execute "Align these sequences, then build a phylogenetic tree." + FASTA | success; tree from same sequences; alignment used |
| **I2** | I2_be | Backend integration | POST /execute "Generate 96 variants of this sequence, select the 12 with lowest gaps, and show them in plasmid maps." + sequence | success; 12 plasmid maps or equivalent; all three steps executed; session state passed |

---

## 11. Frontend-focused test cases

These verify UI behavior; they can be manual checklists or automated (e.g. Playwright) when available.

| Test ID | Description | Key checks |
|---------|-------------|------------|
| **FE_session** | New Session button | Clears history; new session_id; no stale run links |
| **FE_demo_load** | Load demo scenario | Prompt and optional "Load & Run" / "Use example data" appear; scenario title visible |
| **FE_needs_inputs** | needs_inputs response | Markdown tables render (no raw \| \|); "Provide the required inputs above" visible |
| **FE_results_links** | After any run with artifacts | analysis.py and bundle.zip links visible for that response (all demos, not only Demo 1) |
| **FE_prompt_format** | User prompt in conversation | Whitespace and newlines preserved (pre-wrap) |
| **FE_error_display** | On API error | Alert or inline error message shown; no blank response |
| **FE_jobs_panel** | Jobs panel | When job_id in response, job appears in panel; status updates if polled |

---

## 12. Planned / roadmap scenarios (J–N)

Scenarios J1–J3, K2–K4, L2–L3, M1–M3 are **planned**. N1–N3 are **partially supported** (provenance exists). Tests can be added as stubs (e.g. skip with reason "VCF tool not implemented") or omitted until tools exist.

| Scenario group | Test approach |
|----------------|----------------|
| **J** (VCF) | Skip or placeholder test: "J1_be: annotate VCF → skip (no VCF tool)". |
| **K2–K4** | K1_be: already covered by B3/B5 + read_trimming/read_merging. K2–K4: skip until taxonomy/DA/diversity tools exist. |
| **L2–L3** | L1_be: G3 covers lookup_go_term. L2/L3: skip or partial (e.g. single_cell pathways step). |
| **M** | Skip until biomarker tools exist. |
| **N** | N1_be: GET run + manifest; assert run.json and manifest structure. N2_be: bundle contains README. N3_be: manifest has run_id and artifact refs. |

---

## 13. Test implementation layout

- **Backend tests (pytest)**  
  - `tests/testbed/conftest.py` — fixtures: TestClient(app), temporary session dir, demo payloads.  
  - `tests/testbed/test_scenario_a_exploration.py` — A1_be–A4_be.  
  - `tests/testbed/test_scenario_b_first_run.py` — B1_be–B5_be.  
  - `tests/testbed/test_scenario_c_iteration.py` — C1_be–C4_be.  
  - `tests/testbed/test_scenario_d_validation.py` — D1_be–D3_be.  
  - `tests/testbed/test_scenario_e_reproducibility.py` — E1_be–E4_be.  
  - `tests/testbed/test_scenario_f_sequence_plasmid.py` — F1_be–F6_be.  
  - `tests/testbed/test_scenario_g_literature.py` — G1_be–G3_be.  
  - `tests/testbed/test_scenario_h_dna_synthesis.py` — H1_be–H2_be.  
  - `tests/testbed/test_scenario_i_multistep.py` — I1_be–I2_be.  
  - `tests/testbed/test_cross_cutting.py` — CC1–CC5 (backend + frontend checklist).

- **Frontend tests**  
  - No Jest/Vitest in repo yet: document in this TESTBED as manual or future Playwright.  
  - When added: mirror FE_* cases above (session, demo load, needs_inputs rendering, results links, prompt format, error display, jobs panel).

- **Running the testbed**  
  - Backend (mock mode, no live server):  
    `pytest tests/testbed/ -v`  
  - Backend with live server (integration):  
    `pytest tests/testbed/ -v --backend-url http://localhost:8001` (if supported by fixtures).  
  - Frontend (when implemented):  
    `npm run test` (unit) or `npx playwright test` (e2e).

---

## 14. Traceability

| Source | Target |
|--------|--------|
| USER_SCENARIOS.md §3 (A–E) | TESTBED.md §2–6, tests/testbed/test_scenario_*.py |
| USER_SCENARIOS.md §6 (F–I) | TESTBED.md §7–10, tests/testbed/test_scenario_*.py |
| USER_SCENARIOS.md §7–8 (J–N) | TESTBED.md §12 (stubs/skips) |
| USER_SCENARIOS.md §4 (session/state) | TESTBED.md §1 (CC1–CC5), test_scenario_* session usage |
| USER_SCENARIOS “Correctness” bullets | Assertions in each test row (real tool output, no hallucination, parent_run_id, etc.) |
