# Helix.AI Bioinformatics Agent Plan

## Current System Snapshot
- Frontend (`frontend/`): React-based UI with prompt input, file uploads, collapsible examples sidebar, chat-style history, and workflow context panels. Sends commands via `mcpApi` to backend REST endpoints.
- Backend (`backend/main_with_mcp.py`): FastAPI service managing sessions, routing commands, and exposing MCP-backed bioinformatics operations (sequence alignment, mutation, variant selection, plasmid visualization, etc.).
- Tools (`tools/`): Python modules implementing specific analyses such as `alignment`, `mutations`, `variant_selection`, `phylogenetic_tree`, `plasmid_visualizer`, and general analytics in `bio` and `data_science`.
- MCP Integration: `call_mcp_tool` maps logical tool names to Python implementations. History manager tracks session context for downstream steps.
- Gap Analysis: No current tooling for read trimming, read merging/assembly, or dataset introspection prior to planning.

## Proposed End-State Architecture
```
┌──────────────────────────────────────────────────────────────────┐
│                           Frontend (React)                       │
│   - Prompt input, file upload, history                           │
│   - Sends user intent/files → /agent endpoint                    │
│   - Renders agent answers + tool outputs                         │
└──────────────────────────────────────────────────────────────────┘
                │ REST/JSON (prompt, files, session)
                ▼
┌──────────────────────────────┐      ┌────────────────────────────┐
│  FastAPI Backend (Helix.AI)  │◀────▶│  Session / History Store   │
│  - `/agent` endpoint         │      │  - Tracks workflow context │
│  - File handling, auth       │      │  - Supplies past outputs   │
└──────────────┬───────────────┘      └─────────-─────┬────────────┘
               │ delegates to agent orchestrator      │context
               ▼                                      │
┌──────────────────────────────────────────────────────────────────┐
│    Bioinformatics Agent Orchestrator (LLM policy / planner)      │
│  - Interprets prompts, asks clarifying questions                 │
│  - Provides Q&A, workflow design, and tool execution directives  │
│  - Emits narrative responses + structured tool calls             │
└──────────────┬──────────────────────────────────────────────-────┘
               │ MCP tool invocations
               ▼
┌──────────────────────────────────────────────────────────────────┐
│                  MCP Server & Tool Adapters                      │
│  - Existing tools (alignment, mutate_sequence, etc.)             │
│  - Planned additions (read trimming, read merging)               │
│  - Wrap external utilities (Clustal, MUSCLE, QC tools)           │
└──────────────┬──────────────────────────────────────────────-────┘
               │
               ▼
┌──────────────────────────────────────────────────────────────────┐
│           External Bioinformatics Utilities / Libraries          │
└──────────────────────────────────────────────────────────────────┘

## Agent System Prompt (Draft)
```
You are the Helix.AI Bioinformatics Agent, an expert assistant embedded in a workflow platform. Follow these principles on every request:

1. Core Responsibilities
   - Educate: answer conceptual and practical bioinformatics questions with clear, concise explanations and cite relevant tools or algorithms (e.g., Clustal Omega, MUSCLE, MAFFT, FastQC).
   - Strategize: design end-to-end analysis plans. Gather essential context about the dataset (type, format, size, quality), the user’s research objective, and desired outputs before proposing steps.
   - Execute: when asked (or after explicit confirmation), run single pipeline steps by invoking available MCP tools. Share progress, surface outputs, and record context for future steps.

2. Conversation & Clarification
   - Always confirm you understand the biological question and the input data. Ask targeted follow-ups when information is missing (e.g., file format, read length, organism, sequencing depth).
   - Offer tool recommendations mapped to the user’s goal (e.g., “For multiple sequence alignment, consider MUSCLE or Clustal Omega.”).

3. Tool Usage Guidelines
   - Available MCP tools include:
     • `sequence_alignment` — perform MSA (supports Clustal, MUSCLE, MAFFT).
     • `mutate_sequence` — generate sequence variants.
     • `select_variants` / `variant_selection` — prioritize variants from previous runs.
     • `analyze_sequence_data` — run exploratory analysis/visualizations.
     • `visualize_alignment` — render alignment plots.
     • `plasmid_visualization` / `plasmid_for_representatives` — build plasmid maps.
     • `phylogenetic_tree` / `clustering_analysis` — infer phylogeny and clusters.
   - Use new trimming/merging tools once they are available (ask backend for `read_trimming`, `read_merging`, etc.).
   - Always verify that required inputs for the chosen tool are supplied (e.g., FASTA sequences for MSA, aligned data for phylogenetic trees, plasmid context for visualization). Request missing data or halt execution until prerequisites are met.
   - Before executing, confirm parameters with the user unless they provided explicit instructions.
   - After each tool call, summarize results, note any files generated, and update the conversation’s context.

4. Response Format
   - Begin with a concise answer or plan.
   - When executing tools, include a “Tool Actions” section describing calls and inputs.
   - Surface assumptions, required next steps, and recommendations.
   - Encourage users to proceed step-by-step and highlight potential quality-control checks.

5. Safety & Boundaries
   - Do not fabricate tool results. If a tool fails or is unavailable, explain the issue and suggest alternatives.
   - Respect data privacy: never share data between sessions or make assumptions about user data ownership.
```

## Development Roadmap
1. **Backend Enhancements**
   - Add `/agent` FastAPI endpoint that forwards prompts, session context, and uploaded files to the agent orchestrator.
   - Implement read trimming and read merging MCP tools (wrappers around e.g., `cutadapt`, `trimmomatic`, `FLASH` or equivalent Python libraries).
   - Extend MCP tool metadata endpoints (`/mcp/tools`) to reflect new capabilities.

2. **Agent Orchestrator**
   - Integrate an LLM runner (OpenAI, Anthropic, etc.) with the drafted system prompt.
   - Implement structured response parsing (text + tool call directives).
   - Maintain conversation state, including previous tool outputs and clarifying questions.

3. **Frontend Updates**
   - Route all prompt submissions through the new `/agent` endpoint.
  - Display agent reasoning, clarifications, and tool plans clearly (e.g., collapsible “Agent Plan” cards).
   - Surface tool execution progress and outputs in the chat history.

4. **Testing & Validation**
   - Create scenario tests for education, planning, and execution modes.
   - Validate trimming/merging pipelines with sample sequencing datasets.
   - Ensure session context flows end-to-end between frontend, backend, agent, and tools.

5. **Documentation & Onboarding**
   - Update developer docs with agent architecture, prompt contract, and tool catalogue.
   - Provide user-facing guidance on interacting with the agent (commands, required metadata, confirmation workflow).

Once these steps are complete, Helix.AI will deliver a cohesive agent-driven experience that educates users, plans analyses, and executes bioinformatics workflows end-to-end.

