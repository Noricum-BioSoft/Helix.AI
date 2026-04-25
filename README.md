# Helix.AI

**Helix.AI** is an autonomous bioinformatics agent: it takes a natural-language analysis
request, proposes a plan, waits for approval, then **generates code, chooses infrastructure,
and executes — end to end, without the user touching a terminal.**

---

## How it differs from general AI coding assistants

Tools like Cursor or Claude Code help you write and edit code. Helix goes further:

| | Cursor / Claude Code | Helix.AI |
|---|---|---|
| **Execution** | Generates code; you run it | Generates *and runs* code in a sandbox |
| **Infrastructure** | Your local machine only | Chooses Local / EC2 / EMR per run based on data size and cost |
| **Approval loop** | Applies changes immediately | Proposes a plan; nothing runs until you confirm |
| **Domain knowledge** | General-purpose | Bioinformatics-first: 20+ built-in tools, biological file formats, scientific data awareness |
| **Reproducibility** | You manage the outputs | Every run produces `analysis.py` + `bundle.zip` with full provenance |
| **Iteration** | Regenerate from scratch | Patches the actual script and re-runs; every run is linked to its parent |

---

## How it works — Plan → Approve → Execute

```
1. Plan      Helix reads the request, selects tools, and proposes a
             step-by-step analysis plan in plain language.

2. Approve   You review the plan and confirm before anything runs.
             Nothing executes without explicit approval.

3. Execute   Helix generates code, picks the right infrastructure
             (local sandbox / EC2 / EMR), runs it, and returns
             real artifacts — plots, tables, downloadable bundles.
```

---

## Capabilities

| | |
|---|---|
| **Universal biological data** | Upload CSV, Excel, FASTQ, FASTA, BAM/SAM, VCF, H5AD and more; profiled at upload time |
| **Code Interpreter** | Generates and executes sandboxed Python scripts for tabular and custom analyses |
| **Infrastructure routing** | Per-run decision: local / EC2 / EMR based on data size and cost |
| **20+ built-in tools** | scRNA-seq, bulk RNA-seq, FastQC, alignment, phylogenetics, NCBI/UniProt, GO enrichment, plasmid visualization, variant calling, read trimming, and more |
| **Tool generation** | No matching tool? `ToolGeneratorAgent` writes a new one on-demand |
| **Iterative refinement** | "Change the threshold to 0.01" patches the script and re-runs — not a new LLM prompt |
| **Reproducibility bundle** | Every run: `analysis.py` (standalone, editable) + `bundle.zip` (script + plots + tables + full history) |

---

## Architecture

```
User prompt
     │
     ▼
 Plan          LLM selects tools and writes a plain-language analysis plan.
     │
     ▼
 Approve       Plan shown to user. Nothing runs until confirmed.
     │         (Plan persists across browser sessions via WorkflowCheckpoint.)
     ▼
 Generate      BioAgent writes the analysis code (Python).
     │
     ▼
 Route         Infrastructure decision: Local sandbox / EC2 / EMR
     │         chosen per-run based on data size and cost.
     ▼
 Execute       ExecutionBroker runs the code, validates outputs,
               and returns real artifacts — plots, tables, bundles.
```

The LLM decides *what* to do and *how to code it*; `ExecutionBroker` decides
*where to run it* and validates every output — the LLM never touches I/O directly.

For the full request routing detail see [`docs/architecture/ROUTING_FLOW.md`](docs/architecture/ROUTING_FLOW.md).

---

## Quick start (local dev)

**Prerequisites:** Python 3.10+, Node 18+

```bash
# 1. Clone
git clone https://github.com/your-org/Helix.AI && cd Helix.AI

# 2. Configure
cp backend/.env.example .env
# Edit .env — add at least one of: OPENAI_API_KEY, DEEPSEEK_API_KEY

# 3. Start everything
./start.sh
```

Frontend: `http://localhost:5173`  
Backend API + docs: `http://localhost:8001/docs`

---

## Configuration

Copy `backend/.env.example` to `.env` at the repo root and fill in:

| Variable | Required | Description |
|---|---|---|
| `OPENAI_API_KEY` | One of these two | OpenAI API access |
| `DEEPSEEK_API_KEY` | One of these two | DeepSeek model access |
| `HELIX_AGENT_OPENAI_MODEL` | Optional | Main orchestration agent. Default: `openai:gpt-5.5` |
| `HELIX_INTENT_OPENAI_MODEL` | Optional | Intent classifier. Default: `openai:gpt-5.5` |
| `HELIX_INFRASTRUCTURE_OPENAI_MODEL` | Optional | Infrastructure decision agent. Default: `openai:gpt-5.5` |
| `HELIX_TOOLGEN_OPENAI_MODEL` | Optional | Tool-generation agent. Default: `openai:gpt-5.5` |
| `HELIX_TABULAR_QA_MODEL` | Optional | Tabular Q&A agent. Default: `openai:gpt-5.5` |
| `HELIX_ANALYSIS_PLANNER_MODEL` | Optional | Code Interpreter planner. Default: `openai:gpt-5.5` |
| `HELIX_ANALYSIS_EXECUTOR_MODEL` | Optional | Code Interpreter executor. Default: `openai:gpt-5.5` |
| `HELIX_APPROVAL_OPENAI_MODEL` | Optional | Approval intent classifier. Default: `openai:gpt-4.1-mini` (lightweight yes/no) |
| `HELIX_STAGING_OPENAI_MODEL` | Optional | Staging classifier. Default: `openai:gpt-4.1-mini` (lightweight yes/no) |
| `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` | Optional | S3 data access and EC2/EMR routing |
| `HELIX_USE_EC2` | Optional | Set `true` to enable EC2 job routing |
| `NCBI_API_KEY` | Optional | Higher NCBI rate limits for sequence queries |
| `HELIX_MOCK_MODE` | Dev/test only | Set `1` to disable all LLM calls (unit-test mode) |

---

## Running the demos

The frontend ships with five built-in demo scenarios accessible from the sidebar:

| Demo | Tool | What it shows |
|---|---|---|
| **Bulk RNA-seq** | `bulk_rnaseq_analysis` | Factorial DE analysis, volcano plot, MA plot, iterative parameter changes |
| **Single-Cell RNA-seq** | `single_cell_analysis` | QC → normalization → UMAP → marker genes on SLE PBMC dataset |
| **FastQC / Sequencing QC** | `fastqc_analysis` | Per-base quality, adapter contamination, async HTML report |
| **Phylogenetic Tree** | `phylogenetic_tree` | Interactive tree rendering from aligned sequences |
| **Read Trimming** | `read_trimming` | Adapter trimming, quality filtering, paired-end support |

Click **Load & Run** on any scenario to pre-fill the input and execute with one click.

---

## Iterative workflow example

```
User:   Upload counts.csv and metadata.csv — run bulk RNA-seq analysis.
Helix:  Here is my plan:
          Step 1 — Load counts matrix (counts.csv) and metadata (metadata.csv)
          Step 2 — Normalize and run DESeq2-style differential expression
          Step 3 — Generate volcano plot and DE table
        Approve to proceed, or ask me to revise any step.

User:   Proceed.
Helix:  → Executes plan, returns volcano plot + DE table + analysis.py

User:   Change the significance threshold to 0.01.
Helix:  → Patches ALPHA = 0.01 in analysis.py, re-runs, returns updated plots
          with a diff showing exactly what changed between runs.

User:   Switch x-axis to linear fold change.
Helix:  → Patches X_SCALE = "linear", re-runs, shows new plot alongside old.
```

Each step is a new run linked to the previous one. The `bundle.zip` contains the full chain.

---

## Reproducibility

Every completed run produces two downloadable artifacts:

- **`analysis.py`** — a self-contained Python script with a clearly marked `# ── Parameters ──`
  block. Edit the parameters and run it directly: `python analysis.py`
- **`bundle.zip`** — contains `analysis.py`, all plots and tables, `run_manifest.json`,
  `README.md`, and an `iteration_history/` directory with the full ancestor chain

---

## Repo structure

```
backend/          FastAPI backend — agents, tools, orchestration, LLM classifiers
  orchestration/  LLM-based intent, approval, and staging classifiers
  tabular_qa/     Code Interpreter for tabular/CSV/Excel analysis
frontend/         React + Vite UI (TypeScript)
tools/            Bioinformatics tool implementations
tests/
  unit/           pytest unit suite (~700 tests); runs with HELIX_MOCK_MODE=1
  integration/    Integration tests against a live backend
  benchmarks/     P0 benchmark gate cases and scorer
scripts/          Utility and deployment scripts
infrastructure/   AWS CDK stack (EC2, EMR, S3)
docs/             Architecture and deployment documentation
benchmarks/       Benchmark definitions and release thresholds
```

---

## Tests

```bash
# All backend unit tests (no live LLM required)
HELIX_MOCK_MODE=1 pytest tests/unit/backend/ -q

# Specific test groups
pytest tests/unit/backend/test_approval_classifier.py   # LLM approval classifier
pytest tests/unit/backend/test_intent_classifier.py     # LLM intent classifier
pytest tests/unit/benchmarks/                           # P0 benchmark gates

# Frontend type check + build
cd frontend && npm run build
```

Tests run with `HELIX_MOCK_MODE=1` by default (see `pytest.ini`). LLM-dependent
functions are patched via autouse fixtures in `tests/unit/backend/conftest.py` — no
live API keys needed for the unit suite.

---

## Deployment

See [`docs/deployment/AWS_DEPLOYMENT_GUIDE.md`](docs/deployment/AWS_DEPLOYMENT_GUIDE.md) for
EC2 and EMR deployment with the CDK stack.

For the cheapest beta option (single EC2 instance), see
[`docs/deployment/BETA_DEPLOYMENT_CHEAP.md`](docs/deployment/BETA_DEPLOYMENT_CHEAP.md).
