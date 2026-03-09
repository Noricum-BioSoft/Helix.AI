# Helix.AI

**Helix.AI** is a multi-agent web application that turns natural-language biology questions into
orchestrated, auditable, and reproducible bioinformatics analyses.

A scientist types a question in plain English. Helix routes it through intent detection, tool
selection, infrastructure decisions, and execution — without the user touching a command line,
writing a config file, or running a script.

---

## What it does

| Capability | Description |
|---|---|
| **Natural language → workflow** | Describe your analysis goal; Helix selects the right tool, collects missing inputs, and executes |
| **Iterative science** | Every run is a tracked child of the previous one (`parent_run_id`). Ask "change the x-axis to log2" and Helix patches the *actual script* and re-runs — not re-prompts the LLM |
| **Reproducibility bundle** | Every run produces `analysis.py` (standalone, editable) and `bundle.zip` (script + all plots, tables, manifest, full history) |
| **16+ built-in tools** | scRNA-seq, bulk RNA-seq, FastQC/QC, alignment, phylogenetics, NCBI/UniProt, GO enrichment, plasmid visualization, variant selection, read trimming, DNA synthesis, and more |
| **Tool generation** | When a tool is missing, the `ToolGeneratorAgent` creates a new wrapper on-demand rather than failing |
| **Cloud routing** | Infrastructure selection (Local / EC2 / EMR) is a first-class per-run decision based on data size and cost |

---

## Architecture

```
User (natural language prompt)
        │
   IntentDetector
        │
   ┌────┴────────────────────────────┐
   │ ask path         execute path   │
   │                                 │
BioinformaticsGuru           WorkflowPlannerAgent
                                     │
                              ImplementationAgent
                                     │
                        InfrastructureDecisionAgent
                                     │
                            ┌────────┴───────────────┐
                            │    Tool exists?        │
                            │   No → ToolGenerator   │
                            │  Yes → ExecutionBroker |
                            |        (jobs, I/O)     |
                            └────────┬───────────────┘
                               DataVisualizer
                                     │
                                 Response 
                     (report + artifacts + download links)
```

**Key structural properties:**
- The LLM *proposes* actions; typed Pydantic contracts and a non-LLM `ExecutionBroker` *validate and execute* them — the LLM never directly touches I/O
- `HandoffPolicy` validates every agent transition; illegal handoffs raise `PolicyViolationError`
- Tools return real artifacts (plots, tables, Newick trees) against validated schemas — results cannot be hallucinated
- Every run writes `sessions/<sid>/runs/<run_id>/run.json` and `analysis.py`; runs are linked by `parent_run_id` for full iteration history

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
| `OPENAI_API_KEY` | One of these two | OpenAI API access (model configurable, see below) |
| `DEEPSEEK_API_KEY` | One of these two | DeepSeek model access |
| `HELIX_INTENT_OPENAI_MODEL` | Optional | Intent classifier model (e.g. `openai:gpt-5.2`, `openai:gpt-4o`). Default: `openai:gpt-5.2`. |
| `HELIX_INFRASTRUCTURE_OPENAI_MODEL` | Optional | Infrastructure decision model. Default: `openai:gpt-5.2`. |
| `HELIX_TOOLGEN_OPENAI_MODEL` | Optional | Tool-generation model. Default: `openai:gpt-5.2`. |
| `AWS_ACCESS_KEY_ID` / `AWS_SECRET_ACCESS_KEY` | Optional | S3 data access and EC2/EMR routing |
| `HELIX_USE_EC2` | Optional | Set `true` to enable EC2 job routing |
| `NCBI_API_KEY` | Optional | Higher NCBI rate limits for sequence queries |

---

## Running the demos

The frontend ships with five built-in demo scenarios accessible from the sidebar:

| Demo | Tool | What it shows |
|---|---|---|
| **Bulk RNA-seq** | `bulk_rnaseq_analysis` | Factorial DE analysis (DESeq2 equivalent), volcano plot, MA plot, iterative parameter changes |
| **Single-Cell RNA-seq** | `single_cell_analysis` | QC → normalization → UMAP → marker genes on SLE PBMC dataset |
| **FastQC / Sequencing QC** | `fastqc_analysis` | Per-base quality, adapter contamination, async HTML report |
| **Phylogenetic Tree** | `phylogenetic_tree` | Interactive tree rendering from aligned sequences |
| **Read Trimming** | `read_trimming` | Adapter trimming, quality filtering, paired-end support |

Click **Load & Run** on any scenario to pre-fill the input and execute with one click.

---

## Iterative workflow example

```
User:  Run bulk RNA-seq on s3://my-bucket/counts.csv with metadata s3://my-bucket/meta.csv
       → Helix executes, returns volcano plot + DE table + analysis.py download

User:  Change the significance threshold to 0.01
       → Helix patches ALPHA = 0.01 in analysis.py, re-runs, returns updated plots
         with a diff showing what changed between runs

User:  Switch x-axis to linear fold change
       → Helix patches X_SCALE = "linear", re-runs, shows new plot alongside old
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
backend/          FastAPI backend — agents, tools, orchestration, MCP server
frontend/         React + Vite UI
tools/            Bioinformatics tool implementations
scripts/          Utility and deployment scripts
infrastructure/   AWS CDK stack (EC2, EMR, S3)
docs/             Architecture and deployment documentation
tests/            pytest suite + demo scenario fixtures
```

---

## Tests

```bash
# Backend
pytest

# Frontend build check
cd frontend && npm run build
```

---

## Deployment

See [`docs/deployment/AWS_DEPLOYMENT_GUIDE.md`](docs/deployment/AWS_DEPLOYMENT_GUIDE.md) for
EC2 and EMR deployment with the CDK stack.
