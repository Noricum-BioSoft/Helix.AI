# 🧬 Helix.AI - Bioinformatics AI Platform

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Node.js](https://img.shields.io/badge/Node.js-16+-green.svg)](https://nodejs.org/)
[![React](https://img.shields.io/badge/React-18+-blue.svg)](https://reactjs.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.100+-green.svg)](https://fastapi.tiangolo.com/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An AI-powered web application for managing biotechnology workflows via natural language commands, featuring interactive phylogenetic tree visualization, simulated DNA synthesis vendor research, and comprehensive bioinformatics tools with session management and history tracking.

## 📌 Current Status (as of Version 2.0)

- **Architecture**: ✅ **Refactored and improved** — centralized `agent_tools.py` module, simplified agent codebase, better separation of concerns.
- **Execution architecture**: ✅ **Phase 1 complete** — centralized `ExecutionBroker`, routing policy v1, input discovery/size estimation, generalized jobs.
- **Sync execution**: ✅ Supported (local / EC2 depending on deployment & environment).
- **Async execution**: ✅ Supported for **FastQC on EMR** (returns `job_id` immediately).
- **Jobs API**: ✅ `/jobs/{job_id}` + `/jobs/{job_id}/results` working (and additional job helpers exist).
- **Testing**: ✅ Comprehensive test suite with unit tests (mock mode by default), integration tests (opt-in), evaluation framework, and deployed backend tests.
- **Agent Tools**: ✅ 16+ bioinformatics tools centrally defined and organized in `agent_tools.py`.
- **UI/UX**: ✅ Enhanced with ThinkingIndicator component, improved component organization, and better user feedback.

## 🚀 Features

### 🧭 Tool‑first Dispatch (not just an LLM wrapper)
- **Hybrid routing (auto + override)**: Prompts are routed to either:
  - **Direct tool execution (Submit)** when required inputs are present and the intent matches a known operation.
  - **Agent planning (Agent)** for open-ended questions, missing inputs, or multi‑step workflows.
- **Validation before execution**: Required parameters (e.g., sequences for alignment) are checked up front.
- **Provenance & reproducibility**: Each run records tool name, version, parameters, inputs, and outputs.
- **Explainability**: When using the agent, the final answer is rendered in Markdown; intermediate planner chatter is hidden.

### 🧱 Execution Broker + Routing Policy (Phase 1)
- **Single execution entrypoint**: `/execute` routes tool runs through a central **ExecutionBroker**.
- **Consistent execution envelope**: Results include a stable wrapper with `routing`, discovered `inputs`, `artifacts`, and the original tool output preserved under `result`.
- **Routing policy v1 (sync vs async)**:
  - **Bytes threshold** (default **100MB**) with env override: `HELIX_ASYNC_BYTES_THRESHOLD`
  - **Tool overrides** (curated list; currently forces FastQC async)
  - **Timeout promotion hook**: present as a stub for later phases
- **Input discovery + size estimation**:
  - Detects **S3 URIs** (`s3://bucket/key`) and **local paths** in tool arguments
  - Incorporates **session uploads** and **dataset references** from session metadata (when available)
  - Best-effort sizing via S3 `head_object(ContentLength)` and local `stat()`

### 🌳 **Interactive Phylogenetic Tree Visualization**
- **ETE3 Integration**: High-quality phylogenetic tree visualization with SVG rendering
- **D3.js Tree Rendering**: Interactive phylogenetic trees with zoom and pan
- **Newick Format Support**: Parse and display standard tree formats
- **Tree Statistics**: Display comprehensive tree metrics and relationships
- **Real-time Updates**: Trees update as you analyze new sequences
- **Clustering Analysis**: Select representative sequences from phylogenetic clusters

### 🔬 **DNA Synthesis Vendor Research (Simulated)**
- **Vendor Comparison**: Compare 7 major DNA synthesis vendors with simulated data
- **Pricing Analysis**: Simulated pricing ranges and service comparisons
- **Service Matching**: Find vendors based on sequence length and quantity
- **Recommendations**: AI-powered vendor selection advice
- **Note**: Currently uses simulated data. Real integration requires coordination with vendors for API access and pricing agreements.

### 🤖 **Natural Language Commands**
- **Intelligent Command Parsing**: Understand commands like "visualize the phylogenetic tree"
- **Multi-step Workflows**: Chain commands together for complex bioinformatics pipelines
- **Context Awareness**: Maintains session history for continuous workflows
- **Command Routing**: Smart routing to appropriate bioinformatics tools

### 📊 **Interactive Visualizations**
- **Plotly Integration**: Dynamic charts for sequence analysis, mutations, and alignments
- **SeqViz Integration**: Professional plasmid and vector visualization
- **Real-time Updates**: Live visualizations that update as you work
- **Export Capabilities**: Save plots and results in multiple formats

### 🧰 **Comprehensive Tool Suite**
All tools are centrally defined in `backend/agent_tools.py` for better maintainability and organization.

**Core Sequence Analysis:**
- **Sequence Alignment**: Multiple algorithms (ClustalW, Muscle, MAFFT)
- **Mutation Analysis**: Generate and analyze sequence variants
- **Phylogenetic Tree Construction**: Build and visualize evolutionary relationships
- **Sequence Selection**: Smart selection based on diversity, length, conservation, or custom criteria
- **Read Merging**: Merge paired-end reads with quality control

**Visualization:**
- **Plasmid Visualization**: Interactive plasmid and vector visualization with circular/linear views
- **Phylogenetic Tree Rendering**: Interactive tree visualization with D3.js and ETE3
- **Clustering Analysis**: Hierarchical clustering with representative sequence selection

**Genomics & Proteomics:**
- **NCBI Sequence Fetching**: Retrieve sequences from NCBI databases by accession
- **UniProt Query**: Query UniProt protein database for sequences and metadata
- **Gene Ontology (GO) Lookup**: Look up GO term details and annotations

**RNA-seq Analysis:**
- **Single-Cell RNA-seq**: Comprehensive scRNA-seq analysis (Seurat-based pipeline)
  - Preprocessing and quality control
  - Marker gene identification
  - Differential expression analysis
  - Pathway enrichment
  - Cell-type annotation
  - Batch correction
- **Bulk RNA-seq**: Differential expression analysis using DESeq2
- **FastQC Quality Control**: Quality assessment for FASTQ files (async EMR jobs)

**DNA Synthesis & Ordering:**
- **Vendor Research**: Simulated DNA synthesis vendor comparison and selection
- **Synthesis Submission**: Submit sequences for DNA synthesis with pricing quotes

**Session & Workflow Management:**
- **Session Creation**: Create and manage analysis sessions
- **Toolbox Inventory**: Ask "what tools do you have?" to see all available tools, MCP tools, and CLI tools

### 🔄 **Session Management**
- **Persistent Sessions**: Track your workflow across multiple commands
- **History Tracking**: Complete audit trail of all operations
- **Context Preservation**: Maintain state between commands
- **Workflow Context**: Pass data between different analysis steps

### 🧾 Job Management (Async)
- **FastQC jobs on EMR** return immediately with a `job_id`
- **Core job endpoints**:
  - `GET /jobs/{job_id}` — status + metadata
  - `GET /jobs/{job_id}/results` — results pointer(s) for completed jobs
  - Additional helpers exist (e.g., logs/cancel/retry/copy-to-session) depending on deployment configuration

### 🎯 **User Experience**
- **Drag-and-Drop File Upload**: Upload FASTA/CSV files directly (no auto-population)
- **Command Mode Toggle**: Switch between natural language and structured commands
- **Real-time Feedback**: Immediate response and progress indicators
- **Responsive Design**: Works on desktop and mobile devices
- **Professional UI**: Clean layout with collapsible Examples; MCP Tools panel removed for clarity
- **Agent vs Submit**: Two buttons with clear tooltips; default dispatch uses heuristics but users can override

### 🔒 Usage Limits & Safety
- **Max upload size**: 10 MB per file (enforced server-side).
- **Max prompt length**: 5,000 tokens (estimated server-side).
- **Daily prompt cap**: 100 prompts per session/IP per day.
- In production, back these limits with Redis or a database for distributed enforcement.

## 🏗 Architecture Improvements (Version 2.0)

### Centralized Agent Tools Module
- **`backend/agent_tools.py`**: All agent tool definitions (`@tool` decorated functions) are now centralized in a single module
  - Improved maintainability and discoverability
  - Better separation of concerns between agent logic and tool implementations
  - Simplified imports and dependency management
  - 16+ bioinformatics tools organized and documented

### Refactored Agent System
- **Simplified `agent.py`**: Major codebase reduction (1400+ lines refactored)
  - Cleaner `CommandProcessor` class with dependency injection
  - Improved error handling and response formatting
  - Better testability and modularity

### Enhanced Testing Infrastructure
- **Evaluation Framework** (`tests/evals/`): Systematic evaluation of agent components
- **Deployed Backend Tests**: Integration tests for production deployments
- **Comprehensive Test Coverage**: Unit, integration, workflow, and evaluation tests

### Agent Specifications
- **`agents/` directory**: Structured agent specifications and prompts
  - Intent detector agent specification
  - Tool generator agent specification
  - Clear separation of agent roles and responsibilities

### Frontend Improvements
- **ThinkingIndicator Component**: Better user feedback during agent processing
- **Improved Component Organization**: Better structure and maintainability
- **Enhanced Styling**: Component-specific styles in dedicated directory

## 🗂 Project Structure

```
Helix.AI/
├── start.sh                 # 🚀 Single startup command
├── frontend/                # React frontend with natural language support
│   ├── src/
│   │   ├── components/     # React components
│   │   │   ├── PhylogeneticTree.tsx      # Interactive phylogenetic tree visualization
│   │   │   ├── PlasmidVisualizer.tsx     # Plasmid and vector visualization
│   │   │   ├── ThinkingIndicator.tsx     # Loading/thinking indicator component
│   │   │   ├── ExamplesPanel.tsx         # Example commands panel
│   │   │   ├── JobsPanel.tsx             # Job management panel
│   │   │   └── ...                       # Other UI components
│   │   ├── services/       # API services
│   │   │   └── mcpApi.ts   # MCP API service with session management
│   │   ├── styles/         # Component-specific styles
│   │   ├── utils/          # Command parser and utilities
│   │   └── App.tsx        # Main application with drag-and-drop
│   └── package.json
├── backend/                 # FastAPI + Enhanced MCP backend
│   ├── main_with_mcp.py   # Main FastAPI application with MCP integration
│   ├── agent.py           # LangChain agent (refactored, simplified)
│   ├── agent_tools.py     # 🆕 Centralized agent tool definitions (@tool functions)
│   ├── command_router.py  # Command handling and routing
│   ├── context_builder.py # Context building for agent execution
│   ├── execution_broker.py # Execution routing and policy management
│   ├── history_manager.py # Session and history management
│   ├── intent_classifier.py # Intent classification for command routing
│   ├── tool_generator_agent.py # Dynamic tool generation agent
│   ├── tool_inventory.py  # Tool discovery and inventory management
│   ├── tool_schemas.py    # Tool schema definitions
│   ├── job_manager.py     # Async job management (EMR, etc.)
│   ├── plan_ir.py         # Plan intermediate representation
│   ├── prompts/           # Agent prompt templates
│   └── requirements.txt
├── agents/                 # 🆕 Agent specification files
│   ├── intent-detector-agent.md  # Intent classification agent spec
│   └── tool-generator-agent.md   # Tool generation agent spec
├── tools/                   # Bioinformatics tool modules
│   ├── alignment.py        # Sequence alignment tools
│   ├── mutations.py        # Sequence mutation and variant generation
│   ├── phylogenetic_tree.py # Phylogenetic tree construction and analysis
│   ├── plasmid_visualizer.py # Plasmid and vector visualization
│   ├── dna_vendor_research.py # Simulated DNA synthesis vendor research
│   ├── sequence_selection.py # Smart sequence selection algorithms
│   ├── single_cell_analysis.py # Single-cell RNA-seq analysis
│   ├── bulk_rnaseq.py      # Bulk RNA-seq differential expression (DESeq2)
│   ├── quality_assessment.py # FastQC quality control
│   ├── read_merging.py     # Read merging and preprocessing
│   ├── synthesis_submission.py # DNA synthesis submission
│   ├── ncbi_tools.py       # NCBI sequence fetching
│   ├── uniprot_tools.py    # UniProt protein database queries
│   ├── go_tools.py         # Gene Ontology term lookup
│   └── r_scripts/          # R scripts for analysis workflows
├── tests/                   # Comprehensive test suite
│   ├── backend/            # Backend-specific tests
│   │   ├── test_agent_tools.py        # 🆕 Agent tools tests
│   │   ├── test_command_processor.py  # 🆕 Command processor tests
│   │   ├── test_execution_logs.py     # 🆕 Execution logging tests
│   │   └── ...                        # Other backend tests
│   ├── evals/              # 🆕 Evaluation framework
│   │   ├── cases/          # Test case files (JSONL format)
│   │   ├── test_eval_intent_classifier.py
│   │   └── test_eval_router_tool_mapping.py
│   ├── integration/        # End-to-end integration tests
│   │   ├── test_deployed_backend.py   # 🆕 Deployed backend tests
│   │   └── run_deployed_tests.sh      # 🆕 Test runner script
│   ├── workflows/          # Workflow integration tests
│   └── README.md          # Test documentation
├── infrastructure/         # AWS CDK infrastructure as code
│   ├── helix_infrastructure/
│   │   └── helix_stack.py  # CDK stack definition
│   └── app.py              # CDK app entry point
├── scripts/                # Deployment and utility scripts
│   ├── aws/                # AWS deployment scripts
│   │   ├── deploy.sh       # Main deployment script
│   │   ├── setup-*.sh      # Various setup scripts
│   │   └── ...             # Other AWS utilities
│   └── emr/                # EMR-specific scripts
├── docs/                    # Documentation
│   ├── demos/             # Demo and tutorial files
│   ├── reports/           # Test reports and analysis
│   ├── MCP_SERVER_README.md
│   ├── ENHANCED_MCP_README.md
│   ├── HISTORY_TRACKING.md
│   └── NATURAL_LANGUAGE_GUIDE.md
├── shared/                 # Shared utilities and models
├── data/                   # Data files and samples
│   ├── samples/            # Sample sequence files
│   └── rnaseq_demo/        # RNA-seq demo data
└── README.md              # This file
```

## 🚀 Quick Start

### Prerequisites
- Python 3.10+
- Node.js 16+
- npm or yarn
- [uv](https://github.com/astral-sh/uv) (optional, recommended for faster installs)
- AWS credentials configured (for FastQC EMR jobs)

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/your-username/Helix.AI.git
cd Helix.AI
   ```

2. **Install backend dependencies** (choose one):
   
   **Option A: Using uv (recommended - faster)**
   ```bash
   # Install uv if not already installed
   curl -LsSf https://astral.sh/uv/install.sh | sh
   
   # Install dependencies
   cd backend
   uv pip install -r requirements.txt
   cd ..
   ```
   
   **Option B: Using pip**
   ```bash
   cd backend
   pip install -r requirements.txt
   cd ..
   ```

3. **Install frontend dependencies**
   ```bash
   cd frontend
   npm install
   cd ..
   ```

4. **Configure environment variables** (optional)
   
   The startup script automatically sets the EMR cluster ID (`EMR_CLUSTER_ID=j-12QYDE51Q9LDP`) for FastQC jobs.
   
   If you need to **set up AWS infrastructure (EMR + EC2)** or use a different cluster/instance, we have dedicated docs + scripts:
   - **EMR cluster setup**: see `docs/EMR_SETUP_GUIDE.md` (script: `scripts/aws/setup-emr-cluster.sh`)
   - **EC2 execution setup (bioinformatics tools box)**: see `docs/EC2_EXECUTION_SETUP.md` (scripts: `scripts/aws/create-bioinformatics-ec2.sh`, `scripts/aws/setup-existing-instance.sh`)
   - **AWS scripts quick-start**: `scripts/aws/QUICK_START.md`
   
   To override defaults or set additional variables, create a `.env` file in the project root:
   ```bash
   # .env file
   EMR_CLUSTER_ID=j-12QYDE51Q9LDP
   OPENAI_API_KEY=your_openai_api_key_here
   DEEPSEEK_API_KEY=your_deepseek_api_key_here
   ```

5. **Start the unified system**
   ```bash
   ./start.sh
   ```

The application will be available at:
- **Frontend**: http://localhost:5173
- **Backend API**: http://localhost:8001

### Environment Variables

The startup script automatically sets the following environment variables:

- **EMR_CLUSTER_ID**: Set to `j-12QYDE51Q9LDP` for FastQC job submissions to AWS EMR
- **REDIS_URL**: Automatically configured if Redis is available
- **PYTHONPATH**: Configured to include the tools directory

Common optional variables:
- **HELIX_ASYNC_BYTES_THRESHOLD**: Bytes threshold for routing sync → async (default 100MB)
- **HELIX_MOCK_MODE**: Set `1` to run without real LLM/cloud dependencies (used by unit tests)

You can override these by setting them in a `.env` file in the project root or backend directory.

## 🧭 How Dispatch Works (Agent vs Submit)

- **Auto-dispatch defaults**
  - Question/exploration (“what is…”, “how should I…”, “design a workflow…”) → **Agent**.
  - Clear operation with required inputs present (align/trim/merge/mutate/visualize) → **Submit**.
- **Override anytime**: Click the other button to force the route.
- **Input validation**: Direct runs are blocked with actionable hints if required inputs are missing.
- **Agent output**: Only the final response content is shown (Markdown). Raw traces are suppressed.

Examples:
- “What is a sequence alignment?” → Agent
- “Align these sequences” (and sequences uploaded/pasted) → Submit
- “Plan a short-read preprocessing pipeline” → Agent (may propose trim → merge → QC → align)

## 📖 Usage Examples

### Basic Workflow
1. **Upload sequences**: Drag and drop a FASTA file
2. **Visualize tree**: Type "visualize the phylogenetic tree"
3. **Select representatives**: Type "select 10 representative sequences"
4. **Insert into plasmid**: Type "insert each of the sequences into pUC19"

### Advanced Commands
- `"align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATC"`
- `"perform multiple sequence alignment on the uploaded sequences"`
- `"select 5 sequences with the highest mutation rate"`
- `"research DNA synthesis vendors for 1000bp sequences"` (simulated data)
- `"show me the plasmid visualization with circular view"`

## ☁️ Deployment (AWS)

### Quick Start

**Fully automated deployment** is available with our deployment scripts:

1. **Set up deployment configuration**:
   ```bash
   cd scripts/aws
   ./setup-deployment.sh
   ```

2. **Deploy infrastructure** (first time only):
   ```bash
   cd infrastructure
   pip install -r requirements.txt
   cdk bootstrap
   cdk deploy
   ```

3. **Deploy application**:
   ```bash
   cd scripts/aws
   ./deploy.sh
   ```

### Architecture Overview

The deployment includes:
- **Backend**: FastAPI application running on ECS Fargate with ALB
- **Frontend**: React application hosted on S3 with CloudFront CDN
- **Container Registry**: ECR for Docker images
- **Infrastructure as Code**: AWS CDK for complete infrastructure automation

### Deployment Options

1. **Automated Scripts** (Recommended): Use `scripts/aws/deploy.sh` for complete automation
2. **Infrastructure as Code**: AWS CDK in `infrastructure/` directory
3. **CI/CD**: GitHub Actions workflow (`.github/workflows/deploy.yml`)

### Detailed Documentation

For comprehensive deployment instructions, see:
- **[AWS Deployment Guide](docs/AWS_DEPLOYMENT_GUIDE.md)**: Complete step-by-step deployment guide
- **[Infrastructure README](infrastructure/README.md)**: AWS CDK setup and configuration
- **[Deployment Scripts README](scripts/aws/README.md)**: Script usage and options

### GitHub Actions (CI/CD) Setup

The project includes automated CI/CD via GitHub Actions. To enable:

1. **Add GitHub Secrets**:
   - `AWS_REGION` (e.g., us-east-1)
   - `AWS_ACCOUNT_ID`
   - `AWS_OIDC_ROLE_ARN` (IAM role for GitHub OIDC)
   - `ECR_REPOSITORY` (e.g., helix-backend)
   - `S3_BUCKET` (S3 bucket for frontend)
   - `CLOUDFRONT_DISTRIBUTION_ID` (optional)
   - `VITE_API_BASE_URL` (frontend build-time API base URL)
   - `ECS_CLUSTER_NAME`, `ECS_SERVICE_NAME`, `ECS_TASK_DEFINITION_FAMILY` (optional - for auto ECS updates)

2. **Workflow**: `.github/workflows/deploy.yml` automatically:
   - Builds and pushes backend Docker image to ECR
   - Updates ECS service with new image (if configured)
   - Builds and deploys frontend to S3
   - Invalidates CloudFront cache (if configured)

   Triggers on: push to `main` branch or manual workflow dispatch

## 🧪 Testing

Run the comprehensive test suite:

```bash
# Unit tests (default; integration tests are deselected)
# Runs in mock mode (HELIX_MOCK_MODE=1) by default
pytest

# Opt-in integration tests (requires backend running locally and network access)
pytest -m integration

# Run evaluation framework tests
pytest tests/evals/

# Run deployed backend tests (requires deployed backend URL)
pytest tests/integration/test_deployed_backend.py

# Frontend tests
cd frontend && npm test

# Run specific test categories
pytest tests/backend/              # Backend unit tests only
pytest tests/workflows/            # Workflow integration tests
pytest tests/backend/test_agent_tools.py  # Agent tools tests
```

### Test Structure

- **Unit Tests** (`tests/backend/`): Fast, isolated tests that run in mock mode by default
  - Agent tools tests (`test_agent_tools.py`)
  - Command processor tests (`test_command_processor.py`)
  - Execution broker and routing tests
  - Intent classifier tests
  - Tool generator agent tests

- **Evaluation Framework** (`tests/evals/`): Systematic evaluation of agent components
  - Intent classification evaluation
  - Router tool mapping evaluation
  - Test cases in JSONL format for reproducibility

- **Integration Tests** (`tests/integration/`): End-to-end tests with real backend
  - Deployed backend functionality tests
  - Natural language command mapping tests
  - Core functionality validation

- **Workflow Tests** (`tests/workflows/`): Multi-step workflow validation
  - Alignment and phylogenetic tree workflows
  - RNA-seq analysis workflows
  - Variant analysis workflows

## 📚 Documentation

### System Documentation
- **[Bioinformatics Agent Plan](docs/BioinformaticsAgentPlan.md)**: System architecture, agent prompt, roadmap
- **[Natural Language Guide](docs/NATURAL_LANGUAGE_GUIDE.md)**: How to use natural language commands
- **[Development Guide](docs/DEVELOPMENT_GUIDE.md)**: Development setup and guidelines

### Feature Guides
- **[Phylogenetic Tree Guide](docs/PHYLOGENETIC_TREE_GUIDE.md)**: Tree visualization and analysis
- **[DNA Vendor Research](docs/DNA_VENDOR_RESEARCH_GUIDE.md)**: Simulated vendor comparison and selection
- **[EMR Setup Guide](docs/EMR_SETUP_GUIDE.md)**: AWS EMR cluster setup and configuration

### Agent Specifications
- **[Intent Detector Agent](agents/intent-detector-agent.md)**: Intent classification agent specification
- **[Tool Generator Agent](agents/tool-generator-agent.md)**: Dynamic tool generation agent specification

## 🤝 Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct and the process for submitting pull requests.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **ETE3**: Phylogenetic tree visualization
- **SeqViz**: Plasmid visualization
- **D3.js**: Interactive data visualization
- **FastAPI**: Modern web framework
- **React**: Frontend framework
- **BioPython**: Bioinformatics toolkit
- **Cursor.AI**: AI-based code generator and editor

## 📊 Status

- **Architecture**: ✅ Unified monolithic system (primary) - Refactored and improved in v2.0
- **Backend**: ✅ Running (FastAPI + Enhanced MCP + Centralized Agent Tools)
- **Frontend**: ✅ Running (React + TypeScript + Enhanced UI components)
- **Agent Tools**: ✅ 16+ tools centrally organized in `agent_tools.py`
- **Session Management**: ✅ File-based with optional Redis
- **Testing**: ✅ Comprehensive test suite (unit, integration, eval, workflow tests)
- **Documentation**: ✅ Complete documentation including agent specifications
- **Evaluation Framework**: ✅ Systematic evaluation infrastructure
- **Cloud Ready**: ✅ AWS deployment with CDK, ECS, EMR support
