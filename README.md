# 🧬 Helix.AI - Bioinformatics AI Platform

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Node.js](https://img.shields.io/badge/Node.js-16+-green.svg)](https://nodejs.org/)
[![React](https://img.shields.io/badge/React-18+-blue.svg)](https://reactjs.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.100+-green.svg)](https://fastapi.tiangolo.com/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An AI-powered web application for managing biotechnology workflows via natural language commands, featuring interactive phylogenetic tree visualization, simulated DNA synthesis vendor research, and comprehensive bioinformatics tools with session management and history tracking.

## 🚀 Features

### 🧭 Tool‑first Dispatch (not just an LLM wrapper)
- **Hybrid routing (auto + override)**: Prompts are routed to either:
  - **Direct tool execution (Submit)** when required inputs are present and the intent matches a known operation.
  - **Agent planning (Agent)** for open-ended questions, missing inputs, or multi‑step workflows.
- **Validation before execution**: Required parameters (e.g., sequences for alignment) are checked up front.
- **Provenance & reproducibility**: Each run records tool name, version, parameters, inputs, and outputs.
- **Explainability**: When using the agent, the final answer is rendered in Markdown; intermediate planner chatter is hidden.

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
- **Sequence Alignment**: Multiple algorithms (ClustalW, Muscle, MAFFT)
- **Mutation Analysis**: Generate and analyze sequence variants
- **Data Science Tools**: Statistical analysis and feature engineering
- **Variant Selection**: Smart selection based on diversity, length, or custom criteria
- **Plasmid Visualization**: Interactive plasmid and vector visualization with circular/linear views
- **Clustering Analysis**: Hierarchical clustering with representative sequence selection

### 🔄 **Session Management**
- **Persistent Sessions**: Track your workflow across multiple commands
- **History Tracking**: Complete audit trail of all operations
- **Context Preservation**: Maintain state between commands
- **Workflow Context**: Pass data between different analysis steps

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

## 🗂 Project Structure

```
Helix.AI/
├── start.sh                 # 🚀 Single startup command
├── frontend/                # React frontend with natural language support
│   ├── src/
│   │   ├── components/     # React components including PhylogeneticTree, PlasmidVisualizer
│   │   ├── services/       # MCP API service with session management
│   │   ├── utils/          # Command parser and utilities
│   │   └── App.tsx        # Main application with drag-and-drop
│   └── package.json
├── backend/                 # FastAPI + Enhanced MCP backend
│   ├── main_with_mcp.py   # Main FastAPI application with MCP integration
│   ├── history_manager.py  # Session and history management
│   ├── command_router.py  # Command handling and routing
│   ├── agent.py           # LangChain agent with bioinformatics tools
│   └── requirements.txt
├── tools/                   # Bioinformatics tool modules
│   ├── mutations.py        # Sequence mutation and variant generation
│   ├── alignment.py        # Sequence alignment tools
│   ├── data_science.py     # Statistical analysis and visualization
│   ├── variant_selection.py # Smart variant selection algorithms
│   ├── phylogenetic_tree.py # Phylogenetic tree construction and analysis
│   ├── dna_vendor_research.py # Simulated DNA synthesis vendor research
│   ├── command_parser.py   # Natural language command parsing
│   ├── command_executor.py # Command execution engine
│   ├── command_handler.py  # Combined parser and executor
│   └── plasmid_visualizer.py # Plasmid and vector visualization
├── tests/                   # Comprehensive test suite
│   ├── backend/            # Backend-specific tests
│   ├── frontend/           # Frontend-specific tests
│   ├── integration/        # End-to-end integration tests
│   └── README.md          # Test documentation
├── data/                    # Data files and samples
│   ├── samples/            # Sample sequence files
│   ├── phylogenetic/       # Phylogenetic tree datasets
│   └── README.md          # Data documentation
├── docs/                    # Documentation
│   ├── demos/             # Demo and tutorial files
│   ├── reports/           # Test reports and analysis
│   ├── MCP_SERVER_README.md
│   ├── ENHANCED_MCP_README.md
│   ├── HISTORY_TRACKING.md
│   └── NATURAL_LANGUAGE_GUIDE.md
├── shared/                 # Shared utilities and models
└── README.md              # This file
```

## 🚀 Quick Start

### Prerequisites
- Python 3.10+
- Node.js 16+
- npm or yarn
- [uv](https://github.com/astral-sh/uv) (optional, recommended for faster installs)

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

4. **Start the unified system**
   ```bash
   ./start.sh
   ```

The application will be available at:
- **Frontend**: http://localhost:5173
- **Backend API**: http://localhost:8001

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
# Backend tests
python -m pytest tests/backend/

# Frontend tests
cd frontend && npm test

# Integration tests
python -m pytest tests/integration/
```

## 📚 Documentation

- **[Bioinformatics Agent Plan](docs/BioinformaticsAgentPlan.md)**: System architecture, agent prompt, roadmap
- **[Natural Language Guide](docs/NATURAL_LANGUAGE_GUIDE.md)**: How to use natural language commands
- **[Phylogenetic Tree Guide](docs/PHYLOGENETIC_TREE_GUIDE.md)**: Tree visualization and analysis
- **[DNA Vendor Research](docs/DNA_VENDOR_RESEARCH_GUIDE.md)**: Simulated vendor comparison and selection
- **[Development Guide](docs/DEVELOPMENT_GUIDE.md)**: Development setup and guidelines

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

- **Architecture**: ✅ Unified monolithic system (primary)
- **Backend**: ✅ Running (FastAPI + Enhanced MCP)
- **Frontend**: ✅ Running (React + TypeScript)
- **Session Management**: ✅ File-based with optional Redis
- **Testing**: ✅ Comprehensive test suite
- **Documentation**: ✅ Complete documentation
- **Cloud Ready**: 🔄 Microservices option for future scaling
