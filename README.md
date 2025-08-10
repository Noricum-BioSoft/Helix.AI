# 🧬 Helix.AI - Bioinformatics AI Platform

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Node.js](https://img.shields.io/badge/Node.js-16+-green.svg)](https://nodejs.org/)
[![React](https://img.shields.io/badge/React-18+-blue.svg)](https://reactjs.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.100+-green.svg)](https://fastapi.tiangolo.com/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An AI-powered web application for managing biotechnology workflows via natural language commands, featuring interactive phylogenetic tree visualization, DNA synthesis vendor research, and comprehensive bioinformatics tools with session management and history tracking.

## 🚀 Features

### 🌳 **Interactive Phylogenetic Tree Visualization**
- **ETE3 Integration**: High-quality phylogenetic tree visualization with SVG rendering
- **D3.js Tree Rendering**: Interactive phylogenetic trees with zoom and pan
- **Newick Format Support**: Parse and display standard tree formats
- **Tree Statistics**: Display comprehensive tree metrics and relationships
- **Real-time Updates**: Trees update as you analyze new sequences
- **Clustering Analysis**: Select representative sequences from phylogenetic clusters

### 🔬 **DNA Synthesis Vendor Research**
- **Vendor Comparison**: Compare 7 major DNA synthesis vendors
- **Pricing Analysis**: Real-time pricing ranges and service comparisons
- **Service Matching**: Find vendors based on sequence length and quantity
- **Recommendations**: AI-powered vendor selection advice

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
- **Professional UI**: Clean 75/25 layout with optimized spacing and helpful tips

## 🗂 Project Structure

```
Helix.AI/
├── frontend/                 # React frontend with natural language support
│   ├── src/
│   │   ├── components/      # React components including PhylogeneticTree, PlasmidVisualizer
│   │   ├── services/        # MCP API service with session management
│   │   ├── utils/           # Command parser and utilities
│   │   └── App.tsx         # Main application with drag-and-drop
│   └── package.json
├── backend/                  # FastAPI + MCP server with session tracking
│   ├── main_with_mcp.py    # Main FastAPI application with MCP integration
│   ├── history_manager.py   # Session and history management
│   ├── command_router.py   # Command handling and routing
│   ├── agent.py            # LangChain agent with bioinformatics tools
│   └── requirements.txt
├── tools/                    # Bioinformatics tool modules
│   ├── mutations.py         # Sequence mutation and variant generation
│   ├── alignment.py         # Sequence alignment tools
│   ├── data_science.py      # Statistical analysis and visualization
│   ├── variant_selection.py # Smart variant selection algorithms
│   ├── phylogenetic_tree.py # Phylogenetic tree construction and analysis
│   ├── dna_vendor_research.py # DNA synthesis vendor research
│   ├── command_parser.py    # Natural language command parsing
│   ├── command_executor.py  # Command execution engine
│   ├── command_handler.py   # Combined parser and executor
│   └── plasmid_visualizer.py # Plasmid and vector visualization
├── tests/                    # Comprehensive test suite
│   ├── backend/             # Backend-specific tests
│   ├── frontend/            # Frontend-specific tests
│   ├── integration/         # End-to-end integration tests
│   └── README.md           # Test documentation
├── data/                     # Data files and samples
│   ├── samples/             # Sample sequence files
│   ├── phylogenetic/        # Phylogenetic tree datasets
│   └── README.md           # Data documentation
├── docs/                     # Documentation
│   ├── demos/              # Demo and tutorial files
│   ├── reports/            # Test reports and analysis
│   ├── MCP_SERVER_README.md
│   ├── ENHANCED_MCP_README.md
│   ├── HISTORY_TRACKING.md
│   └── NATURAL_LANGUAGE_GUIDE.md
├── shared/                  # Shared utilities and models
└── start.sh                 # Unified startup script
```

## 🚀 Quick Start

### Prerequisites
- Python 3.10+
- Node.js 16+
- npm or yarn

### Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/your-username/Helix.AI.git
cd Helix.AI
   ```

2. **Install backend dependencies**
   ```bash
   pip install -r requirements.txt
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
- `"research DNA synthesis vendors for 1000bp sequences"`
- `"show me the plasmid visualization with circular view"`

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

- **[Natural Language Guide](docs/NATURAL_LANGUAGE_GUIDE.md)**: How to use natural language commands
- **[Phylogenetic Tree Guide](docs/PHYLOGENETIC_TREE_GUIDE.md)**: Tree visualization and analysis
- **[DNA Vendor Research](docs/DNA_VENDOR_RESEARCH_GUIDE.md)**: Vendor comparison and selection
- **[MCP Server Documentation](docs/MCP_SERVER_README.md)**: Model Context Protocol integration
- **[Development Guide](docs/DEVELOPMENT_SUMMARY.md)**: Development setup and guidelines

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

## 📊 Status

- **Architecture**: ✅ Unified monolithic system
- **Backend**: ✅ Running (FastAPI + Enhanced MCP)
- **Frontend**: ✅ Running (React + TypeScript)
- **Session Management**: ✅ File-based with optional Redis
- **Testing**: ✅ Comprehensive test suite
- **Documentation**: ✅ Complete documentation
