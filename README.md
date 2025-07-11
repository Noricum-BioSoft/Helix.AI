# 🧬 DataBloom.AI - Bioinformatics AI Platform

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Node.js](https://img.shields.io/badge/Node.js-16+-green.svg)](https://nodejs.org/)
[![React](https://img.shields.io/badge/React-18+-blue.svg)](https://reactjs.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.100+-green.svg)](https://fastapi.tiangolo.com/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

An AI-powered web application for managing biotechnology workflows via natural language commands, featuring a Model Context Protocol (MCP) server for standardized bioinformatics operations with session management and history tracking.

## 🚀 Features

### 🤖 **Natural Language Commands**
- **Intelligent Command Parsing**: Understand commands like "from the sequence variants, pick 10 sequences randomly and output them"
- **Multi-step Workflows**: Chain commands together for complex bioinformatics pipelines
- **Context Awareness**: Maintains session history for continuous workflows

### 📊 **Interactive Visualizations**
- **Plotly Integration**: Dynamic charts for sequence analysis, mutations, and alignments
- **Real-time Updates**: Live visualizations that update as you work
- **Export Capabilities**: Save plots and results in multiple formats

### 🧰 **Comprehensive Tool Suite**
- **Sequence Alignment**: Multiple algorithms (ClustalW, Muscle, MAFFT)
- **Mutation Analysis**: Generate and analyze sequence variants
- **Data Science Tools**: Statistical analysis and feature engineering
- **Variant Selection**: Smart selection based on diversity, length, or custom criteria
- **DNA Synthesis Research**: Vendor research and testing options
- **Plasmid Visualization**: Interactive plasmid and vector visualization

### 🔄 **Session Management**
- **Persistent Sessions**: Track your workflow across multiple commands
- **History Tracking**: Complete audit trail of all operations
- **Context Preservation**: Maintain state between commands

### 🎯 **User Experience**
- **Drag-and-Drop File Upload**: Upload FASTA/CSV files directly
- **Command Mode Toggle**: Switch between natural language and structured commands
- **Real-time Feedback**: Immediate response and progress indicators
- **Responsive Design**: Works on desktop and mobile devices

## 🗂 Project Structure

```
DataBloom.AI/
├── frontend/                 # React frontend with natural language support
│   ├── src/
│   │   ├── components/      # React components including PlasmidVisualizer
│   │   ├── services/        # MCP API service with session management
│   │   ├── utils/           # Command parser and utilities
│   │   └── App.tsx         # Main application with drag-and-drop
│   └── package.json
├── backend/                  # FastAPI + MCP server with session tracking
│   ├── main_with_mcp.py    # FastAPI with natural language integration
│   ├── history_manager.py   # Session and history management
│   ├── agent.py            # Command handling and routing
│   └── requirements.txt
├── tools/                    # Bioinformatics tool modules
│   ├── mutations.py         # Sequence mutation and variant generation
│   ├── alignment.py         # Sequence alignment tools
│   ├── data_science.py      # Statistical analysis and visualization
│   ├── variant_selection.py # Smart variant selection algorithms
│   ├── command_parser.py    # Natural language command parsing
│   ├── command_executor.py  # Command execution engine
│   ├── command_handler.py   # Combined parser and executor
│   ├── plasmid_visualizer.py # Plasmid and vector visualization
│   └── dna_vendor_research.py # DNA synthesis vendor research
├── tests/                    # Test files and sample data
│   └── sample_files/        # Example FASTA and CSV files
├── docs/                     # Documentation
│   ├── MCP_SERVER_README.md
│   ├── ENHANCED_MCP_README.md
│   ├── HISTORY_TRACKING.md
│   └── NATURAL_LANGUAGE_GUIDE.md
└── start-app.sh             # Complete startup script
```

## 📦 Installation & Setup

### Prerequisites

- **Node.js** >= 16
- **Python** >= 3.10
- **Conda** (recommended for environment management)
- **Git** (for cloning the repository)

### 1. Clone and Setup Repository

```bash
# Clone the repository
git clone https://github.com/yourusername/DataBloom.AI.git
cd DataBloom.AI

# Make startup script executable
chmod +x start-app.sh
```

### 2. Backend Setup

```bash
# Navigate to backend directory
cd backend

# Create conda environment (recommended)
conda create -n databloom python=3.10
conda activate databloom

# Or use virtual environment
python -m venv venv
source venv/bin/activate  # macOS/Linux
# venv\Scripts\activate   # Windows

# Install Python dependencies
pip install -r requirements.txt

# Install additional bioinformatics tools
conda install -c bioconda clustalo muscle mafft
# Or on macOS:
# brew install clustalo muscle mafft
```

### 3. Frontend Setup

```bash
# Navigate to frontend directory
cd frontend

# Install Node.js dependencies
npm install
```

### 4. Environment Configuration

Create a `.env` file in the backend directory:

```bash
cd backend
touch .env
```

Add the following to `.env`:

```env
# API Keys (if using external services)
OPENAI_API_KEY=your-openai-api-key
DEEPSEEK_API_KEY=your-deepseek-api-key

# Server Configuration
MCP_LOG_LEVEL=INFO
MCP_VALIDATION_STRICT=true

# Session Management
SESSION_TIMEOUT=3600
MAX_SESSION_SIZE=1000

# Database (if using)
DATABASE_URL=sqlite:///./data.db
```

## 🚀 Running the Application

### Option 1: Complete Startup (Recommended)

Use the provided startup script to run everything at once:

```bash
# From the project root
./start-app.sh
```

This script will:
- ✅ Check prerequisites and dependencies
- ✅ Start the backend MCP server with session management
- ✅ Start the frontend development server
- ✅ Wait for services to be ready
- ✅ Open the application in your browser

### Option 2: Manual Startup

#### Start Backend with Session Management

```bash
# Navigate to backend
cd backend

# Activate conda environment
conda activate databloom

# Start the FastAPI server with MCP integration
python main_with_mcp.py
```

#### Start Frontend

```bash
# Navigate to frontend
cd frontend

# Start development server
npm run dev
```

## 🌐 Accessing the Application

Once running, access the application at:

- **Frontend**: http://localhost:5173
- **Backend API**: http://localhost:8001
- **API Documentation**: http://localhost:8001/docs
- **Health Check**: http://localhost:8001/health

## 🧪 Testing the Natural Language Commands

### Complete Workflow Example

```bash
# 1. Create a session and generate mutations
curl -X POST http://localhost:8001/session/create \
  -H "Content-Type: application/json" \
  -d '{}'

# 2. Generate mutations (replace SESSION_ID with actual session ID)
curl -X POST http://localhost:8001/mcp/mutate-sequence \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ACTGTTGAC",
    "num_variants": 10,
    "session_id": "SESSION_ID"
  }'

# 3. Use natural language to select variants
curl -X POST http://localhost:8001/mcp/handle-natural-command \
  -H "Content-Type: application/json" \
  -d '{
    "command": "from the sequence variants, pick 5 sequences randomly and output them",
    "session_id": "SESSION_ID"
  }'
```

### Frontend Testing

1. **Open the web interface** at http://localhost:5173
2. **Switch to "Natural Language" mode**
3. **Run these example commands**:

```
# Step 1: Generate variants
mutate sequence ACTGTTGAC with 10 variants

# Step 2: Select variants using natural language
from the sequence variants, pick 5 sequences randomly and output them

# Step 3: Analyze the selected variants
analyze the selected variants and show me the most diverse ones

# Step 4: Create a multi-step workflow
mutate this sequence, then align the variants and pick the best ones
```

## 🧠 Example Commands

### Natural Language Commands

| Command | Action | Description |
|---------|--------|-------------|
| `"from the sequence variants, pick 10 sequences randomly and output them"` | Variant Selection | Selects 10 random variants from previous mutation results |
| `"select 5 sequences with the highest mutation rate"` | Smart Selection | Uses custom criteria for variant selection |
| `"analyze the alignment and show me the most conserved regions"` | Analysis | Performs statistical analysis on alignment data |
| `"mutate this sequence, then align the variants and pick the best ones"` | Multi-step Workflow | Chains multiple operations together |
| `"I want to order these sequences from a DNA synthesis vendor"` | Vendor Research | Research DNA synthesis vendors and pricing |
| `"what testing options are available for my sequences?"` | Testing Research | Find validation and testing options |

### Structured Commands

| Command | Action | MCP Tool Used |
|---------|--------|---------------|
| `"align sequences ACTGTTGAC ACTGCATCC"` | Sequence alignment | `sequence_alignment` |
| `"mutate sequence ACTGTTGAC with 10 variants"` | Generate mutations | `mutate_sequence` |
| `"analyze sequence data for phylogeny"` | Sequence analysis | `analyze_sequence_data` |
| `"visualize alignment in PNG format"` | Create visualization | `visualize_alignment` |
| `"express ATGCGATCG in pTet vector"` | Plasmid visualization | `plasmid_visualization` |

### File Upload Examples

1. **Upload FASTA file** and choose action:
   - Align sequences
   - Mutate sequences
   - Analyze data

2. **Upload CSV file** for data analysis:
   - Statistical analysis
   - Feature engineering
   - Visualization

## 🔧 API Endpoints

### Natural Language Commands

```http
POST /mcp/handle-natural-command
POST /mcp/parse-command
POST /mcp/execute-command
```

### Session Management

```http
POST /session/create
GET /session/{session_id}
```

### Traditional MCP Endpoints

```http
POST /mcp/sequence-alignment
POST /mcp/mutate-sequence
POST /mcp/analyze-sequence-data
POST /mcp/visualize-alignment
POST /mcp/select-variants
POST /mcp/plasmid-visualization
GET /mcp/tools
```

### Health and Status

```http
GET /health
GET /execute
```

## 🛠️ Development

### Adding New Natural Language Commands

1. **Update Command Parser** in `tools/command_parser.py`:

```python
def parse_command_raw(command: str, session_id: str = None) -> Dict[str, Any]:
    # Add new command patterns
    if "new command pattern" in command.lower():
        return {
            "action": "new_action",
            "tool": "new_tool",
            "parameters": {"param": "value"},
            "session_id": session_id
        }
```

2. **Add Tool Implementation** in `tools/new_tool.py`:

```python
def run_new_tool_raw(param: str) -> Dict[str, Any]:
    """Execute the new tool with raw parameters."""
    # Implementation here
    return {
        "status": "success",
        "result": {"data": "processed_data"}
    }
```

3. **Update Frontend** in `frontend/src/App.tsx`:

```typescript
// Add to renderOutput function
{actualResult.new_data && (
  <div className="bg-light p-3 border rounded mb-3">
    <h5>New Tool Results</h5>
    <pre>{JSON.stringify(actualResult.new_data, null, 2)}</pre>
  </div>
)}
```

### Testing New Features

```bash
# Test command parsing
python -c "
import sys; sys.path.append('../tools')
from command_parser import parse_command_raw
result = parse_command_raw('your new command')
print(result)
"

# Test complete workflow
python test_command_handling.py

# Test API endpoints
curl -X POST http://localhost:8001/mcp/handle-natural-command \
  -H "Content-Type: application/json" \
  -d '{"command": "your new command", "session_id": "test"}'
```

## 🐛 Troubleshooting

### Common Issues

#### Natural Language Commands Not Working

1. **Check Session Management**:
   ```bash
   # Verify session creation
   curl -X POST http://localhost:8001/session/create
   
   # Check session history
   curl http://localhost:8001/session/{session_id}
   ```

2. **Check Command Parsing**:
   ```bash
   # Test command parser directly
   python -c "
   import sys; sys.path.append('../tools')
   from command_parser import parse_command_raw
   result = parse_command_raw('your command')
   print(result)
   "
   ```

#### Backend Issues

1. **Import Errors**:
   ```bash
   # Ensure conda environment is activated
   conda activate databloom
   
   # Reinstall dependencies
   pip install -r requirements.txt
   ```

2. **MCP Server Not Starting**:
   ```bash
   # Check Python path
   export PYTHONPATH="."
   
   # Check dependencies
   python -c "import plotly, pandas, biopython"
   ```

#### Frontend Issues

1. **API Connection Failed**:
   ```bash
   # Check backend is running
   curl http://localhost:8001/health
   
   # Check CORS settings
   # Ensure backend allows frontend origin
   ```

2. **Natural Language Mode Not Working**:
   ```bash
   # Check browser console for errors
   # Verify session ID is being set
   # Check network tab for API calls
   ```

### Debugging

1. **Enable Debug Logging**:
   ```bash
   export MCP_LOG_LEVEL="DEBUG"
   ```

2. **Check Session Logs**:
   ```bash
   # Check session files
   ls -la backend/sessions/
   
   # Check session content
   cat backend/sessions/{session_id}.json
   ```

3. **Test Individual Components**:
   ```bash
   # Test command parsing
   python test_command_handling.py
   
   # Test API
   curl -X POST http://localhost:8001/mcp/handle-natural-command \
     -H "Content-Type: application/json" \
     -d '{"command": "test command"}'
   ```

## 📚 Documentation

- **Backend MCP**: `backend/MCP_SERVER_README.md`
- **Enhanced MCP**: `backend/ENHANCED_MCP_README.md`
- **History Tracking**: `backend/HISTORY_TRACKING.md`
- **Natural Language Guide**: `NATURAL_LANGUAGE_GUIDE.md`
- **Complete Workflow Example**: `COMPLETE_WORKFLOW_EXAMPLE.md`
- **API Documentation**: http://localhost:8001/docs

## 🤝 Contributing

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/natural-language`
3. **Add your changes**: Follow the development guidelines above
4. **Test thoroughly**: Run all tests and verify functionality
5. **Submit a pull request**: Include detailed description of changes

### Development Guidelines

- **Follow existing patterns** for natural language command development
- **Add comprehensive tests** for new features
- **Update documentation** for new commands
- **Validate inputs** using existing validation functions
- **Include error handling** for all new commands
- **Maintain session compatibility** for multi-step workflows

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details.

## 🆘 Support

- **Issues**: Create an issue on GitHub
- **Documentation**: Check the README files in each directory
- **API Docs**: Visit http://localhost:8001/docs when running
- **Testing**: Use the provided test scripts to verify functionality
- **Examples**: See `tests/sample_files/` for example data files

## 🙏 Acknowledgments

- **BioPython** for bioinformatics tools
- **Plotly** for interactive visualizations
- **FastAPI** for the backend framework
- **React** for the frontend framework
- **MCP** for the Model Context Protocol implementation
