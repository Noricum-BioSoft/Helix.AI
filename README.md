# 🧬 DataBloom.AI - Bioinformatics MCP Server

An AI-powered web application for managing biotechnology workflows via natural language commands, featuring a Model Context Protocol (MCP) server for standardized bioinformatics operations.

---

## 🚀 Features

- ✍️ **Prompt-Based Interface** — Users enter commands like:
  - "Align the given sequences"
  - "Mutate this sequence and generate 96 variants"
  - "Calculate sequence statistics"
  - "Generate reverse complement"
- 🤖 **MCP-Powered Backend** — Standardized bioinformatics tools via Model Context Protocol
- 📊 **Interactive Visualizations** — Supports graphical results (e.g., Plotly plots for alignments)
- 🧰 **Modular Tool System** — Easily extend with new bioinformatics tools
- 🔄 **React Frontend** — Built with Bootstrap for a responsive, clean UI
- ⚡ **FastAPI Backend** — Lightweight and scalable RESTful API
- 🔍 **Enhanced Validation** — Comprehensive input validation and error handling

---

## 🗂 Project Structure

```
DataBloom.AI/
├── frontend/                 # React frontend with MCP integration
│   ├── src/
│   │   ├── services/        # MCP API service
│   │   ├── utils/           # Command parser
│   │   └── App.tsx         # Main application
│   └── package.json
├── backend/                  # FastAPI + MCP server
│   ├── mcp_server.py        # Basic MCP server
│   ├── mcp_server_enhanced.py # Enhanced MCP server
│   ├── main_with_mcp.py    # FastAPI with MCP integration
│   ├── tools/               # Bioinformatics tools
│   └── requirements.txt
├── tools/                    # Bioinformatics tool modules
└── start-app.sh             # Complete startup script
```

---

## 📦 Installation & Setup

### Prerequisites

- **Node.js** >= 16
- **Python** >= 3.9
- **Git** (for cloning the repository)

### 1. Clone and Setup Repository

```bash
# Clone the repository
git clone <repository-url>
cd DataBloom.AI

# Make startup script executable
chmod +x start-app.sh
```

### 2. Backend Setup

```bash
# Navigate to backend directory
cd backend

# Create virtual environment
python -m venv venv

# Activate virtual environment
# On macOS/Linux:
source venv/bin/activate
# On Windows:
# venv\Scripts\activate

# Install Python dependencies
pip install -r requirements.txt

# Install bioinformatics tools (optional but recommended)
# For Clustal Omega (sequence alignment)
# macOS:
brew install clustalo
# Ubuntu/Debian:
# sudo apt-get install clustalo
```

### 3. Frontend Setup

```bash
# Navigate to frontend directory
cd frontend

# Install Node.js dependencies
npm install

# Build the project (optional)
npm run build
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

# Database (if using)
DATABASE_URL=sqlite:///./data.db
```

---

## 🚀 Running the Application

### Option 1: Complete Startup (Recommended)

Use the provided startup script to run everything at once:

```bash
# From the project root
./start-app.sh
```

This script will:
- ✅ Check prerequisites
- ✅ Install dependencies
- ✅ Start the backend MCP server
- ✅ Start the frontend development server
- ✅ Wait for services to be ready
- ✅ Open the application in your browser

### Option 2: Manual Startup

#### Start Backend MCP Server

```bash
# Navigate to backend
cd backend

# Activate virtual environment
source venv/bin/activate  # macOS/Linux
# venv\Scripts\activate   # Windows

# Start the enhanced MCP server
python mcp_server_enhanced.py

# Or start the basic MCP server
python mcp_server.py

# Or start the FastAPI server with MCP integration
python main_with_mcp.py
```

#### Start Frontend

```bash
# Navigate to frontend
cd frontend

# Start development server
npm run dev
```

### Option 3: Individual Component Testing

#### Test Backend MCP Server

```bash
cd backend

# Test the enhanced MCP server
python test_enhanced_mcp.py

# Test basic integration
python test_mcp_integration.py

# Test API endpoints
python -m pytest tests/  # if tests exist
```

#### Test Frontend Integration

```bash
cd frontend

# Test API integration
node test-integration.js
```

---

## 🌐 Accessing the Application

Once running, access the application at:

- **Frontend**: http://localhost:5173
- **Backend API**: http://localhost:8001
- **API Documentation**: http://localhost:8001/docs
- **Health Check**: http://localhost:8001/health

---

## 🧪 Testing the MCP Server

### Run Comprehensive Tests

```bash
cd backend
python test_enhanced_mcp.py
```

### Test API Endpoints

```bash
# Test health endpoint
curl http://localhost:8001/health

# Test MCP tools listing
curl http://localhost:8001/mcp/tools

# Test sequence alignment
curl -X POST http://localhost:8001/mcp/sequence-alignment \
  -H "Content-Type: application/json" \
  -d '{
    "sequences": ">seq1\nACTGTTGAC\n>seq2\nACTGCATCC",
    "algorithm": "clustal"
  }'

# Test sequence statistics
curl -X POST http://localhost:8001/mcp/sequence-statistics \
  -H "Content-Type: application/json" \
  -d '{
    "sequence": "ACTGTTGAC",
    "include_composition": true
  }'
```

---

## 🧠 Example Commands

### Frontend Interface Commands

| Command | Action | MCP Tool Used |
|---------|--------|---------------|
| `"align sequences ACTGTTGAC ACTGCATCC"` | Sequence alignment | `sequence_alignment` |
| `"mutate sequence ACTGTTGAC with 10 variants"` | Generate mutations | `mutate_sequence` |
| `"calculate statistics for ACTGTTGAC"` | Sequence analysis | `sequence_statistics` |
| `"reverse complement of ACTGTTGAC"` | DNA manipulation | `reverse_complement` |
| `"validate sequence ACTGTTGAC"` | Input validation | `sequence_validation` |

### API Endpoints

#### MCP-Specific Endpoints

```http
POST /mcp/sequence-alignment
POST /mcp/mutate-sequence
POST /mcp/analyze-sequence-data
POST /mcp/visualize-alignment
POST /mcp/sequence-statistics
POST /mcp/reverse-complement
POST /mcp/sequence-validation
GET /mcp/tools
```

#### Legacy Endpoints

```http
POST /execute  # General command execution
GET /health    # Health check
```

---

## 🔧 MCP Server Configuration

### Enhanced MCP Server Features

The enhanced MCP server (`mcp_server_enhanced.py`) includes:

- **Input Validation**: Comprehensive validation of DNA/RNA sequences
- **Error Handling**: Structured error responses with timestamps
- **Performance Tracking**: Execution time measurement
- **Additional Tools**: Sequence statistics, reverse complement, validation

### Configuration Files

- `mcp_config.json` - Basic MCP server configuration
- `mcp_config_enhanced.json` - Enhanced configuration with validation settings

### Environment Variables

```bash
export MCP_LOG_LEVEL="INFO"
export MCP_VALIDATION_STRICT="true"
export PYTHONPATH="."
```

---

## 🛠️ Development

### Adding New Bioinformatics Tools

1. **Create Tool Function** in `backend/tools/`:

```python
# backend/tools/new_tool.py
from langchain.agents import tool

@tool
def run_new_tool(sequence: str):
    """Description of the new tool."""
    return {
        "text": "Tool executed successfully.",
        "result": {"data": "processed_data"}
    }
```

2. **Add to MCP Server** in `mcp_server_enhanced.py`:

```python
# Add tool to handle_list_tools()
Tool(
    name="new_tool",
    description="Description of the new tool",
    inputSchema={
        "type": "object",
        "properties": {
            "sequence": {
                "type": "string",
                "description": "Input sequence"
            }
        },
        "required": ["sequence"]
    }
)

# Add handler function
async def handle_new_tool(arguments: Dict[str, Any]) -> Dict[str, Any]:
    sequence = arguments.get("sequence", "")
    if not validate_sequence(sequence):
        raise ValidationError("Invalid sequence")
    
    result = run_new_tool(sequence)
    return {
        "status": "success",
        "tool": "new_tool",
        "result": result
    }

# Add to tool router
elif name == "new_tool":
    result = await handle_new_tool(arguments)
```

3. **Update Frontend** in `frontend/src/utils/commandParser.ts`:

```typescript
// Add keywords
private static newToolKeywords = [
    'new tool', 'process sequence'
];

// Add parsing method
private static parseNewToolCommand(command: string): ParsedCommand {
    // Implementation
}
```

### Testing New Features

```bash
# Test backend tools
cd backend
python test_enhanced_mcp.py

# Test frontend integration
cd frontend
npm test

# Test API endpoints
curl -X POST http://localhost:8001/mcp/new-tool \
  -H "Content-Type: application/json" \
  -d '{"sequence": "ACTGTTGAC"}'
```

---

## 🐛 Troubleshooting

### Common Issues

#### Backend Issues

1. **Import Errors**:
   ```bash
   # Ensure virtual environment is activated
   source venv/bin/activate
   
   # Reinstall dependencies
   pip install -r requirements.txt
   ```

2. **MCP Server Not Starting**:
   ```bash
   # Check Python path
   export PYTHONPATH="."
   
   # Check dependencies
   python -c "import mcp"
   ```

3. **Tool Not Found**:
   ```bash
   # Check tool registration
   curl http://localhost:8001/mcp/tools
   
   # Check logs
   tail -f logs/mcp_server.log
   ```

#### Frontend Issues

1. **API Connection Failed**:
   ```bash
   # Check backend is running
   curl http://localhost:8001/health
   
   # Check CORS settings
   # Ensure backend allows frontend origin
   ```

2. **Build Errors**:
   ```bash
   # Clear node modules
   rm -rf node_modules package-lock.json
   npm install
   
   # Check TypeScript
   npm run build
   ```

### Debugging

1. **Enable Debug Logging**:
   ```bash
   export MCP_LOG_LEVEL="DEBUG"
   ```

2. **Check Logs**:
   ```bash
   # Backend logs
   tail -f logs/mcp_server.log
   
   # Frontend logs
   # Check browser developer tools
   ```

3. **Test Individual Components**:
   ```bash
   # Test MCP server
   python test_enhanced_mcp.py
   
   # Test API
   node test-integration.js
   ```

---

## 📚 Documentation

- **Backend MCP**: `backend/MCP_README.md`
- **Enhanced MCP**: `backend/ENHANCED_MCP_README.md`
- **Frontend**: `frontend/README.md`
- **API Documentation**: http://localhost:8001/docs

---

## 🤝 Contributing

1. **Fork the repository**
2. **Create a feature branch**: `git checkout -b feature/new-tool`
3. **Add your changes**: Follow the development guidelines above
4. **Test thoroughly**: Run all tests and verify functionality
5. **Submit a pull request**: Include detailed description of changes

### Development Guidelines

- **Follow existing patterns** for tool development
- **Add comprehensive tests** for new features
- **Update documentation** for new tools
- **Validate inputs** using existing validation functions
- **Include error handling** for all new tools

---

## 📄 License

MIT License - see LICENSE file for details.

---

## 🆘 Support

- **Issues**: Create an issue on GitHub
- **Documentation**: Check the README files in each directory
- **API Docs**: Visit http://localhost:8001/docs when running
- **Testing**: Use the provided test scripts to verify functionality
