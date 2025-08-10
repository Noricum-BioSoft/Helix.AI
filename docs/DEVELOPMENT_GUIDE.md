# 🛠️ Helix.AI Development Guide

This guide provides comprehensive information for developers contributing to the Helix.AI bioinformatics platform.

## 🚀 Development Setup

### Prerequisites

- **Python 3.10+**: For backend development
- **Node.js 16+**: For frontend development
- **Git**: Version control
- **Docker** (optional): For containerized development

### Environment Setup

1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-username/Helix.AI.git
cd Helix.AI
   ```

2. **Set up Python environment**:
   ```bash
   # Create virtual environment
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   
   # Install dependencies
   pip install -r requirements.txt
   ```

3. **Set up Node.js environment**:
   ```bash
   cd frontend
   npm install
   cd ..
   ```

4. **Set up development environment**:
   ```bash
   # Make startup script executable
   chmod +x start.sh
   
   # Start the unified system
   ./start.sh
   ```

## 🏗️ Project Architecture

### Backend Architecture

```
backend/
├── main_with_mcp.py          # FastAPI application with MCP integration
├── command_router.py         # Natural language command routing
├── history_manager.py        # Session and history management
├── agent.py                  # LangChain agent with bioinformatics tools
├── utils.py                  # Shared utilities
└── requirements.txt          # Python dependencies
```

### Frontend Architecture

```
frontend/src/
├── App.tsx                   # Main application component
├── components/               # React components
│   ├── PhylogeneticTree.tsx  # Tree visualization
│   ├── PlasmidVisualizer.tsx # Plasmid visualization
│   └── DirectedEvolutionDemo.tsx
├── services/                 # API services
│   └── mcpApi.ts            # MCP API client
├── utils/                    # Utilities
│   └── commandParser.ts     # Command parsing utilities
└── main.tsx                 # Application entry point
```

### Tools Architecture

```
tools/
├── phylogenetic_tree.py      # Phylogenetic analysis
├── plasmid_visualizer.py     # Plasmid visualization
├── alignment.py              # Sequence alignment
├── mutations.py              # Sequence mutations
├── data_science.py          # Statistical analysis
├── variant_selection.py     # Variant selection
├── dna_vendor_research.py   # Vendor research
├── command_parser.py        # Natural language parsing
├── command_executor.py      # Command execution
└── command_handler.py       # Combined parser/executor
```

## 🔧 Development Workflow

### 1. **Feature Development**

```bash
# Create feature branch
git checkout -b feature/new-bioinformatics-tool

# Make changes
# ... implement feature ...

# Run tests
python -m pytest tests/backend/ -v
cd frontend && npm test

# Commit changes
git add .
git commit -m "feat: add new bioinformatics tool"

# Push and create PR
git push origin feature/new-bioinformatics-tool
```

### 2. **Backend Development**

#### Adding New Tools

1. **Create tool module** in `tools/`:
   ```python
   # tools/new_tool.py
   from typing import Dict, Any
   
   def run_new_tool_raw(param: str) -> Dict[str, Any]:
       """Execute the new tool with raw parameters."""
       # Implementation here
       return {
           "status": "success",
           "result": {"data": "processed_data"}
       }
   ```

2. **Add MCP endpoint** in `backend/main_with_mcp.py`:
   ```python
   @app.post("/mcp/new-tool")
   async def new_tool_mcp(req: NewToolRequest):
       """Execute new tool via MCP."""
       try:
           result = run_new_tool_raw(req.param)
           return CustomJSONResponse({
               "success": True,
               "result": result,
               "session_id": req.session_id
           })
       except Exception as e:
           return CustomJSONResponse({
               "success": False,
               "error": str(e),
               "session_id": req.session_id
           })
   ```

3. **Update command router** in `backend/command_router.py`:
   ```python
   # Add routing logic
   if any(phrase in command_lower for phrase in ['new tool', 'new analysis']):
       return 'new_tool', self._extract_parameters(command, 'new_tool', session_context)
   ```

4. **Add tests** in `tests/backend/`:
   ```python
   def test_new_tool_functionality():
       result = run_new_tool_raw("test_param")
       assert result["status"] == "success"
   ```

#### Adding New API Endpoints

1. **Define request/response models**:
   ```python
   from pydantic import BaseModel
   
   class NewToolRequest(BaseModel):
       param: str
       session_id: str = None
   ```

2. **Add endpoint**:
   ```python
   @app.post("/api/new-tool")
   async def new_tool_api(req: NewToolRequest):
       # Implementation
       pass
   ```

3. **Add tests**:
   ```python
   def test_new_tool_api():
       client = TestClient(app)
       response = client.post("/api/new-tool", json={"param": "test"})
       assert response.status_code == 200
   ```

### 3. **Frontend Development**

#### Adding New Components

1. **Create component** in `frontend/src/components/`:
   ```typescript
   // NewTool.tsx
   import React from 'react';
   
   interface NewToolProps {
       data: any;
   }
   
   export const NewTool: React.FC<NewToolProps> = ({ data }) => {
       return (
           <div className="card">
               <div className="card-header">
                   <h5>New Tool Results</h5>
               </div>
               <div className="card-body">
                   {/* Component content */}
               </div>
           </div>
       );
   };
   ```

2. **Add to main App** in `frontend/src/App.tsx`:
   ```typescript
   import { NewTool } from './components/NewTool';
   
   // In renderOutput function
   if (actualResult.new_tool_data) {
       return <NewTool data={actualResult.new_tool_data} />;
   }
   ```

3. **Add tests**:
   ```typescript
   // NewTool.test.tsx
   import { render, screen } from '@testing-library/react';
   import { NewTool } from './NewTool';
   
   test('renders new tool component', () => {
       render(<NewTool data={{}} />);
       expect(screen.getByText('New Tool Results')).toBeInTheDocument();
   });
   ```

#### Adding New API Services

1. **Add service method** in `frontend/src/services/mcpApi.ts`:
   ```typescript
   export const mcpApi = {
       // ... existing methods ...
       
       async newTool(param: string, sessionId?: string): Promise<any> {
           const response = await fetch(`${API_BASE_URL}/mcp/new-tool`, {
               method: 'POST',
               headers: { 'Content-Type': 'application/json' },
               body: JSON.stringify({ param, session_id: sessionId })
           });
           return response.json();
       }
   };
   ```

2. **Add to command handling** in `frontend/src/App.tsx`:
   ```typescript
   // In handleSubmit function
   if (tool === 'new_tool') {
       const result = await mcpApi.newTool(params.param, sessionId);
       // Handle result
   }
   ```

### 4. **Testing Strategy**

#### Backend Testing

```bash
# Run all backend tests
python -m pytest tests/backend/ -v

# Run specific test
python -m pytest tests/backend/test_command_router.py::test_specific_function -v

# Run with coverage
python -m pytest tests/backend/ --cov=backend --cov-report=html
```

#### Frontend Testing

```bash
# Run all frontend tests
cd frontend && npm test

# Run specific test
cd frontend && npm test -- --testNamePattern="NewTool"

# Run with coverage
cd frontend && npm test -- --coverage
```

#### Integration Testing

```bash
# Run integration tests
python -m pytest tests/integration/ -v

# Run specific workflow
python -m pytest tests/integration/test_workflows/test_new_workflow.py -v
```

## 🧪 Testing Guidelines

### Unit Tests

- **Backend**: Test individual functions and classes
- **Frontend**: Test individual components and utilities
- **Coverage**: Aim for >90% line coverage

### Integration Tests

- **API Tests**: Test complete API endpoints
- **Workflow Tests**: Test complete bioinformatics workflows
- **UI Tests**: Test user interactions and state management

### Test Data

- Use synthetic but realistic test data
- Keep test data minimal but representative
- Document data sources and formats

## 🔍 Code Quality

### Python Guidelines

- **PEP 8**: Follow Python style guide
- **Type Hints**: Use type annotations
- **Docstrings**: Document all functions and classes
- **Error Handling**: Proper exception handling

```python
def process_sequence(sequence: str) -> Dict[str, Any]:
    """
    Process a DNA sequence and return analysis results.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Dictionary containing analysis results
        
    Raises:
        ValueError: If sequence is invalid
    """
    if not sequence or not all(base in 'ATCG' for base in sequence.upper()):
        raise ValueError("Invalid DNA sequence")
    
    # Implementation
    return {"length": len(sequence), "gc_content": calculate_gc(sequence)}
```

### TypeScript Guidelines

- **Strict Mode**: Enable strict TypeScript settings
- **Interfaces**: Define proper interfaces for all data structures
- **Error Handling**: Proper error handling and user feedback
- **Component Props**: Use proper prop types

```typescript
interface SequenceData {
    name: string;
    sequence: string;
    length: number;
    gcContent: number;
}

interface SequenceVisualizerProps {
    data: SequenceData;
    onSelect?: (sequence: SequenceData) => void;
}

export const SequenceVisualizer: React.FC<SequenceVisualizerProps> = ({ 
    data, 
    onSelect 
}) => {
    // Implementation
};
```

## 🚀 Deployment

### Development Deployment

```bash
# Start the unified system
./start.sh

# Access application
# Frontend: http://localhost:5173
# Backend: http://localhost:8001
```

### Production Deployment

```bash
# Build frontend
cd frontend && npm run build

# Start production backend
python backend/main_with_mcp.py --production
```

## 🐛 Debugging

### Backend Debugging

```bash
# Enable debug logging
export MCP_LOG_LEVEL="DEBUG"

# Run with debug output
python -m pytest tests/ -v -s

# Check logs
tail -f backend/logs/app.log
```

### Frontend Debugging

```bash
# Enable React DevTools
# Install browser extension

# Check browser console
# Open Developer Tools (F12)

# Debug specific component
console.log('Debug data:', data);
```

### Common Issues

1. **Import Errors**:
   ```bash
   # Add project root to Python path
   export PYTHONPATH="${PYTHONPATH}:$(pwd)"
   ```

2. **CORS Issues**:
   ```python
   # In backend/main_with_mcp.py
   from fastapi.middleware.cors import CORSMiddleware
   
   app.add_middleware(
       CORSMiddleware,
       allow_origins=["http://localhost:5173"],
       allow_credentials=True,
       allow_methods=["*"],
       allow_headers=["*"],
   )
   ```

3. **Session Issues**:
   ```bash
   # Check session files
   ls -la backend/sessions/
   
   # Clear sessions
   rm -rf backend/sessions/*
   ```

## 📚 Documentation

### Code Documentation

- **Docstrings**: Document all functions and classes
- **Comments**: Explain complex logic
- **README**: Keep README files updated
- **API Docs**: Document all API endpoints

### User Documentation

- **User Guides**: Create guides for new features
- **Examples**: Provide usage examples
- **Troubleshooting**: Document common issues

## 🤝 Contributing

### Pull Request Process

1. **Create feature branch**
2. **Implement feature with tests**
3. **Update documentation**
4. **Run full test suite**
5. **Create pull request**
6. **Code review and merge**

### Code Review Checklist

- [ ] Code follows style guidelines
- [ ] Tests are included and passing
- [ ] Documentation is updated
- [ ] No breaking changes (or documented)
- [ ] Performance impact considered

## 📊 Performance

### Backend Performance

- **Async/Await**: Use async functions for I/O operations
- **Caching**: Cache expensive computations
- **Database**: Use connection pooling
- **Memory**: Monitor memory usage

### Frontend Performance

- **React Optimization**: Use React.memo, useMemo, useCallback
- **Bundle Size**: Monitor bundle size
- **Lazy Loading**: Implement code splitting
- **Caching**: Cache API responses

## 🔒 Security

### Backend Security

- **Input Validation**: Validate all inputs
- **Authentication**: Implement proper authentication
- **CORS**: Configure CORS properly
- **HTTPS**: Use HTTPS in production

### Frontend Security

- **XSS Prevention**: Sanitize user inputs
- **CSRF Protection**: Implement CSRF tokens
- **Content Security Policy**: Configure CSP headers

## 📈 Monitoring

### Backend Monitoring

```python
# Add logging
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info("Processing sequence: %s", sequence)
```

### Frontend Monitoring

```typescript
// Add error boundary
class ErrorBoundary extends React.Component {
    componentDidCatch(error: Error, errorInfo: React.ErrorInfo) {
        console.error('Error:', error, errorInfo);
        // Send to monitoring service
    }
}
```

## 🎯 Best Practices

### General

- **Version Control**: Use meaningful commit messages
- **Code Review**: Always review code before merging
- **Testing**: Write tests for all new features
- **Documentation**: Keep documentation updated
- **Performance**: Consider performance impact

### Bioinformatics Specific

- **Data Validation**: Validate biological data formats
- **Error Handling**: Handle biological data errors gracefully
- **User Feedback**: Provide clear feedback for bioinformatics operations
- **Data Privacy**: Handle sensitive biological data appropriately

This development guide should help you contribute effectively to the Helix.AI project. For specific questions, please refer to the documentation or create an issue. 