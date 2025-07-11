# Contributing to DataBloom.AI

Thank you for your interest in contributing to DataBloom.AI! This document provides guidelines for contributing to the project.

## 🚀 Getting Started

### Prerequisites

- Python 3.10+
- Node.js 16+
- Git
- Conda (recommended) or virtual environment

### Development Setup

1. **Fork and Clone**
   ```bash
   git clone https://github.com/yourusername/DataBloom.AI.git
   cd DataBloom.AI
   ```

2. **Setup Backend**
   ```bash
   cd backend
   conda create -n databloom-dev python=3.10
   conda activate databloom-dev
   pip install -r requirements.txt
   ```

3. **Setup Frontend**
   ```bash
   cd frontend
   npm install
   ```

## 🛠️ Development Guidelines

### Code Style

#### Python
- Follow PEP 8 style guidelines
- Use type hints for function parameters and return values
- Add docstrings for all functions and classes
- Keep functions focused and under 50 lines when possible

#### TypeScript/React
- Use TypeScript for all new code
- Follow React best practices
- Use functional components with hooks
- Add proper type definitions

### Adding New Features

#### 1. Natural Language Commands

When adding new natural language commands:

1. **Update Command Parser** (`tools/command_parser.py`):
   ```python
   def parse_command_raw(command: str, session_id: str = None) -> Dict[str, Any]:
       # Add new command patterns
       if "your new pattern" in command.lower():
           return {
               "action": "new_action",
               "tool": "new_tool",
               "parameters": {"param": "value"},
               "session_id": session_id
           }
   ```

2. **Create Tool Implementation** (`tools/new_tool.py`):
   ```python
   def run_new_tool_raw(param: str) -> Dict[str, Any]:
       """Execute the new tool with raw parameters."""
       try:
           # Implementation here
           return {
               "status": "success",
               "result": {"data": "processed_data"}
           }
       except Exception as e:
           return {
               "status": "error",
               "error": str(e)
           }
   ```

3. **Update Frontend** (`frontend/src/App.tsx`):
   ```typescript
   // Add to renderOutput function
   {actualResult.new_data && (
     <div className="bg-light p-3 border rounded mb-3">
       <h5>New Tool Results</h5>
       <pre>{JSON.stringify(actualResult.new_data, null, 2)}</pre>
     </div>
   )}
   ```

#### 2. API Endpoints

When adding new API endpoints:

1. **Add to FastAPI** (`backend/main_with_mcp.py`):
   ```python
   @app.post("/mcp/new-endpoint")
   async def new_endpoint(request: Dict[str, Any]):
       # Implementation
       return {"status": "success", "result": data}
   ```

2. **Add to Frontend API** (`frontend/src/services/mcpApi.ts`):
   ```typescript
   export const newEndpoint = async (params: any, sessionId?: string) => {
     return await apiCall('/mcp/new-endpoint', params, sessionId);
   };
   ```

### Testing

#### Backend Testing
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
```

#### Frontend Testing
```bash
cd frontend
npm test
```

#### API Testing
```bash
# Test new endpoint
curl -X POST http://localhost:8001/mcp/new-endpoint \
  -H "Content-Type: application/json" \
  -d '{"param": "value"}'
```

### Documentation

When adding new features:

1. **Update README.md** with new examples
2. **Add to appropriate docs/** files
3. **Update API documentation** comments
4. **Add inline comments** for complex logic

## 🐛 Bug Reports

When reporting bugs:

1. **Use the GitHub issue template**
2. **Include steps to reproduce**
3. **Provide error messages and logs**
4. **Specify your environment** (OS, Python/Node versions)

## 🔄 Pull Request Process

1. **Create a feature branch**:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make your changes** following the guidelines above

3. **Test thoroughly**:
   - Run backend tests
   - Run frontend tests
   - Test the complete workflow
   - Verify natural language commands work

4. **Commit with clear messages**:
   ```bash
   git commit -m "feat: add new natural language command for sequence analysis"
   git commit -m "fix: resolve session management issue in command parser"
   ```

5. **Push and create PR**:
   ```bash
   git push origin feature/your-feature-name
   ```

6. **Fill out the PR template** with:
   - Description of changes
   - Testing performed
   - Screenshots (if UI changes)
   - Related issues

## 📋 Code Review Guidelines

### For Reviewers

- **Check functionality**: Does the feature work as expected?
- **Review code quality**: Is the code readable and maintainable?
- **Verify testing**: Are there adequate tests?
- **Check documentation**: Is the feature documented?

### For Contributors

- **Respond to feedback** promptly
- **Make requested changes** clearly
- **Add tests** if requested
- **Update documentation** if needed

## 🎯 Areas for Contribution

### High Priority
- **New bioinformatics tools**: Add more analysis capabilities
- **Enhanced natural language parsing**: Improve command understanding
- **Better error handling**: More user-friendly error messages
- **Performance optimization**: Faster processing of large datasets

### Medium Priority
- **UI improvements**: Better visualizations and user experience
- **Documentation**: More examples and tutorials
- **Testing**: Additional test coverage
- **CI/CD**: Automated testing and deployment

### Low Priority
- **Code refactoring**: Clean up existing code
- **Minor bug fixes**: Small issues and edge cases
- **Documentation updates**: Typos and clarifications

## 📞 Getting Help

- **GitHub Issues**: For bugs and feature requests
- **GitHub Discussions**: For questions and general discussion
- **Documentation**: Check the docs/ directory
- **Code Examples**: See tests/sample_files/ for examples

## 📄 License

By contributing to DataBloom.AI, you agree that your contributions will be licensed under the MIT License.

Thank you for contributing to DataBloom.AI! 🧬 