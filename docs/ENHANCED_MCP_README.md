# Enhanced Bioinformatics MCP Server

This enhanced version of the Bioinformatics MCP Server provides improved validation, error handling, and additional bioinformatics tools.

## 🚀 New Features

### **Enhanced Validation**
- **Sequence Validation**: Validates DNA/RNA sequences for proper characters
- **FASTA Format Validation**: Ensures proper FASTA format structure
- **Parameter Validation**: Validates all input parameters with proper ranges
- **Custom ValidationError**: Specific error handling for validation failures

### **New Bioinformatics Tools**

#### **1. Sequence Validation Tool**
```python
# Validate individual sequences
result = await handle_sequence_validation({
    "sequence": "ACTGTTGAC"
})

# Validate FASTA content
result = await handle_sequence_validation({
    "fasta_content": ">seq1\nACTGTTGAC\n>seq2\nACTGCATCC"
})
```

#### **2. Sequence Statistics Tool**
```python
# Calculate comprehensive sequence statistics
result = await handle_sequence_statistics({
    "sequence": "ACTGTTGAC",
    "include_composition": True
})
```

**Features:**
- Sequence length calculation
- GC content analysis
- Nucleotide composition breakdown
- Gap content analysis

#### **3. Reverse Complement Tool**
```python
# Generate reverse complement of DNA sequence
result = await handle_reverse_complement({
    "sequence": "ACTGTTGAC"
})
```

### **Enhanced Error Handling**
- **Structured Error Responses**: Consistent error format with timestamps
- **Validation Errors**: Specific error messages for invalid inputs
- **Execution Errors**: Detailed error information for tool failures
- **Performance Tracking**: Execution time measurement for all operations

### **Performance Improvements**
- **Execution Timing**: All operations include execution time tracking
- **Input Validation**: Pre-validation prevents unnecessary processing
- **Optimized Parsing**: Improved FASTA parsing with error handling
- **Memory Efficiency**: Better memory management for large sequences

## 📁 Files

- `mcp_server_enhanced.py` - Enhanced MCP server implementation
- `mcp_config_enhanced.json` - Enhanced configuration file
- `test_enhanced_mcp.py` - Comprehensive test suite
- `start_enhanced_mcp.sh` - Enhanced startup script

## 🛠️ Installation

1. **Install Dependencies**:
```bash
pip install -r requirements.txt
```

2. **Make Startup Script Executable**:
```bash
chmod +x start_enhanced_mcp.sh
```

## 🚀 Usage

### **Starting the Enhanced Server**

```bash
# Using the startup script
./start_enhanced_mcp.sh

# Or directly
python3 mcp_server_enhanced.py
```

### **Configuration**

The enhanced server uses `mcp_config_enhanced.json`:

```json
{
  "mcpServers": {
    "enhanced-bioinformatics": {
      "command": "python",
      "args": ["mcp_server_enhanced.py"],
      "env": {
        "PYTHONPATH": ".",
        "MCP_LOG_LEVEL": "INFO",
        "MCP_VALIDATION_STRICT": "true"
      }
    }
  },
  "settings": {
    "validation": {
      "strict_mode": true,
      "max_sequence_length": 100000,
      "max_variants": 1000
    },
    "logging": {
      "level": "INFO",
      "file": "mcp_server.log"
    }
  }
}
```

## 🧪 Testing

### **Run Comprehensive Tests**

```bash
python3 test_enhanced_mcp.py
```

### **Test Categories**

1. **Import and Initialization**: Tests server module imports
2. **Sequence Validation**: Tests DNA/RNA sequence validation
3. **FASTA Validation**: Tests FASTA format validation
4. **FASTA Parsing**: Tests FASTA content parsing
5. **Enhanced Tools**: Tests new bioinformatics tools
6. **Error Handling**: Tests validation and error scenarios
7. **Performance**: Tests execution timing and performance
8. **Integration**: Tests integration with existing tools

## 🔧 API Endpoints

### **Enhanced Tools**

#### **Sequence Validation**
```http
POST /mcp/sequence-validation
Content-Type: application/json

{
  "sequence": "ACTGTTGAC",
  "fasta_content": ">seq1\nACTGTTGAC"
}
```

#### **Sequence Statistics**
```http
POST /mcp/sequence-statistics
Content-Type: application/json

{
  "sequence": "ACTGTTGAC",
  "include_composition": true
}
```

#### **Reverse Complement**
```http
POST /mcp/reverse-complement
Content-Type: application/json

{
  "sequence": "ACTGTTGAC"
}
```

## 📊 Response Format

### **Success Response**
```json
{
  "status": "success",
  "tool": "sequence_statistics",
  "sequence_length": 9,
  "statistics": {
    "length": 9,
    "gc_content": 0.44,
    "at_content": 0.56,
    "composition": {
      "A": 2,
      "T": 3,
      "C": 2,
      "G": 2
    }
  },
  "execution_time": 0.0023,
  "timestamp": "2024-01-15T10:30:45.123456"
}
```

### **Error Response**
```json
{
  "error": "validation_error",
  "message": "Invalid DNA/RNA sequence: contains invalid characters",
  "tool": "sequence_statistics",
  "timestamp": "2024-01-15T10:30:45.123456"
}
```

## 🔍 Validation Rules

### **Sequence Validation**
- **Valid Characters**: A, T, C, G, N, -
- **Case Insensitive**: Converts to uppercase
- **Whitespace Handling**: Removes spaces and newlines
- **Empty Sequences**: Rejected

### **FASTA Validation**
- **Header Format**: Must start with '>'
- **Sequence Content**: Must contain valid DNA/RNA characters
- **Structure**: Must have at least one sequence
- **Format**: Proper FASTA structure

### **Parameter Validation**
- **num_variants**: 1-1000 range
- **mutation_rate**: 0.0-1.0 range
- **algorithm**: Must be one of ["clustal", "muscle", "mafft"]
- **output_format**: Must be one of ["png", "svg", "pdf"]

## 📈 Performance Features

### **Execution Timing**
- All operations include execution time measurement
- Timing information included in response
- Performance monitoring for optimization

### **Memory Management**
- Efficient sequence parsing
- Optimized FASTA handling
- Memory-conscious data structures

### **Logging**
- Structured logging with timestamps
- Error tracking and debugging
- Performance metrics logging

## 🔧 Development

### **Adding New Tools**

1. **Define Tool Schema**:
```python
Tool(
    name="new_tool",
    description="Description of the new tool",
    inputSchema={
        "type": "object",
        "properties": {
            "parameter": {
                "type": "string",
                "description": "Parameter description"
            }
        },
        "required": ["parameter"]
    }
)
```

2. **Add Handler Function**:
```python
async def handle_new_tool(arguments: Dict[str, Any]) -> Dict[str, Any]:
    # Validate input
    parameter = arguments.get("parameter", "")
    if not parameter:
        raise ValidationError("Parameter is required")
    
    # Process request
    result = process_parameter(parameter)
    
    return {
        "status": "success",
        "tool": "new_tool",
        "result": result
    }
```

3. **Update Tool Router**:
```python
elif name == "new_tool":
    result = await handle_new_tool(arguments)
```

### **Testing New Features**

1. **Add Unit Tests**:
```python
# Test validation
assert validate_sequence("ACTGTTGAC") == True
assert validate_sequence("ACTGTTGAC123") == False
```

2. **Add Integration Tests**:
```python
# Test tool functionality
result = await handle_new_tool({"parameter": "test"})
assert result["status"] == "success"
```

## 🐛 Troubleshooting

### **Common Issues**

1. **Validation Errors**:
   - Check input format and characters
   - Verify parameter ranges
   - Ensure required parameters are provided

2. **Import Errors**:
   - Verify PYTHONPATH is set correctly
   - Check all dependencies are installed
   - Ensure file paths are correct

3. **Performance Issues**:
   - Check execution time in response
   - Monitor memory usage for large sequences
   - Review logs for bottlenecks

### **Debugging**

1. **Enable Debug Logging**:
```bash
export MCP_LOG_LEVEL="DEBUG"
```

2. **Check Logs**:
```bash
tail -f logs/mcp_server.log
```

3. **Run Tests**:
```bash
python3 test_enhanced_mcp.py
```

## 📝 Changelog

### **Version 1.1.0 (Enhanced)**
- ✅ Added sequence validation tool
- ✅ Added sequence statistics tool
- ✅ Added reverse complement tool
- ✅ Enhanced error handling with structured responses
- ✅ Added execution time tracking
- ✅ Improved input validation
- ✅ Added comprehensive test suite
- ✅ Enhanced logging and debugging

### **Version 1.0.0 (Original)**
- ✅ Basic MCP server implementation
- ✅ Sequence alignment tool
- ✅ Mutation analysis tool
- ✅ Data analysis tool
- ✅ Visualization tool

## 🤝 Contributing

When contributing to the enhanced MCP server:

1. **Follow Validation Patterns**: Use the existing validation functions
2. **Add Error Handling**: Include proper error handling and validation
3. **Write Tests**: Add comprehensive tests for new features
4. **Update Documentation**: Keep documentation current
5. **Performance**: Consider performance implications of changes

## 📄 License

This enhanced MCP server implementation is part of the Helix.AI project and follows the same license terms. 