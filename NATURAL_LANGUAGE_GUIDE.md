# Natural Language Command Guide

This guide provides comprehensive examples and best practices for using natural language commands in DataBloom.AI.

## 🚀 Quick Start

### 1. Enable Natural Language Mode
- Open the web interface at http://localhost:5173
- Click the "Natural Language" radio button
- You'll see the input placeholder change to reflect natural language mode

### 2. Basic Workflow
```
# Step 1: Generate some variants
mutate sequence ACTGTTGAC with 10 variants

# Step 2: Select variants using natural language
from the sequence variants, pick 5 sequences randomly and output them

# Step 3: Analyze the results
analyze the selected variants and show me the statistics
```

## 📝 Command Examples

### Variant Selection Commands

| Command | Action | Description |
|---------|--------|-------------|
| `"from the sequence variants, pick 10 sequences randomly and output them"` | Random Selection | Selects 10 random variants from previous mutation results |
| `"select 5 sequences with the highest mutation rate"` | Smart Selection | Uses mutation rate as selection criteria |
| `"pick the most diverse variants from the results"` | Diversity Selection | Selects variants based on sequence diversity |
| `"choose 3 sequences with the longest length"` | Length Selection | Selects variants based on sequence length |
| `"select variants that have the highest GC content"` | GC Content Selection | Selects variants based on GC content |

### Analysis Commands

| Command | Action | Description |
|---------|--------|-------------|
| `"analyze the alignment and show me the most conserved regions"` | Conservation Analysis | Identifies conserved regions in alignment |
| `"calculate statistics for the selected variants"` | Statistical Analysis | Performs comprehensive statistical analysis |
| `"show me the mutation patterns in the variants"` | Pattern Analysis | Analyzes mutation patterns and frequencies |
| `"compare the variants and find the most similar ones"` | Similarity Analysis | Compares variants for similarity |

### Multi-step Workflow Commands

| Command | Action | Description |
|---------|--------|-------------|
| `"mutate this sequence, then align the variants and pick the best ones"` | Complete Pipeline | Chains mutation, alignment, and selection |
| `"generate variants, analyze them, and select the top 5"` | Analysis Pipeline | Combines generation, analysis, and selection |
| `"create mutations, align them, and show me the most diverse set"` | Diversity Pipeline | Focuses on creating diverse variant sets |

### Data Science Commands

| Command | Action | Description |
|---------|--------|-------------|
| `"perform statistical analysis on the sequence data"` | Statistical Analysis | Comprehensive statistical analysis |
| `"create a correlation matrix for the variants"` | Correlation Analysis | Analyzes relationships between variants |
| `"generate feature engineering suggestions"` | Feature Engineering | Suggests new features for analysis |
| `"create visualizations for the analysis results"` | Visualization | Generates plots and charts |

## 🔧 Advanced Usage

### Session Management

Natural language commands maintain context across multiple operations:

```
# Session 1: Generate and analyze variants
mutate sequence ACTGTTGAC with 20 variants
from the sequence variants, pick 10 sequences randomly
analyze the selected variants for diversity

# Session 2: Work with different data
mutate sequence GCTAGCTA with 15 variants
select 5 variants with highest mutation rate
compare these variants with the previous results
```

### Complex Commands

You can combine multiple operations in a single command:

```
# Complex analysis with multiple criteria
mutate sequence ACTGTTGAC, then align the variants, 
select the 5 most diverse ones, and analyze their properties

# Multi-step analysis
generate 20 variants, analyze their mutation patterns, 
select the top 10 based on diversity, and create visualizations
```

### File Upload Integration

Combine file uploads with natural language commands:

1. **Upload a FASTA file** by dragging it to the input area
2. **Choose an action** from the modal (align, mutate, analyze)
3. **Use natural language** to further process the results:
   ```
   # After uploading and aligning sequences
   from the alignment results, pick the most conserved sequences
   
   # After uploading and mutating sequences
   select 5 variants with the highest quality scores
   ```

## 🎯 Best Practices

### 1. Clear and Specific Commands
```
✅ Good: "from the sequence variants, pick 10 sequences randomly and output them"
❌ Vague: "pick some variants"

✅ Good: "select 5 sequences with the highest mutation rate"
❌ Unclear: "get the best ones"
```

### 2. Use Descriptive Language
```
✅ Descriptive: "analyze the alignment and show me the most conserved regions"
❌ Generic: "analyze this"

✅ Specific: "select variants that have the highest GC content"
❌ Vague: "pick good ones"
```

### 3. Chain Operations Logically
```
✅ Logical flow: "mutate sequence, align variants, select best ones"
❌ Confusing: "select variants, then mutate sequence"
```

### 4. Provide Context
```
✅ With context: "from the previous mutation results, select 5 variants"
❌ Without context: "select 5 variants"
```

## 🐛 Troubleshooting

### Common Issues

#### 1. "No previous mutation results found in session"
**Cause**: Trying to select variants without first generating them
**Solution**: 
```
# First, generate variants
mutate sequence ACTGTTGAC with 10 variants

# Then select from them
from the sequence variants, pick 5 sequences randomly
```

#### 2. "Session not found" error
**Cause**: Session expired or invalid session ID
**Solution**: 
- Refresh the page to create a new session
- Or restart the backend server

#### 3. Command not recognized
**Cause**: Command format not supported
**Solution**: 
- Use more specific language
- Check the example commands above
- Try structured command mode for complex operations

#### 4. No output displayed
**Cause**: Frontend rendering issue or API error
**Solution**: 
- Check browser console for errors
- Verify backend server is running
- Check network tab for API responses

### Debugging Commands

#### Test Command Parsing
```bash
# Test command parser directly
python -c "
import sys; sys.path.append('../tools')
from command_parser import parse_command_raw
result = parse_command_raw('your command here')
print(result)
"
```

#### Test Complete Workflow
```bash
# Test natural language command handling
curl -X POST http://localhost:8001/mcp/handle-natural-command \
  -H "Content-Type: application/json" \
  -d '{"command": "your command here", "session_id": "test"}'
```

## 📊 Expected Outputs

### Variant Selection Output
```json
{
  "status": "success",
  "selected_variants": [
    {"name": "variant_1", "sequence": "ACTGTTGAC"},
    {"name": "variant_2", "sequence": "ACTGATGAC"}
  ],
  "analysis": {
    "total_variants": 10,
    "selected_variants": 5,
    "selection_ratio": 0.5,
    "criteria_used": "random"
  }
}
```

### Analysis Output
```json
{
  "status": "success",
  "analysis": {
    "statistics": {
      "mean_length": 9.2,
      "gc_content": 0.45,
      "mutation_rate": 0.12
    },
    "visualization": {
      "plot_data": {...},
      "plot_layout": {...}
    }
  }
}
```

## 🔄 Switching Between Modes

### Natural Language Mode
- **Use for**: Complex workflows, multi-step operations, context-aware commands
- **Best for**: "from the sequence variants, pick 10 sequences randomly"

### Structured Command Mode
- **Use for**: Simple operations, precise control, direct tool access
- **Best for**: "align sequences ACTGTTGAC ACTGCATCC"

## 📈 Performance Tips

### For Large Datasets
1. **Use smaller initial sets**: Start with 10-20 variants
2. **Chain operations efficiently**: Combine related operations
3. **Clear session periodically**: Restart for very long workflows

### For Complex Workflows
1. **Break into steps**: Use multiple commands instead of one complex command
2. **Use session context**: Leverage the session to maintain state
3. **Check results**: Verify each step before proceeding

## 🎓 Learning Path

### Beginner
1. Start with simple commands: "mutate sequence ACTGTTGAC"
2. Learn selection: "from the variants, pick 5 randomly"
3. Try analysis: "analyze the selected variants"

### Intermediate
1. Combine operations: "mutate, align, and select best variants"
2. Use file uploads with natural language
3. Explore different selection criteria

### Advanced
1. Create complex workflows with multiple steps
2. Use session context for long-running analyses
3. Combine with structured commands for precise control

## 📚 Additional Resources

- **API Documentation**: http://localhost:8001/docs
- **Example Files**: `tests/sample_files/`
- **Backend Logs**: Check console for detailed error messages
- **Session Files**: Located in `backend/sessions/` (not committed to git)

---

*For more help, check the main README.md or create an issue on GitHub.* 