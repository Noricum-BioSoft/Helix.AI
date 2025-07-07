# History Tracking System

The DataBloom.AI backend now includes a comprehensive history tracking system that maintains context across user interactions. This allows the system to remember previous results and use them in subsequent operations.

## Overview

The history tracking system enables the user story you described:

1. **Step 1**: User uploads a FASTA file and asks to generate mutations of the given sequences. The result is a large set of sequence variants.
2. **Step 2**: The user instructs to pick N variants according to some criteria (e.g., diversity). The backend knows about the sequence variants from Step 1.

## Components

### 1. History Manager (`history_manager.py`)

The core component that manages user sessions and their history:

- **Session Management**: Creates and manages user sessions
- **History Storage**: Persists operation history to disk
- **Result Retrieval**: Retrieves previous results for reuse
- **Session Summary**: Provides overview of session activity

### 2. Variant Selection Tool (`tools/variant_selection.py`)

A new tool that can select variants from previous mutation results:

- **Diversity-based selection**: Selects diverse variants using sequence similarity
- **Random selection**: Randomly selects variants
- **Length-based selection**: Selects variants based on sequence length
- **Custom filtering**: Allows custom filtering criteria

### 3. Enhanced API Endpoints

The FastAPI server now includes session-aware endpoints:

- `/session/create`: Create a new session
- `/session/{session_id}`: Get session information
- `/mcp/select-variants`: Select variants from previous results
- All existing endpoints now support session tracking

## Usage Examples

### Creating a Session

```python
from history_manager import history_manager

# Create a new session
session_id = history_manager.create_session("user123")
print(f"Session created: {session_id}")
```

### Tracking Operations

```python
# Track a mutation operation
mutation_result = {
    "status": "success",
    "variants": [...],
    "statistics": {...}
}

history_manager.add_history_entry(
    session_id,
    "Generate 50 variants with mutation rate 0.15",
    "mutate_sequence",
    mutation_result,
    {"num_variants": 50, "mutation_rate": 0.15}
)
```

### Retrieving Previous Results

```python
# Get the latest mutation results
previous_mutation = history_manager.get_latest_result(session_id, "mutate_sequence")

if previous_mutation:
    variants = previous_mutation.get("variants", [])
    print(f"Found {len(variants)} previous variants")
```

### Selecting Variants

```python
from tools.variant_selection import run_variant_selection_raw

# Select diverse variants from previous results
result = run_variant_selection_raw(
    session_id=session_id,
    selection_criteria="diversity",
    num_variants=10
)

selected_variants = result.get("selected_variants", [])
print(f"Selected {len(selected_variants)} diverse variants")
```

## API Endpoints

### Session Management

#### Create Session
```http
POST /session/create
Content-Type: application/json

{
    "user_id": "optional_user_id"
}
```

Response:
```json
{
    "success": true,
    "session_id": "uuid-string",
    "message": "Session created successfully"
}
```

#### Get Session Info
```http
GET /session/{session_id}
```

Response:
```json
{
    "success": true,
    "session": {
        "session_id": "uuid-string",
        "user_id": "user123",
        "created_at": "2024-01-01T12:00:00",
        "history": [...],
        "results": {...}
    },
    "summary": {
        "total_operations": 2,
        "tool_usage": {
            "mutate_sequence": 1,
            "select_variants": 1
        }
    }
}
```

### Variant Selection

#### Select Variants
```http
POST /mcp/select-variants
Content-Type: application/json

{
    "session_id": "uuid-string",
    "selection_criteria": "diversity",
    "num_variants": 10,
    "custom_filters": {
        "min_length": 8,
        "max_length": 12
    }
}
```

Response:
```json
{
    "success": true,
    "result": {
        "status": "success",
        "session_id": "uuid-string",
        "selection_criteria": "diversity",
        "num_variants_selected": 10,
        "selected_variants": [...],
        "analysis": {...},
        "plot": {...}
    },
    "session_id": "uuid-string"
}
```

## Selection Criteria

### Diversity
Selects variants that are most different from each other using sequence similarity metrics.

### Random
Randomly selects variants from the available pool.

### Length
Selects variants based on sequence length criteria.

### Custom
Uses custom filtering criteria:
- `min_length`: Minimum sequence length
- `max_length`: Maximum sequence length
- `gc_content_range`: [min_gc, max_gc] for GC content
- `substring`: Must contain specific substring

## User Story Implementation

The system now supports the complete user story:

1. **Step 1**: User uploads FASTA file and requests mutations
   ```python
   # This creates a session and tracks the mutation operation
   response = await mutate_sequence_mcp(MutationRequest(
       sequence="ACTGTTGAC",
       num_variants=50,
       mutation_rate=0.15,
       session_id=session_id
   ))
   ```

2. **Step 2**: User selects variants from previous results
   ```python
   # This retrieves previous mutation results and selects variants
   response = await select_variants_mcp(VariantSelectionRequest(
       session_id=session_id,
       selection_criteria="diversity",
       num_variants=10
   ))
   ```

## Testing

Run the test script to verify the history tracking system:

```bash
cd backend
python test_history.py
```

This will demonstrate:
- Session creation and management
- Operation tracking
- Result retrieval
- Variant selection
- Session summaries

## File Structure

```
backend/
├── history_manager.py          # Core history management
├── main_with_mcp.py           # Enhanced FastAPI server
├── test_history.py            # Test script
└── sessions/                  # Session storage directory
    └── *.json                 # Session files

tools/
└── variant_selection.py       # Variant selection tool
```

## Benefits

1. **Context Preservation**: Maintains context across multiple operations
2. **Result Reuse**: Allows subsequent operations to use previous results
3. **Workflow Support**: Enables complex multi-step workflows
4. **Session Management**: Organizes operations by user sessions
5. **Analytics**: Provides insights into tool usage and operation patterns

The history tracking system transforms DataBloom.AI from a simple tool executor into a comprehensive bioinformatics workflow platform that can handle complex, multi-step analyses while maintaining full context and history. 