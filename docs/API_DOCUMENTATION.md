# 📚 Helix.AI API Documentation

This document provides comprehensive documentation for all Helix.AI API endpoints.

## 🏠 Base URL

- **Development**: `http://localhost:8001`
- **Production**: `https://your-domain.com`

## 🔐 Authentication

Currently, the API uses session-based authentication. Each request should include a `session_id` in the request body or query parameters.

## 📋 API Endpoints

### Health Check

#### GET `/health`

Check the health status of the API server.

**Response**:
```json
{
  "status": "healthy",
  "timestamp": "2024-01-XX",
  "version": "0.2.0"
}
```

### Session Management

#### POST `/create_session`

Create a new session for tracking workflow state.

**Request**:
```json
{}
```

**Response**:
```json
{
  "session_id": "uuid-string",
  "created_at": "2024-01-XX",
  "status": "created"
}
```

#### GET `/session/{session_id}`

Get session information and history.

**Response**:
```json
{
  "session_id": "uuid-string",
  "created_at": "2024-01-XX",
  "updated_at": "2024-01-XX",
  "history": [
    {
      "command": "visualize the phylogenetic tree",
      "timestamp": "2024-01-XX",
      "result": {...}
    }
  ]
}
```

### Natural Language Commands

#### POST `/execute`

Execute a natural language command.

**Request**:
```json
{
  "command": "visualize the phylogenetic tree",
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "tool": "phylogenetic_tree",
    "parameters": {...},
    "output": {...}
  },
  "session_id": "uuid-string"
}
```

### MCP (Model Context Protocol) Endpoints

#### GET `/mcp/tools`

Get list of available MCP tools.

**Response**:
```json
{
  "tools": [
    {
      "name": "sequence_alignment",
      "description": "Align DNA sequences",
      "parameters": {...}
    }
  ]
}
```

#### POST `/mcp/sequence-alignment`

Perform sequence alignment.

**Request**:
```json
{
  "sequences": ">seq1 ATGCGATCG\n>seq2 ATGCGATC",
  "method": "clustalw",
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "alignment": [...],
    "statistics": {
      "num_sequences": 2,
      "alignment_length": 9,
      "average_identity": 88.89
    }
  },
  "session_id": "uuid-string"
}
```

#### POST `/mcp/phylogenetic-tree`

Generate phylogenetic tree from aligned sequences.

**Request**:
```json
{
  "aligned_sequences": ">seq1 ATGCGATCG\n>seq2 ATGCGATC",
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "tree_newick": "(seq1:0.1,seq2:0.1);",
    "ete_visualization": {
      "svg": "<svg>...</svg>"
    },
    "statistics": {
      "num_sequences": 2,
      "tree_length": 0.2
    }
  },
  "session_id": "uuid-string"
}
```

#### POST `/mcp/clustering-analysis`

Perform clustering analysis on sequences.

**Request**:
```json
{
  "aligned_sequences": ">seq1 ATGCGATCG\n>seq2 ATGCGATC",
  "num_clusters": 5,
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "clusters": [
      {
        "cluster_id": 0,
        "size": 2,
        "representative": "seq1",
        "average_distance": 0.1234,
        "sequences": ["seq1", "seq2"]
      }
    ],
    "representatives": ["seq1"],
    "total_sequences": 2,
    "num_clusters": 1
  },
  "session_id": "uuid-string"
}
```

#### POST `/mcp/plasmid-visualization`

Generate plasmid visualization.

**Request**:
```json
{
  "vector_name": "pUC19",
  "cloning_sites": "EcoRI, BamHI, HindIII",
  "insert_sequence": "ATGCGATCG",
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "name": "pUC19",
    "sequence": "ATGCGATCG...",
    "features": [
      {
        "name": "EcoRI",
        "start": 100,
        "end": 106,
        "type": "restriction_site"
      }
    ],
    "description": "pUC19 plasmid with insert"
  },
  "session_id": "uuid-string"
}
```

#### POST `/mcp/plasmid-for-representatives`

Generate plasmid visualization for representative sequences.

**Request**:
```json
{
  "vector_name": "pUC19",
  "cloning_sites": "EcoRI, BamHI, HindIII",
  "representative_sequences": ["seq1", "seq2"],
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "representatives": [
      {
        "name": "seq1",
        "sequence": "ATGCGATCG",
        "plasmid_data": {...}
      }
    ]
  },
  "session_id": "uuid-string"
}
```

#### POST `/mcp/mutate-sequence`

Generate sequence mutations.

**Request**:
```json
{
  "sequence": "ATGCGATCG",
  "num_variants": 10,
  "mutation_rate": 0.1,
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "original_sequence": "ATGCGATCG",
    "variants": [
      {
        "sequence": "ATGCGATCA",
        "mutations": [{"position": 8, "from": "G", "to": "A"}]
      }
    ],
    "statistics": {
      "num_variants": 10,
      "mutation_rate": 0.1
    }
  },
  "session_id": "uuid-string"
}
```

#### POST `/mcp/select-variants`

Select variants based on criteria.

**Request**:
```json
{
  "variants": [...],
  "criteria": "diversity",
  "num_selected": 5,
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "selected_variants": [...],
    "selection_criteria": "diversity",
    "num_selected": 5
  },
  "session_id": "uuid-string"
}
```

#### POST `/mcp/analyze-sequence-data`

Perform statistical analysis on sequence data.

**Request**:
```json
{
  "sequences": [...],
  "analysis_type": "statistics",
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "statistics": {
      "length_stats": {...},
      "gc_content_stats": {...},
      "conservation_stats": {...}
    },
    "visualizations": {...}
  },
  "session_id": "uuid-string"
}
```

#### POST `/mcp/dna-vendor-research`

Research DNA synthesis vendors.

**Request**:
```json
{
  "sequence_length": 1000,
  "quantity": 10,
  "requirements": ["high_quality", "fast_turnaround"],
  "session_id": "uuid-string"
}
```

**Response**:
```json
{
  "success": true,
  "result": {
    "vendors": [
      {
        "name": "Vendor A",
        "pricing": "$0.25/bp",
        "turnaround": "5-7 days",
        "quality": "high",
        "recommendation_score": 0.95
      }
    ],
    "recommendations": [...],
    "pricing_comparison": {...}
  },
  "session_id": "uuid-string"
}
```

## 🔄 Error Handling

### Error Response Format

```json
{
  "success": false,
  "error": "Error message",
  "error_code": "ERROR_CODE",
  "session_id": "uuid-string"
}
```

### Common Error Codes

- `INVALID_SESSION`: Session not found or expired
- `INVALID_COMMAND`: Command could not be parsed
- `TOOL_NOT_FOUND`: Requested tool not available
- `INVALID_PARAMETERS`: Invalid parameters provided
- `PROCESSING_ERROR`: Error during tool execution
- `FILE_UPLOAD_ERROR`: Error processing uploaded file

### HTTP Status Codes

- `200 OK`: Request successful
- `400 Bad Request`: Invalid request parameters
- `404 Not Found`: Endpoint or resource not found
- `422 Unprocessable Entity`: Validation error
- `500 Internal Server Error`: Server error

## 📊 Rate Limiting

Currently, no rate limiting is implemented. However, for production deployments, consider implementing:

- **Per-session limits**: Limit requests per session
- **Per-endpoint limits**: Different limits for different endpoints
- **Global limits**: Overall API usage limits

## 🔒 Security Considerations

### Input Validation

All endpoints validate input parameters:

- **Sequence validation**: Ensure valid DNA sequences
- **File validation**: Validate uploaded file formats
- **Parameter validation**: Validate all input parameters

### CORS Configuration

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:5175"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```

## 📈 Performance

### Response Times

- **Health check**: < 100ms
- **Simple commands**: < 1s
- **Complex analysis**: < 30s
- **File processing**: < 60s

### Optimization Tips

- **Use sessions**: Maintain session state for better performance
- **Batch operations**: Combine multiple operations when possible
- **Cache results**: Results are cached within sessions
- **Async processing**: Long-running operations are async

## 🧪 Testing

### Test Endpoints

```bash
# Health check
curl http://localhost:8001/health

# Create session
curl -X POST http://localhost:8001/create_session

# Execute command
curl -X POST http://localhost:8001/execute \
  -H "Content-Type: application/json" \
  -d '{"command": "visualize the phylogenetic tree", "session_id": "test"}'
```

### Example Workflows

#### Complete Phylogenetic Analysis

```bash
# 1. Create session
SESSION_ID=$(curl -s -X POST http://localhost:8001/create_session | jq -r '.session_id')

# 2. Upload sequences
curl -X POST http://localhost:8001/mcp/phylogenetic-tree \
  -H "Content-Type: application/json" \
  -d "{\"aligned_sequences\": \">seq1 ATGCGATCG\n>seq2 ATGCGATC\", \"session_id\": \"$SESSION_ID\"}"

# 3. Select representatives
curl -X POST http://localhost:8001/mcp/clustering-analysis \
  -H "Content-Type: application/json" \
  -d "{\"num_clusters\": 5, \"session_id\": \"$SESSION_ID\"}"

# 4. Insert into plasmid
curl -X POST http://localhost:8001/mcp/plasmid-visualization \
  -H "Content-Type: application/json" \
  -d "{\"vector_name\": \"pUC19\", \"session_id\": \"$SESSION_ID\"}"
```

## 📚 SDK Examples

### Python SDK

```python
import requests

class HelixAI:
    def __init__(self, base_url="http://localhost:8001"):
        self.base_url = base_url
        self.session_id = None
    
    def create_session(self):
        response = requests.post(f"{self.base_url}/create_session")
        self.session_id = response.json()["session_id"]
        return self.session_id
    
    def execute_command(self, command):
        response = requests.post(f"{self.base_url}/execute", json={
            "command": command,
            "session_id": self.session_id
        })
        return response.json()
    
    def visualize_tree(self, sequences):
        response = requests.post(f"{self.base_url}/mcp/phylogenetic-tree", json={
            "aligned_sequences": sequences,
            "session_id": self.session_id
        })
        return response.json()

# Usage
client = HelixAI()
client.create_session()
result = client.execute_command("visualize the phylogenetic tree")
```

### JavaScript SDK

```javascript
class HelixAI {
    constructor(baseUrl = 'http://localhost:8001') {
        this.baseUrl = baseUrl;
        this.sessionId = null;
    }
    
    async createSession() {
        const response = await fetch(`${this.baseUrl}/create_session`, {
            method: 'POST'
        });
        const data = await response.json();
        this.sessionId = data.session_id;
        return this.sessionId;
    }
    
    async executeCommand(command) {
        const response = await fetch(`${this.baseUrl}/execute`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                command,
                session_id: this.sessionId
            })
        });
        return response.json();
    }
}

// Usage
const client = new HelixAI();
await client.createSession();
const result = await client.executeCommand('visualize the phylogenetic tree');
```

## 🔄 WebSocket Support

For real-time updates, consider implementing WebSocket endpoints:

```python
@app.websocket("/ws/{session_id}")
async def websocket_endpoint(websocket: WebSocket, session_id: str):
    await websocket.accept()
    try:
        while True:
            data = await websocket.receive_text()
            # Process command and send updates
            await websocket.send_text(json.dumps(result))
    except WebSocketDisconnect:
        pass
```

## 📝 Changelog

### Version 0.2.0
- Added clustering analysis endpoints
- Enhanced plasmid visualization
- Improved error handling
- Added comprehensive documentation

### Version 0.1.0
- Initial API release
- Basic bioinformatics tools
- Session management
- Natural language commands

For more information, see the [main documentation](README.md) and [development guide](DEVELOPMENT_GUIDE.md). 