# Helix.AI Microservices Architecture

## Overview

This is a refactored version of Helix.AI using a modern microservices architecture. The system is designed for scalability, maintainability, and high performance.

## Architecture Components

### 🏗️ Infrastructure Services

- **Redis**: Event bus, caching, and session storage
- **PostgreSQL**: Persistent data storage for sessions, workflows, and steps
- **MinIO**: Object storage for files, sequences, and visualizations

### 🔧 Core Services

#### API Gateway (`api-gateway`)
- **Port**: 8000
- **Purpose**: Single entry point for all client requests
- **Features**: 
  - Request routing and load balancing
  - Circuit breaker pattern
  - Request/response transformation
  - Authentication and authorization
  - Rate limiting

#### Workflow Engine (`workflow-engine`)
- **Port**: 8001
- **Purpose**: Orchestrates bioinformatics workflows
- **Features**:
  - Workflow definition and execution
  - Step dependency management
  - Event-driven workflow progression
  - State persistence
  - Error handling and recovery

### 🧬 Bioinformatics Services

#### Alignment Service (`alignment-service`)
- **Port**: 8002
- **Purpose**: Sequence alignment operations
- **Features**:
  - Multiple sequence alignment
  - Pairwise alignment
  - Consensus sequence generation
  - Alignment scoring

#### Mutation Service (`mutation-service`)
- **Port**: 8003
- **Purpose**: Sequence mutation and variant generation
- **Features**:
  - Random mutation generation
  - Site-directed mutagenesis
  - Mutation effect prediction
  - Variant library creation

#### Selection Service (`selection-service`)
- **Port**: 8004
- **Purpose**: Sequence selection and filtering
- **Features**:
  - Criteria-based selection
  - Random selection
  - Quality-based filtering
  - Diversity analysis

#### Plasmid Service (`plasmid-service`)
- **Port**: 8005
- **Purpose**: Plasmid design and visualization
- **Features**:
  - Vector design
  - Cloning site analysis
  - Restriction mapping
  - Plasmid visualization

#### Synthesis Service (`synthesis-service`)
- **Port**: 8006
- **Purpose**: DNA synthesis submission
- **Features**:
  - Vendor integration
  - Quote generation
  - Order management
  - Delivery tracking

### 🤖 Supporting Services

#### NLP Service (`nlp-service`)
- **Port**: 8007
- **Purpose**: Natural language processing
- **Features**:
  - Command parsing
  - Intent recognition
  - Entity extraction
  - Context understanding

#### Notification Service (`notification-service`)
- **Port**: 8008
- **Purpose**: Real-time notifications
- **Features**:
  - WebSocket connections
  - Event broadcasting
  - Status updates
  - Progress notifications

## Quick Start

### Prerequisites

- Docker and Docker Compose
- Node.js (for frontend development)
- Python 3.11+ (for local development)

### 1. Start the System

```bash
# Start all services
./start-microservices.sh

# Or manually with Docker Compose
docker-compose up -d
```

### 2. Verify Services

```bash
# Check service health
curl http://localhost:8000/health  # API Gateway
curl http://localhost:8001/health  # Workflow Engine
curl http://localhost:8002/health  # Alignment Service
```

### 3. Access the Application

- **Frontend**: http://localhost:5173
- **API Gateway**: http://localhost:8000
- **MinIO Console**: http://localhost:9001

## Development

### Local Development

```bash
# Start only infrastructure
docker-compose up -d redis postgres minio

# Run services locally
cd api-gateway && python main.py
cd workflow-engine && python main.py
cd bioinformatics-services/alignment-service && python main.py
```

### Adding New Services

1. Create service directory:
```bash
mkdir -p bioinformatics-services/new-service
```

2. Copy template files:
```bash
cp requirements.txt bioinformatics-services/new-service/
cp shared/models.py bioinformatics-services/new-service/shared/
cp shared/utils.py bioinformatics-services/new-service/shared/
```

3. Create service implementation:
```python
# bioinformatics-services/new-service/main.py
from fastapi import FastAPI
from shared.models import ServiceRequest, ServiceResponse
from shared.utils import setup_logging

app = FastAPI()
logger = setup_logging("new-service")

@app.post("/execute")
async def execute_service(request: ServiceRequest):
    # Service implementation
    return ServiceResponse(success=True, data={})
```

4. Add to docker-compose.yml:
```yaml
new-service:
  build: ./bioinformatics-services/new-service
  ports:
    - "8009:8009"
  environment:
    - REDIS_URL=redis://redis:6379
  depends_on:
    - redis
```

## API Documentation

### Workflow API

#### Create Session
```bash
POST /api/v1/sessions
Response: {"session_id": "uuid", "workflow_id": "uuid"}
```

#### Start Workflow
```bash
POST /api/v1/workflows/{workflow_id}/start
Response: {"status": "started", "workflow_id": "uuid"}
```

#### Get Workflow Status
```bash
GET /api/v1/workflows/{workflow_id}
Response: {"workflow": {...}, "steps": [...]}
```

### Tool Execution

#### Execute Tool
```bash
POST /api/v1/tools/{tool_name}/execute
Body: {
  "session_id": "uuid",
  "workflow_id": "uuid", 
  "step_id": "uuid",
  "input_data": {...}
}
```

## Event System

### Event Types

- `workflow_started`: Workflow execution begins
- `step_started`: Individual step begins
- `step_completed`: Step execution completes
- `step_failed`: Step execution fails
- `workflow_completed`: Entire workflow completes

### Event Structure

```json
{
  "session_id": "uuid",
  "workflow_id": "uuid",
  "step_id": "uuid",
  "event_type": "step_completed",
  "data": {...},
  "timestamp": "2024-01-01T00:00:00Z",
  "correlation_id": "uuid"
}
```

## Monitoring and Observability

### Metrics

- Request count and duration
- Error rates
- Active workflows and sessions
- Service health status

### Logging

- Structured logging with correlation IDs
- Service-specific log files
- Centralized log aggregation

### Health Checks

All services expose `/health` endpoints for monitoring:

```bash
curl http://localhost:8000/health
curl http://localhost:8001/health
# ... etc
```

## Security

### Authentication

- JWT-based authentication
- Role-based access control (RBAC)
- Session management

### Data Protection

- Input validation and sanitization
- SQL injection prevention
- XSS protection
- Rate limiting

## Scaling

### Horizontal Scaling

```bash
# Scale specific services
docker-compose up -d --scale alignment-service=3
docker-compose up -d --scale mutation-service=2
```

### Load Balancing

- API Gateway handles load balancing
- Circuit breakers prevent cascade failures
- Retry mechanisms with exponential backoff

## Troubleshooting

### Common Issues

1. **Service not starting**: Check Docker logs
```bash
docker-compose logs service-name
```

2. **Database connection issues**: Verify PostgreSQL is running
```bash
docker-compose ps postgres
```

3. **Redis connection issues**: Check Redis health
```bash
docker exec -it helix-ai-redis-1 redis-cli ping
```

### Debug Mode

```bash
# Run with debug logging
docker-compose up -d --build
docker-compose logs -f
```

## Migration from Monolithic

### Step-by-Step Migration

1. **Phase 1**: Extract shared components
2. **Phase 2**: Create API Gateway
3. **Phase 3**: Implement Workflow Engine
4. **Phase 4**: Extract bioinformatics services
5. **Phase 5**: Add supporting services
6. **Phase 6**: Update frontend

### Data Migration

```sql
-- Migrate existing sessions
INSERT INTO sessions (id, workflow_state)
SELECT session_id, workflow_state FROM old_sessions;

-- Migrate workflows
INSERT INTO workflows (session_id, name, definition, status)
SELECT session_id, 'migrated_workflow', '{}', 'completed' FROM sessions;
```

## Performance Optimization

### Caching Strategy

- Redis for session data
- Service response caching
- Database query caching

### Database Optimization

- Connection pooling
- Query optimization
- Index optimization

### Service Optimization

- Async/await patterns
- Background task processing
- Resource pooling

## Contributing

1. Fork the repository
2. Create feature branch
3. Implement changes
4. Add tests
5. Submit pull request

## License

MIT License - see LICENSE file for details. 