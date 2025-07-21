from fastapi import FastAPI, HTTPException, Depends
from fastapi.middleware.cors import CORSMiddleware
from fastapi.middleware.trustedhost import TrustedHostMiddleware
import httpx
import asyncio
from typing import Dict, Any
import os
import time
from contextlib import asynccontextmanager

from shared.models import Session, Workflow, ServiceRequest, HealthCheck
from shared.utils import get_redis_client, EventBus, ServiceRegistry, CircuitBreaker, setup_logging, track_request_metrics

# Setup logging
logger = setup_logging("api-gateway")

# Circuit breakers for each service
circuit_breakers = {
    "workflow-engine": CircuitBreaker(),
    "alignment-service": CircuitBreaker(),
    "mutation-service": CircuitBreaker(),
    "selection-service": CircuitBreaker(),
    "plasmid-service": CircuitBreaker(),
    "synthesis-service": CircuitBreaker(),
    "nlp-service": CircuitBreaker(),
}

# Service registry
service_registry = ServiceRegistry()
event_bus = EventBus()

@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    logger.info("Starting API Gateway")
    
    # Register this service
    await service_registry.register_service("api-gateway", "http://localhost:8000")
    
    yield
    
    # Shutdown
    logger.info("Shutting down API Gateway")

app = FastAPI(
    title="Helix.AI API Gateway",
    description="API Gateway for Helix.AI microservices",
    version="1.0.0",
    lifespan=lifespan
)

# Middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

app.add_middleware(TrustedHostMiddleware, allowed_hosts=["*"])

# Request timing middleware
@app.middleware("http")
async def add_process_time_header(request, call_next):
    start_time = time.time()
    response = await call_next(request)
    process_time = time.time() - start_time
    response.headers["X-Process-Time"] = str(process_time)
    
    # Track metrics
    track_request_metrics(
        method=request.method,
        endpoint=request.url.path,
        status=response.status_code,
        duration=process_time
    )
    
    return response

# Health check
@app.get("/health")
async def health_check():
    """Health check endpoint"""
    return HealthCheck(
        status="healthy",
        service="api-gateway",
        details={"timestamp": time.time()}
    )

# Session management
@app.post("/api/v1/sessions")
async def create_session():
    """Create a new session"""
    try:
        workflow_url = await service_registry.get_service_url("workflow-engine")
        if not workflow_url:
            raise HTTPException(status_code=503, detail="Workflow engine not available")
        
        async with httpx.AsyncClient() as client:
            response = await circuit_breakers["workflow-engine"].call(
                client.post, f"{workflow_url}/sessions"
            )
            return response.json()
    except Exception as e:
        logger.error(f"Error creating session: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/v1/sessions/{session_id}")
async def get_session(session_id: str):
    """Get session information"""
    try:
        workflow_url = await service_registry.get_service_url("workflow-engine")
        if not workflow_url:
            raise HTTPException(status_code=503, detail="Workflow engine not available")
        
        async with httpx.AsyncClient() as client:
            response = await circuit_breakers["workflow-engine"].call(
                client.get, f"{workflow_url}/sessions/{session_id}"
            )
            return response.json()
    except Exception as e:
        logger.error(f"Error getting session: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Workflow management
@app.post("/api/v1/workflows")
async def create_workflow(workflow: Workflow):
    """Create a new workflow"""
    try:
        workflow_url = await service_registry.get_service_url("workflow-engine")
        if not workflow_url:
            raise HTTPException(status_code=503, detail="Workflow engine not available")
        
        async with httpx.AsyncClient() as client:
            response = await circuit_breakers["workflow-engine"].call(
                client.post, f"{workflow_url}/workflows", json=workflow.dict()
            )
            return response.json()
    except Exception as e:
        logger.error(f"Error creating workflow: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/api/v1/workflows/{workflow_id}")
async def get_workflow(workflow_id: str):
    """Get workflow information"""
    try:
        workflow_url = await service_registry.get_service_url("workflow-engine")
        if not workflow_url:
            raise HTTPException(status_code=503, detail="Workflow engine not available")
        
        async with httpx.AsyncClient() as client:
            response = await circuit_breakers["workflow-engine"].call(
                client.get, f"{workflow_url}/workflows/{workflow_id}"
            )
            return response.json()
    except Exception as e:
        logger.error(f"Error getting workflow: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/api/v1/workflows/{workflow_id}/start")
async def start_workflow(workflow_id: str):
    """Start a workflow"""
    try:
        workflow_url = await service_registry.get_service_url("workflow-engine")
        if not workflow_url:
            raise HTTPException(status_code=503, detail="Workflow engine not available")
        
        async with httpx.AsyncClient() as client:
            response = await circuit_breakers["workflow-engine"].call(
                client.post, f"{workflow_url}/workflows/{workflow_id}/start"
            )
            return response.json()
    except Exception as e:
        logger.error(f"Error starting workflow: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Tool execution
@app.post("/api/v1/tools/{tool_name}/execute")
async def execute_tool(tool_name: str, request: ServiceRequest):
    """Execute a specific tool"""
    try:
        # Map tool names to services
        service_mapping = {
            "sequence_alignment": "alignment-service",
            "sequence_selection": "selection-service",
            "mutate_sequence": "mutation-service",
            "plasmid_visualization": "plasmid-service",
            "synthesis_submission": "synthesis-service",
        }
        
        service_name = service_mapping.get(tool_name)
        if not service_name:
            raise HTTPException(status_code=404, detail=f"Tool {tool_name} not found")
        
        service_url = await service_registry.get_service_url(service_name)
        if not service_url:
            raise HTTPException(status_code=503, detail=f"Service {service_name} not available")
        
        async with httpx.AsyncClient() as client:
            response = await circuit_breakers[service_name].call(
                client.post, f"{service_url}/execute", json=request.dict()
            )
            return response.json()
    except Exception as e:
        logger.error(f"Error executing tool {tool_name}: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Natural language processing
@app.post("/api/v1/nlp/parse")
async def parse_natural_language(command: str, session_id: str):
    """Parse natural language command"""
    try:
        nlp_url = await service_registry.get_service_url("nlp-service")
        if not nlp_url:
            raise HTTPException(status_code=503, detail="NLP service not available")
        
        async with httpx.AsyncClient() as client:
            response = await circuit_breakers["nlp-service"].call(
                client.post, f"{nlp_url}/parse", json={"command": command, "session_id": session_id}
            )
            return response.json()
    except Exception as e:
        logger.error(f"Error parsing natural language: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Service discovery
@app.get("/api/v1/services")
async def list_services():
    """List all available services"""
    try:
        services = await service_registry.list_services()
        return {"services": services}
    except Exception as e:
        logger.error(f"Error listing services: {e}")
        raise HTTPException(status_code=500, detail=str(e))

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000) 