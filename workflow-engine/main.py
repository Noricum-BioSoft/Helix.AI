from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
import asyncio
import httpx
from typing import Dict, Any, List
import os
import time
from contextlib import asynccontextmanager

from shared.models import Session, Workflow, WorkflowStep, WorkflowEvent, WorkflowStatus, StepStatus, ServiceRequest, ServiceResponse
from shared.utils import get_redis_client, EventBus, ServiceRegistry, CircuitBreaker, setup_logging, track_request_metrics, get_db_connection

# Setup logging
logger = setup_logging("workflow-engine")

# Circuit breakers for services
circuit_breakers = {
    "alignment-service": CircuitBreaker(),
    "mutation-service": CircuitBreaker(),
    "selection-service": CircuitBreaker(),
    "plasmid-service": CircuitBreaker(),
    "synthesis-service": CircuitBreaker(),
}

# Service registry and event bus
service_registry = ServiceRegistry()
event_bus = EventBus()

# Workflow definitions
STANDARD_WORKFLOW = {
    "name": "standard_bioinformatics_pipeline",
    "steps": [
        {
            "name": "align",
            "step_type": "sequence_alignment",
            "service_name": "alignment-service",
            "description": "Align multiple DNA/RNA sequences"
        },
        {
            "name": "select",
            "step_type": "sequence_selection",
            "service_name": "selection-service",
            "depends_on": "align",
            "description": "Select sequences from alignment"
        },
        {
            "name": "mutate",
            "step_type": "mutate_sequence",
            "service_name": "mutation-service",
            "depends_on": "select",
            "description": "Generate sequence variants"
        },
        {
            "name": "visualize",
            "step_type": "plasmid_visualization",
            "service_name": "plasmid-service",
            "depends_on": "mutate",
            "description": "Visualize plasmid constructs"
        },
        {
            "name": "synthesize",
            "step_type": "synthesis_submission",
            "service_name": "synthesis-service",
            "depends_on": "mutate",
            "description": "Submit for DNA synthesis"
        }
    ]
}

@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    logger.info("Starting Workflow Engine")
    
    # Register this service
    await service_registry.register_service("workflow-engine", "http://localhost:8001")
    
    # Start event listener
    asyncio.create_task(event_listener())
    
    yield
    
    # Shutdown
    logger.info("Shutting down Workflow Engine")

app = FastAPI(
    title="DataBloom.AI Workflow Engine",
    description="Workflow orchestration engine for bioinformatics pipelines",
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
    return {
        "status": "healthy",
        "service": "workflow-engine",
        "timestamp": time.time()
    }

# Session management
@app.post("/sessions")
async def create_session():
    """Create a new session"""
    try:
        async with get_db_connection() as conn:
            session_id = await conn.fetchval(
                "INSERT INTO sessions (id) VALUES (gen_random_uuid()) RETURNING id"
            )
            
            # Create default workflow
            workflow_id = await conn.fetchval(
                """
                INSERT INTO workflows (session_id, name, definition, status)
                VALUES ($1, $2, $3, $4)
                RETURNING id
                """,
                session_id, "standard_pipeline", STANDARD_WORKFLOW, WorkflowStatus.PENDING
            )
            
            # Create workflow steps
            for step in STANDARD_WORKFLOW["steps"]:
                await conn.execute(
                    """
                    INSERT INTO workflow_steps (workflow_id, step_name, step_type, service_name)
                    VALUES ($1, $2, $3, $4)
                    """,
                    workflow_id, step["name"], step["step_type"], step["service_name"]
                )
            
            return {"session_id": str(session_id), "workflow_id": str(workflow_id)}
    except Exception as e:
        logger.error(f"Error creating session: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/sessions/{session_id}")
async def get_session(session_id: str):
    """Get session information"""
    try:
        async with get_db_connection() as conn:
            session = await conn.fetchrow(
                "SELECT * FROM sessions WHERE id = $1",
                session_id
            )
            
            if not session:
                raise HTTPException(status_code=404, detail="Session not found")
            
            workflows = await conn.fetch(
                "SELECT * FROM workflows WHERE session_id = $1",
                session_id
            )
            
            return {
                "session": dict(session),
                "workflows": [dict(w) for w in workflows]
            }
    except Exception as e:
        logger.error(f"Error getting session: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Workflow management
@app.post("/workflows")
async def create_workflow(workflow: Workflow):
    """Create a new workflow"""
    try:
        async with get_db_connection() as conn:
            workflow_id = await conn.fetchval(
                """
                INSERT INTO workflows (session_id, name, definition, status)
                VALUES ($1, $2, $3, $4)
                RETURNING id
                """,
                workflow.session_id, workflow.name, workflow.definition, workflow.status
            )
            
            # Create workflow steps
            for step in workflow.definition.get("steps", []):
                await conn.execute(
                    """
                    INSERT INTO workflow_steps (workflow_id, step_name, step_type, service_name)
                    VALUES ($1, $2, $3, $4)
                    """,
                    workflow_id, step["name"], step["step_type"], step["service_name"]
                )
            
            return {"workflow_id": str(workflow_id)}
    except Exception as e:
        logger.error(f"Error creating workflow: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get("/workflows/{workflow_id}")
async def get_workflow(workflow_id: str):
    """Get workflow information"""
    try:
        async with get_db_connection() as conn:
            workflow = await conn.fetchrow(
                "SELECT * FROM workflows WHERE id = $1",
                workflow_id
            )
            
            if not workflow:
                raise HTTPException(status_code=404, detail="Workflow not found")
            
            steps = await conn.fetch(
                "SELECT * FROM workflow_steps WHERE workflow_id = $1 ORDER BY created_at",
                workflow_id
            )
            
            return {
                "workflow": dict(workflow),
                "steps": [dict(s) for s in steps]
            }
    except Exception as e:
        logger.error(f"Error getting workflow: {e}")
        raise HTTPException(status_code=500, detail=str(e))

@app.post("/workflows/{workflow_id}/start")
async def start_workflow(workflow_id: str):
    """Start a workflow"""
    try:
        async with get_db_connection() as conn:
            # Update workflow status
            await conn.execute(
                "UPDATE workflows SET status = $1 WHERE id = $2",
                WorkflowStatus.RUNNING, workflow_id
            )
            
            # Get workflow steps
            steps = await conn.fetch(
                "SELECT * FROM workflow_steps WHERE workflow_id = $1 ORDER BY created_at",
                workflow_id
            )
            
            # Start first step
            if steps:
                first_step = steps[0]
                await execute_workflow_step(str(first_step["id"]))
            
            return {"status": "started", "workflow_id": workflow_id}
    except Exception as e:
        logger.error(f"Error starting workflow: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# Step execution
async def execute_workflow_step(step_id: str):
    """Execute a workflow step"""
    try:
        async with get_db_connection() as conn:
            # Get step information
            step = await conn.fetchrow(
                "SELECT * FROM workflow_steps WHERE id = $1",
                step_id
            )
            
            if not step:
                raise Exception(f"Step {step_id} not found")
            
            # Update step status to running
            await conn.execute(
                """
                UPDATE workflow_steps 
                SET status = $1, started_at = NOW() 
                WHERE id = $2
                """,
                StepStatus.RUNNING, step_id
            )
            
            # Get service URL
            service_url = await service_registry.get_service_url(step["service_name"])
            if not service_url:
                raise Exception(f"Service {step['service_name']} not available")
            
            # Prepare request
            request = ServiceRequest(
                session_id=step["workflow_id"],  # Using workflow_id as session_id for now
                workflow_id=step["workflow_id"],
                step_id=step_id,
                input_data=step["input_data"] or {}
            )
            
            # Execute step
            async with httpx.AsyncClient() as client:
                response = await circuit_breakers[step["service_name"]].call(
                    client.post, f"{service_url}/execute", json=request.dict()
                )
                
                result = response.json()
                
                # Update step with result
                await conn.execute(
                    """
                    UPDATE workflow_steps 
                    SET status = $1, output_data = $2, completed_at = NOW()
                    WHERE id = $3
                    """,
                    StepStatus.COMPLETED, result, step_id
                )
                
                # Publish event
                await event_bus.publish_event("workflow_events", {
                    "session_id": step["workflow_id"],
                    "workflow_id": step["workflow_id"],
                    "step_id": step_id,
                    "event_type": "step_completed",
                    "data": result
                })
                
                # Check if there are more steps to execute
                await check_and_execute_next_step(step["workflow_id"], step_id)
                
                return result
    except Exception as e:
        logger.error(f"Error executing workflow step {step_id}: {e}")
        
        # Update step status to failed
        async with get_db_connection() as conn:
            await conn.execute(
                """
                UPDATE workflow_steps 
                SET status = $1, error_message = $2, completed_at = NOW()
                WHERE id = $3
                """,
                StepStatus.FAILED, str(e), step_id
            )
        
        raise e

async def check_and_execute_next_step(workflow_id: str, current_step_id: str):
    """Check if there are more steps to execute and execute the next one"""
    try:
        async with get_db_connection() as conn:
            # Get all steps for this workflow
            steps = await conn.fetch(
                "SELECT * FROM workflow_steps WHERE workflow_id = $1 ORDER BY created_at",
                workflow_id
            )
            
            # Find current step index
            current_index = None
            for i, step in enumerate(steps):
                if str(step["id"]) == current_step_id:
                    current_index = i
                    break
            
            if current_index is None or current_index >= len(steps) - 1:
                # No more steps, workflow is complete
                await conn.execute(
                    "UPDATE workflows SET status = $1 WHERE id = $2",
                    WorkflowStatus.COMPLETED, workflow_id
                )
                return
            
            # Execute next step
            next_step = steps[current_index + 1]
            await execute_workflow_step(str(next_step["id"]))
    except Exception as e:
        logger.error(f"Error checking next step: {e}")

# Event listener
async def event_listener():
    """Listen for workflow events"""
    async def handle_event(event_data):
        logger.info(f"Received event: {event_data}")
        
        if event_data.get("event_type") == "step_completed":
            # Handle step completion
            workflow_id = event_data.get("workflow_id")
            step_id = event_data.get("step_id")
            
            if workflow_id and step_id:
                await check_and_execute_next_step(workflow_id, step_id)
    
    await event_bus.subscribe_to_events("workflow_events", handle_event)

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001) 