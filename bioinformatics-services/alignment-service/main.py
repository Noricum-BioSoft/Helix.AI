from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import Dict, Any, List, Optional
import os
import time
import json
from contextlib import asynccontextmanager

from shared.models import ServiceRequest, ServiceResponse, HealthCheck
from shared.utils import setup_logging, track_request_metrics, get_redis_client, EventBus

# Setup logging
logger = setup_logging("alignment-service")

# Event bus
event_bus = EventBus()

class AlignmentRequest(BaseModel):
    sequences: str
    algorithm: Optional[str] = "clustalw"
    gap_penalty: Optional[float] = 10.0
    match_score: Optional[float] = 1.0
    mismatch_penalty: Optional[float] = -1.0

class AlignmentResult(BaseModel):
    aligned_sequences: List[Dict[str, str]]
    alignment_score: float
    consensus_sequence: str
    metadata: Dict[str, Any]

@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    logger.info("Starting Alignment Service")
    
    # Register this service
    redis_client = get_redis_client()
    redis_client.hset("services", "alignment-service", "http://localhost:8002")
    
    yield
    
    # Shutdown
    logger.info("Shutting down Alignment Service")

app = FastAPI(
    title="DataBloom.AI Alignment Service",
    description="Sequence alignment service for bioinformatics workflows",
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
    return HealthCheck(
        status="healthy",
        service="alignment-service",
        details={"timestamp": time.time()}
    )

# Service information
@app.get("/info")
async def get_service_info():
    """Get service information"""
    return {
        "name": "alignment-service",
        "version": "1.0.0",
        "description": "Sequence alignment service",
        "capabilities": [
            "Multiple sequence alignment",
            "Pairwise alignment",
            "Consensus sequence generation",
            "Alignment scoring"
        ],
        "supported_algorithms": [
            "clustalw",
            "muscle",
            "tcoffee",
            "mafft"
        ]
    }

# Main execution endpoint
@app.post("/execute")
async def execute_alignment(request: ServiceRequest):
    """Execute sequence alignment"""
    try:
        logger.info(f"Executing alignment for step {request.step_id}")
        
        # Parse input data
        input_data = request.input_data
        sequences_text = input_data.get("sequences", "")
        algorithm = input_data.get("algorithm", "clustalw")
        gap_penalty = input_data.get("gap_penalty", 10.0)
        match_score = input_data.get("match_score", 1.0)
        mismatch_penalty = input_data.get("mismatch_penalty", -1.0)
        
        # Parse sequences
        sequences = parse_sequences(sequences_text)
        if not sequences:
            raise ValueError("No valid sequences provided")
        
        # Perform alignment
        alignment_result = await perform_alignment(
            sequences, algorithm, gap_penalty, match_score, mismatch_penalty
        )
        
        # Store result in session
        await store_alignment_result(request.session_id, alignment_result)
        
        # Publish event
        await event_bus.publish_event("alignment_events", {
            "session_id": request.session_id,
            "workflow_id": request.workflow_id,
            "step_id": request.step_id,
            "event_type": "alignment_completed",
            "data": alignment_result.dict()
        })
        
        return ServiceResponse(
            success=True,
            data=alignment_result.dict(),
            correlation_id=request.correlation_id
        )
        
    except Exception as e:
        logger.error(f"Error executing alignment: {e}")
        return ServiceResponse(
            success=False,
            error=str(e),
            correlation_id=request.correlation_id
        )

def parse_sequences(sequences_text: str) -> List[Dict[str, str]]:
    """Parse sequences from text format"""
    sequences = []
    current_name = None
    current_sequence = ""
    
    for line in sequences_text.split('\n'):
        line = line.strip()
        if not line:
            continue
            
        if line.startswith('>'):
            # Save previous sequence
            if current_name and current_sequence:
                sequences.append({
                    "name": current_name,
                    "sequence": current_sequence
                })
            
            # Start new sequence
            current_name = line[1:].strip()
            current_sequence = ""
        else:
            # Add to current sequence
            current_sequence += line.upper()
    
    # Add last sequence
    if current_name and current_sequence:
        sequences.append({
            "name": current_name,
            "sequence": current_sequence
        })
    
    return sequences

async def perform_alignment(
    sequences: List[Dict[str, str]], 
    algorithm: str, 
    gap_penalty: float,
    match_score: float,
    mismatch_penalty: float
) -> AlignmentResult:
    """Perform sequence alignment"""
    try:
        # Simple alignment algorithm (in production, use BioPython or similar)
        aligned_sequences = simple_multiple_alignment(sequences, gap_penalty, match_score, mismatch_penalty)
        
        # Calculate consensus
        consensus = calculate_consensus(aligned_sequences)
        
        # Calculate alignment score
        score = calculate_alignment_score(aligned_sequences, match_score, mismatch_penalty, gap_penalty)
        
        return AlignmentResult(
            aligned_sequences=aligned_sequences,
            alignment_score=score,
            consensus_sequence=consensus,
            metadata={
                "algorithm": algorithm,
                "num_sequences": len(sequences),
                "gap_penalty": gap_penalty,
                "match_score": match_score,
                "mismatch_penalty": mismatch_penalty
            }
        )
    except Exception as e:
        logger.error(f"Error in alignment: {e}")
        raise e

def simple_multiple_alignment(sequences: List[Dict[str, str]], gap_penalty: float, match_score: float, mismatch_penalty: float) -> List[Dict[str, str]]:
    """Simple multiple sequence alignment algorithm"""
    if len(sequences) == 0:
        return []
    
    if len(sequences) == 1:
        return sequences
    
    # Find the longest sequence as reference
    reference = max(sequences, key=lambda x: len(x["sequence"]))
    ref_seq = reference["sequence"]
    
    aligned_sequences = []
    
    for seq in sequences:
        if seq["sequence"] == ref_seq:
            # Reference sequence doesn't need alignment
            aligned_sequences.append({
                "name": seq["name"],
                "sequence": ref_seq
            })
        else:
            # Simple pairwise alignment with reference
            aligned_seq = simple_pairwise_alignment(ref_seq, seq["sequence"], gap_penalty, match_score, mismatch_penalty)
            aligned_sequences.append({
                "name": seq["name"],
                "sequence": aligned_seq
            })
    
    return aligned_sequences

def simple_pairwise_alignment(seq1: str, seq2: str, gap_penalty: float, match_score: float, mismatch_penalty: float) -> str:
    """Simple pairwise alignment"""
    # Simple algorithm: align shorter sequence to longer one
    if len(seq1) >= len(seq2):
        # Pad seq2 with gaps
        padding = len(seq1) - len(seq2)
        return seq2 + "-" * padding
    else:
        # Pad seq1 with gaps
        padding = len(seq2) - len(seq1)
        return seq1 + "-" * padding

def calculate_consensus(aligned_sequences: List[Dict[str, str]]) -> str:
    """Calculate consensus sequence"""
    if not aligned_sequences:
        return ""
    
    # Get the length of aligned sequences
    seq_length = len(aligned_sequences[0]["sequence"])
    consensus = ""
    
    for i in range(seq_length):
        # Count bases at this position
        bases = {}
        for seq in aligned_sequences:
            base = seq["sequence"][i]
            if base != "-":
                bases[base] = bases.get(base, 0) + 1
        
        if bases:
            # Find most common base
            most_common = max(bases.items(), key=lambda x: x[1])[0]
            consensus += most_common
        else:
            consensus += "-"
    
    return consensus

def calculate_alignment_score(aligned_sequences: List[Dict[str, str]], match_score: float, mismatch_penalty: float, gap_penalty: float) -> float:
    """Calculate alignment score"""
    if not aligned_sequences:
        return 0.0
    
    total_score = 0.0
    seq_length = len(aligned_sequences[0]["sequence"])
    
    for i in range(seq_length):
        # Get all bases at this position
        bases = [seq["sequence"][i] for seq in aligned_sequences]
        
        # Calculate score for this position
        for j in range(len(bases)):
            for k in range(j + 1, len(bases)):
                if bases[j] == "-" or bases[k] == "-":
                    total_score += gap_penalty
                elif bases[j] == bases[k]:
                    total_score += match_score
                else:
                    total_score += mismatch_penalty
    
    return total_score

async def store_alignment_result(session_id: str, result: AlignmentResult):
    """Store alignment result in session"""
    try:
        redis_client = get_redis_client()
        key = f"session:{session_id}:alignment"
        redis_client.set(key, json.dumps(result.dict()))
        logger.info(f"Stored alignment result for session {session_id}")
    except Exception as e:
        logger.error(f"Error storing alignment result: {e}")

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8002) 