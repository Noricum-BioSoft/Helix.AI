#!/usr/bin/env python3
"""
Generate a sample OrchestratorTrace JSON without requiring LLM dependencies.

This creates a realistic trace showing the full structure including hashes,
timing, and complete contracts.
"""

import json
import hashlib
import time
from datetime import datetime


def compute_contract_hash(contract_json: str) -> str:
    """Compute SHA256 hash of contract JSON."""
    if not contract_json:
        return ""
    hash_obj = hashlib.sha256(contract_json.encode('utf-8'))
    return hash_obj.hexdigest()


def generate_sample_trace():
    """Generate a realistic OrchestratorTrace."""
    
    # Base timestamps
    start_time = time.time()
    
    # Sample workflow plan (input to InfrastructureDecisionAgent)
    workflow_plan = {
        "description": "Quality control analysis of small RNA-seq data using FastQC",
        "operations": [
            {
                "operation_name": "fastqc",
                "tool_name": "fastqc",
                "parameters": {
                    "threads": 2,
                    "output_format": "html"
                }
            }
        ],
        "data_inputs": [
            {
                "uri": "/data/rnaseq_demo/sample_R1_realistic.fastq",
                "size_bytes": 50000000,
                "location_type": "Local",
                "metadata": {
                    "read": "R1",
                    "format": "fastq",
                    "quality_encoding": "Phred33"
                }
            },
            {
                "uri": "/data/rnaseq_demo/sample_R2_realistic.fastq",
                "size_bytes": 48000000,
                "location_type": "Local",
                "metadata": {
                    "read": "R2",
                    "format": "fastq",
                    "quality_encoding": "Phred33"
                }
            }
        ],
        "expected_outputs": [
            {
                "uri": "/output/sample_R1_fastqc.html",
                "type": "html",
                "description": "FastQC quality report for R1"
            },
            {
                "uri": "/output/sample_R2_fastqc.html",
                "type": "html",
                "description": "FastQC quality report for R2"
            }
        ]
    }
    
    # InfraDecision contract (output from InfrastructureDecisionAgent)
    infra_decision = {
        "infrastructure": "Local",
        "reasoning": "Small local files (total 98MB) with no network transfer costs. Local execution is optimal for this workload. Files are already on the local filesystem, avoiding S3 transfer latency and costs. FastQC has modest memory requirements (<2GB) that can be satisfied locally. Estimated runtime: 2-3 minutes on modern CPU.",
        "confidence_score": 0.92,
        "warnings": [],
        "alternatives": [
            {
                "environment": "EC2",
                "reasoning": "Over-provisioned for small files; would incur EC2 startup overhead (~30s) and instance costs",
                "relative_cost": 1.5,
                "tradeoffs": "Faster CPU (t3.medium) but unnecessary network transfer and instance startup time. Total time would be similar or worse."
            },
            {
                "environment": "EMR",
                "reasoning": "Significantly over-provisioned; EMR designed for distributed processing of large datasets (>10GB)",
                "relative_cost": 3.0,
                "tradeoffs": "EMR cluster startup (5-10 min) far exceeds local execution time. Only beneficial for data >100GB."
            }
        ],
        "cost_analysis": {
            "estimated_cost_range_usd": [0.0, 0.01],
            "cost_assumptions": "Local execution has negligible cost (only electricity, ~$0.001/hour). Based on estimated 2-3 minute runtime with FastQC on 2 cores. No cloud infrastructure costs (EC2, S3 transfer, EMR).",
            "cost_confidence": 0.95
        },
        "resource_requirements": {
            "min_cpu_cores": 2,
            "min_memory_gb": 2,
            "min_storage_gb": 1,
            "estimated_runtime_minutes": 2,
            "notes": "FastQC memory usage scales with file size; 2GB sufficient for 50MB FASTQ files. Runtime: ~1min per file with 2 threads."
        }
    }
    
    # ExecutionToolSpec contract (output from ImplementationAgent)
    execution_spec = {
        "tool_name": "fastqc",
        "infrastructure": "Local",
        "commands": [
            {
                "command": "fastqc --threads 2 --outdir /output /data/rnaseq_demo/sample_R1_realistic.fastq /data/rnaseq_demo/sample_R2_realistic.fastq",
                "description": "Run FastQC quality control on paired-end FASTQ files with 2 threads",
                "timeout_seconds": 300,
                "retry_policy": {
                    "max_retries": 2,
                    "backoff_multiplier": 2.0,
                    "retry_on_exit_codes": [1, 2, 137]
                }
            }
        ],
        "container_spec": {
            "image": "biocontainers/fastqc:v0.11.9_cv8",
            "entrypoint": None,
            "working_dir": "/workspace",
            "environment": {
                "FASTQC_THREADS": "2",
                "JAVA_OPTS": "-Xmx2g"
            },
            "volumes": [
                "/data/rnaseq_demo:/data/rnaseq_demo:ro",
                "/output:/output:rw"
            ]
        },
        "inputs": [
            {
                "uri": "/data/rnaseq_demo/sample_R1_realistic.fastq",
                "type": "fastq",
                "size_bytes": 50000000,
                "source": "Local",
                "required": True
            },
            {
                "uri": "/data/rnaseq_demo/sample_R2_realistic.fastq",
                "type": "fastq",
                "size_bytes": 48000000,
                "source": "Local",
                "required": True
            }
        ],
        "expected_outputs": [
            {
                "uri": "/output/sample_R1_realistic_fastqc.html",
                "type": "html",
                "description": "FastQC quality report for R1",
                "required": True
            },
            {
                "uri": "/output/sample_R2_realistic_fastqc.html",
                "type": "html",
                "description": "FastQC quality report for R2",
                "required": True
            },
            {
                "uri": "/output/sample_R1_realistic_fastqc.zip",
                "type": "zip",
                "description": "Detailed FastQC data archive for R1",
                "required": False
            },
            {
                "uri": "/output/sample_R2_realistic_fastqc.zip",
                "type": "zip",
                "description": "Detailed FastQC data archive for R2",
                "required": False
            }
        ],
        "validation_criteria": {
            "required_outputs": [
                "/output/sample_R1_realistic_fastqc.html",
                "/output/sample_R2_realistic_fastqc.html"
            ],
            "success_conditions": [
                "Exit code 0",
                "All required output files exist",
                "HTML reports contain 'FastQC Report'",
                "No 'FAIL' markers in summary.txt (if present)"
            ]
        },
        "resource_requirements": {
            "min_cpu_cores": 2,
            "min_memory_gb": 2,
            "min_storage_gb": 1,
            "estimated_runtime_minutes": 2
        },
        "metadata": {
            "tool_version": "0.11.9",
            "container_registry": "biocontainers",
            "documentation_url": "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/",
            "citation": "Andrews, S. (2010). FastQC: a quality control tool for high throughput sequence data."
        },
        "confidence_score": 0.88,
        "warnings": [
            "Using latest stable version (0.11.9); consider pinning specific version for reproducibility"
        ]
    }
    
    # Compute hashes
    workflow_plan_json = json.dumps(workflow_plan, sort_keys=True)
    infra_decision_json = json.dumps(infra_decision, sort_keys=True)
    execution_spec_json = json.dumps(execution_spec, sort_keys=True)
    
    infra_input_hash = compute_contract_hash(workflow_plan_json)
    infra_output_hash = compute_contract_hash(infra_decision_json)
    
    impl_input_json = json.dumps({
        "workflow_plan": workflow_plan,
        "infra_decision": infra_decision
    }, sort_keys=True)
    impl_input_hash = compute_contract_hash(impl_input_json)
    impl_output_hash = compute_contract_hash(execution_spec_json)
    
    # Build trace
    request_id = f"demo_small_local_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    
    # InfrastructureDecisionAgent invocation (2.3s)
    infra_start = start_time
    infra_end = infra_start + 2.333
    
    infra_invocation = {
        "agent_name": "InfrastructureDecisionAgent",
        "request_id": request_id,
        "start_time": infra_start,
        "end_time": infra_end,
        "duration_ms": (infra_end - infra_start) * 1000,
        "success": True,
        "error": None,
        "input_hash": infra_input_hash,
        "output_hash": infra_output_hash,
        "confidence_score": 0.92,
        "warnings": [],
        "timestamp": datetime.fromtimestamp(infra_start).isoformat()
    }
    
    # ImplementationAgent invocation (3.2s)
    impl_start = infra_end + 0.111
    impl_end = impl_start + 3.221
    
    impl_invocation = {
        "agent_name": "ImplementationAgent",
        "request_id": request_id,
        "start_time": impl_start,
        "end_time": impl_end,
        "duration_ms": (impl_end - impl_start) * 1000,
        "success": True,
        "error": None,
        "input_hash": impl_input_hash,
        "output_hash": impl_output_hash,
        "confidence_score": 0.88,
        "warnings": [
            "Using latest stable version (0.11.9); consider pinning specific version for reproducibility"
        ],
        "timestamp": datetime.fromtimestamp(impl_start).isoformat()
    }
    
    # Complete trace
    end_time = impl_end
    
    trace = {
        "request_id": request_id,
        "user_command": "Run FastQC quality control on local FASTQ files",
        "start_time": start_time,
        "end_time": end_time,
        "duration_ms": (end_time - start_time) * 1000,
        "success": True,
        "error": None,
        "timestamp": datetime.fromtimestamp(start_time).isoformat(),
        "invocations": [infra_invocation, impl_invocation],
        "infra_decision": infra_decision,
        "execution_spec": execution_spec
    }
    
    return trace


def main():
    print("=" * 70)
    print("GENERATING SAMPLE ORCHESTRATOR TRACE")
    print("=" * 70)
    print()
    
    trace = generate_sample_trace()
    
    # Save to file
    output_path = "docs/SAMPLE_TRACE_GENERATED.json"
    with open(output_path, 'w') as f:
        json.dump(trace, f, indent=2)
    
    print(f"✅ Generated trace saved to: {output_path}")
    print()
    print("📊 Trace Summary:")
    print("=" * 70)
    print(f"Request ID: {trace['request_id']}")
    print(f"Command: {trace['user_command']}")
    print(f"Duration: {trace['duration_ms']/1000:.2f}s")
    print(f"Success: {trace['success']}")
    print()
    print(f"Agents Called: {len(trace['invocations'])}")
    for inv in trace['invocations']:
        print(f"  - {inv['agent_name']}: {inv['duration_ms']/1000:.2f}s (confidence={inv['confidence_score']})")
        print(f"    Input hash:  {inv['input_hash'][:16]}...")
        print(f"    Output hash: {inv['output_hash'][:16]}...")
        if inv['warnings']:
            print(f"    Warnings: {len(inv['warnings'])}")
    print()
    print(f"Infrastructure: {trace['infra_decision']['infrastructure']}")
    print(f"Cost Range: ${trace['infra_decision']['cost_analysis']['estimated_cost_range_usd'][0]:.2f}-${trace['infra_decision']['cost_analysis']['estimated_cost_range_usd'][1]:.2f}")
    print(f"Tool: {trace['execution_spec']['tool_name']}")
    print(f"Container: {trace['execution_spec']['container_spec']['image']}")
    print()
    print("=" * 70)
    print("✅ Sample trace generated successfully!")
    print()
    print("This trace demonstrates:")
    print("  ✅ Provenance hashing (P1-1 fix) - SHA256 for all contracts")
    print("  ✅ Complete observability - timing, confidence, warnings")
    print("  ✅ Full contracts - InfraDecision + ExecutionToolSpec")
    print("  ✅ Baseline testing ready - hashes enable fast comparison")
    print()


if __name__ == "__main__":
    main()
