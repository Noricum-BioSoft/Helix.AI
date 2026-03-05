#!/usr/bin/env python3
"""
CLI Demo for Multi-Agent Refactored Architecture (Phases 1-4)

Demonstrates the full pipeline:
1. WorkflowPlan (input contract)
2. Infrastructure Decision Agent → InfraDecision
3. Implementation Agent → ExecutionToolSpec
4. Orchestrator coordination
5. Phase 2 tools integration (FileMetadataInspector, etc.)

Usage:
    python backend/cli_demo.py --scenario small_local
    python backend/cli_demo.py --scenario large_s3
    python backend/cli_demo.py --scenario unknown_sizes
    python backend/cli_demo.py --custom --command "Run FastQC on my data" --files "s3://bucket/file.fq"
"""

import sys
import asyncio
import argparse
import json
from pathlib import Path
from typing import Dict, Any, Optional
import logging

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from backend.contracts.workflow_plan import WorkflowPlan, OperationSpec, DataInput
from backend.contracts.dataset_spec import DatasetSpec
from backend.orchestrator import Orchestrator
from backend.tools import FileMetadataInspector, EnvironmentCapabilityCatalog, CostHeuristicTable

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


# ============================================================================
# Predefined Scenarios
# ============================================================================

SCENARIOS = {
    "small_local": {
        "description": "Small local files (<100MB) → Expect Local execution",
        "command": "Run FastQC quality control on local FASTQ files",
        "workflow_plan": {
            "description": "Quality control analysis of small RNA-seq data using FastQC",
            "operations": [
                {
                    "operation_name": "fastqc",
                    "tool_name": "fastqc"
                }
            ],
            "data_inputs": [
                {
                    "uri": "/data/sample_R1.fastq",
                    "size_bytes": 50000000,  # 50MB
                    "location_type": "Local",
                    "metadata": {"read": "R1", "format": "fastq"}
                },
                {
                    "uri": "/data/sample_R2.fastq",
                    "size_bytes": 48000000,  # 48MB
                    "location_type": "Local",
                    "metadata": {"read": "R2", "format": "fastq"}
                }
            ]
        },
        "expected": {
            "infrastructure": "Local",
            "confidence_range": [0.8, 1.0],
            "reasoning_keywords": ["local", "small", "no transfer"]
        }
    },
    
    "large_s3": {
        "description": "Large S3 files (>100MB) → Expect EMR execution",
        "command": "Merge paired-end reads from S3",
        "workflow_plan": {
            "description": "Read merging operation on large RNA-seq paired-end files stored in S3",
            "operations": [
                {
                    "operation_name": "read_merging",
                    "tool_name": "bbmerge"
                }
            ],
            "data_inputs": [
                {
                    "uri": "s3://helix-biodata/rnaseq/sample_R1.fastq",
                    "size_bytes": 262144000,  # 250MB
                    "location_type": "S3",
                    "metadata": {"read": "R1", "format": "fastq"}
                },
                {
                    "uri": "s3://helix-biodata/rnaseq/sample_R2.fastq",
                    "size_bytes": 251658240,  # 240MB
                    "location_type": "S3",
                    "metadata": {"read": "R2", "format": "fastq"}
                }
            ]
        },
        "expected": {
            "infrastructure": "EMR",
            "confidence_range": [0.8, 1.0],
            "reasoning_keywords": ["EMR", "S3", "transfer", "490MB", "100MB"]
        }
    },
    
    "unknown_sizes": {
        "description": "Unknown file sizes → Expect reduced confidence, assumptions",
        "command": "Analyze files with unknown sizes",
        "workflow_plan": {
            "description": "Quality assessment of files with unknown sizes (permissions error or non-existent)",
            "operations": [
                {
                    "operation_name": "quality_control",
                    "tool_name": "fastqc"
                }
            ],
            "data_inputs": [
                {
                    "uri": "s3://restricted-bucket/file1.fastq",
                    "size_bytes": None,  # Unknown
                    "location_type": "S3",
                    "metadata": {"note": "Size unavailable due to permissions", "format": "fastq"}
                },
                {
                    "uri": "s3://restricted-bucket/file2.fastq",
                    "size_bytes": None,  # Unknown
                    "location_type": "S3",
                    "metadata": {"note": "Size unavailable due to permissions", "format": "fastq"}
                }
            ]
        },
        "expected": {
            "infrastructure": "EMR",  # Assume large for S3 unknowns
            "confidence_range": [0.5, 0.7],
            "reasoning_keywords": ["unknown", "assume", "confidence"]
        }
    },
    
    "mixed_locations": {
        "description": "Mixed S3 + Local files → Expect EMR with upload recommendation",
        "command": "Process files from mixed locations",
        "workflow_plan": {
            "description": "Alignment workflow with mixed S3 and local input files",
            "operations": [
                {
                    "operation_name": "alignment",
                    "tool_name": "bwa"
                }
            ],
            "data_inputs": [
                {
                    "uri": "s3://helix-biodata/reference/hg38.fa",
                    "size_bytes": 314572800,  # 300MB
                    "location_type": "S3",
                    "metadata": {"type": "reference", "format": "fasta"}
                },
                {
                    "uri": "/local/data/sample.fastq",
                    "size_bytes": 52428800,  # 50MB
                    "location_type": "Local",
                    "metadata": {"type": "reads", "format": "fastq"}
                }
            ]
        },
        "expected": {
            "infrastructure": "EMR",
            "confidence_range": [0.7, 0.9],
            "reasoning_keywords": ["mixed", "upload", "S3 files dominate"]
        }
    },
    
    "gpu_required": {
        "description": "GPU-required tool → Expect EC2 with GPU or warning",
        "command": "Run AlphaFold protein structure prediction",
        "workflow_plan": {
            "description": "Protein structure prediction using AlphaFold (requires GPU)",
            "operations": [
                {
                    "operation_name": "structure_prediction",
                    "tool_name": "alphafold"
                }
            ],
            "data_inputs": [
                {
                    "uri": "/data/protein.fasta",
                    "size_bytes": 1000,  # 1KB (small protein sequence)
                    "location_type": "Local",
                    "metadata": {"type": "protein sequence", "format": "fasta"}
                }
            ]
        },
        "expected": {
            "infrastructure": "EC2",  # Or Batch with GPU
            "confidence_range": [0.6, 0.8],
            "reasoning_keywords": ["GPU", "AlphaFold"]
        }
    }
}


# ============================================================================
# CLI Demo Functions
# ============================================================================

def print_header(title: str):
    """Print formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80 + "\n")


def print_json(data: Dict[str, Any], title: str = "JSON Output"):
    """Pretty-print JSON data."""
    print(f"\n{title}:")
    print("-" * 80)
    print(json.dumps(data, indent=2))
    print("-" * 80)


def validate_expectations(result: Dict[str, Any], expected: Dict[str, Any]) -> bool:
    """Validate that results match expectations."""
    print("\n📊 Validation Results:")
    print("-" * 80)
    
    all_passed = True
    
    # Check infrastructure
    if "infrastructure" in expected:
        actual = result.get("infrastructure")
        expected_infra = expected["infrastructure"]
        passed = actual == expected_infra
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status} Infrastructure: expected={expected_infra}, actual={actual}")
        all_passed = all_passed and passed
    
    # Check confidence range
    if "confidence_range" in expected:
        actual = result.get("confidence_score", 0.0)
        min_conf, max_conf = expected["confidence_range"]
        passed = min_conf <= actual <= max_conf
        status = "✅ PASS" if passed else "❌ FAIL"
        print(f"{status} Confidence: expected={min_conf}-{max_conf}, actual={actual:.2f}")
        all_passed = all_passed and passed
    
    # Check reasoning keywords
    if "reasoning_keywords" in expected:
        reasoning = result.get("reasoning", "").lower()
        keywords = expected["reasoning_keywords"]
        found = [kw for kw in keywords if kw.lower() in reasoning]
        missing = [kw for kw in keywords if kw.lower() not in reasoning]
        
        if found:
            print(f"✅ PASS Reasoning keywords found: {', '.join(found)}")
        if missing:
            print(f"⚠️  WARN Reasoning keywords missing: {', '.join(missing)}")
            # Don't fail on missing keywords, just warn
    
    print("-" * 80)
    return all_passed


async def run_scenario(scenario_name: str, scenario_data: Dict[str, Any], mock_mode: bool = True):
    """Run a predefined scenario through the full pipeline."""
    
    print_header(f"Scenario: {scenario_name}")
    print(f"Description: {scenario_data['description']}")
    print(f"Command: {scenario_data['command']}")
    
    # Step 1: Create WorkflowPlan from scenario data
    print("\n📋 Step 1: Create WorkflowPlan")
    workflow_plan = WorkflowPlan(**scenario_data["workflow_plan"])
    print(f"  ✓ Created WorkflowPlan with {len(workflow_plan.data_inputs)} inputs, {len(workflow_plan.operations)} operations")
    
    # Step 2: (Optional) Demonstrate Phase 2 tools
    if not mock_mode:
        print("\n🔧 Step 2: Query Phase 2 Tools (FileMetadataInspector)")
        inspector = FileMetadataInspector(mock_mode=False)
        uris = [inp.uri for inp in workflow_plan.data_inputs]
        file_metadata = inspector.inspect_files(uris)
        for fm in file_metadata:
            print(f"  - {fm.uri}: {fm.size_mb:.2f} MB (source: {fm.source}, confidence: {fm.size_confidence})")
    else:
        print("\n🔧 Step 2: Skipping Phase 2 tools demo (mock_mode=True)")
    
    # Step 3: Initialize Orchestrator
    print("\n🎯 Step 3: Initialize Orchestrator")
    orchestrator = Orchestrator()
    print(f"  ✓ Orchestrator initialized (request_id generation enabled)")
    
    # Step 4: Run the full pipeline
    print("\n🚀 Step 4: Run Multi-Agent Pipeline")
    try:
        trace = await orchestrator.run_pipeline(
            command=scenario_data["command"],
            workflow_plan=workflow_plan,
            request_id=f"demo_{scenario_name}"
        )
        
        print(f"  ✓ Pipeline completed successfully")
        print(f"  ✓ Duration: {trace.duration_ms:.2f} ms")
        print(f"  ✓ Agents invoked: {len(trace.agent_invocations)}")
        
        # Step 5: Display InfraDecision
        print_header("Infrastructure Decision (Agent 1 Output)")
        infra_invocation = [inv for inv in trace.agent_invocations if inv.agent_name == "InfrastructureDecisionAgent"][0]
        infra_decision = infra_invocation.output
        
        print(f"🎯 Recommended Infrastructure: {infra_decision.infrastructure}")
        print(f"🎯 Confidence Score: {infra_decision.confidence_score:.2f}")
        print(f"🎯 Decision Summary: {infra_decision.decision_summary}")
        print(f"\n📊 File Analysis:")
        print(f"  - Total Size: {infra_decision.file_analysis.total_size_mb:.2f} MB")
        print(f"  - File Count: {infra_decision.file_analysis.file_count}")
        print(f"  - Unknown Sizes: {infra_decision.file_analysis.unknown_sizes}")
        print(f"  - All in S3: {infra_decision.file_analysis.all_in_s3}")
        
        print(f"\n💰 Cost Analysis:")
        cost_min, cost_max = infra_decision.cost_analysis.estimated_cost_range_usd
        print(f"  - Estimated Cost Range: ${cost_min:.2f} - ${cost_max:.2f}")
        print(f"  - Cost Class: {infra_decision.cost_analysis.cost_class}")
        print(f"  - Assumptions: {infra_decision.cost_analysis.cost_assumptions}")
        print(f"  - Cost Confidence: {infra_decision.cost_analysis.cost_confidence:.2f}")
        
        if infra_decision.warnings:
            print(f"\n⚠️  Warnings:")
            for warning in infra_decision.warnings:
                print(f"  - {warning}")
        
        if infra_decision.alternatives:
            print(f"\n🔄 Alternatives:")
            for alt in infra_decision.alternatives:
                print(f"  - {alt.infrastructure}: {alt.reasoning}")
                print(f"    Trade-offs: {alt.tradeoffs}")
        
        # Step 6: Display ExecutionToolSpec
        print_header("Execution Tool Spec (Agent 2 Output)")
        impl_invocation = [inv for inv in trace.agent_invocations if inv.agent_name == "ImplementationAgent"][0]
        execution_spec = impl_invocation.output
        
        print(f"🔧 Tool: {execution_spec.tool_name}")
        print(f"🔧 Infrastructure: {execution_spec.infrastructure}")
        print(f"🔧 Confidence Score: {execution_spec.confidence_score:.2f}")
        
        if execution_spec.container_spec:
            print(f"\n📦 Container:")
            print(f"  - Image: {execution_spec.container_spec.image}")
            print(f"  - Type: {execution_spec.container_spec.image_type}")
        else:
            print(f"\n📦 Container: Native execution (no container)")
        
        print(f"\n💻 Commands ({len(execution_spec.commands)}):")
        for i, cmd in enumerate(execution_spec.commands, 1):
            print(f"  {i}. {cmd.name}")
            print(f"     Command: {cmd.command[:80]}{'...' if len(cmd.command) > 80 else ''}")
            print(f"     Timeout: {cmd.timeout_minutes} min")
        
        print(f"\n💾 Resources:")
        print(f"  - CPU: {execution_spec.resource_requirements.min_cpu_cores} cores")
        print(f"  - Memory: {execution_spec.resource_requirements.min_memory_gb} GB")
        print(f"  - Disk: {execution_spec.resource_requirements.min_disk_gb} GB")
        print(f"  - GPU: {'Yes' if execution_spec.resource_requirements.gpu_required else 'No'}")
        
        print(f"\n🔁 Retry Policy:")
        print(f"  - Max Retries: {execution_spec.retry_policy.max_retries}")
        print(f"  - Retry On: {', '.join(execution_spec.retry_policy.retry_on)}")
        
        if execution_spec.warnings:
            print(f"\n⚠️  Warnings:")
            for warning in execution_spec.warnings:
                print(f"  - {warning}")
        
        # Step 7: Validate against expectations (if provided)
        if "expected" in scenario_data:
            print_header("Expectation Validation")
            infra_dict = infra_decision.model_dump()
            validation_passed = validate_expectations(infra_dict, scenario_data["expected"])
            
            if validation_passed:
                print("\n✅ All validations PASSED")
            else:
                print("\n⚠️  Some validations FAILED (see details above)")
        
        return trace
        
    except Exception as e:
        print(f"\n❌ Pipeline failed: {e}")
        logger.exception("Pipeline execution failed")
        return None


async def run_custom_scenario(command: str, files: list[str], mock_mode: bool = True):
    """Run a custom user-defined scenario."""
    
    print_header("Custom Scenario")
    print(f"Command: {command}")
    print(f"Files: {', '.join(files)}")
    
    # Create a basic WorkflowPlan
    data_inputs = []
    for uri in files:
        # Infer format from extension
        ext = Path(uri).suffix.lower().lstrip('.')
        format_map = {
            'fq': 'fastq', 'fastq': 'fastq',
            'fa': 'fasta', 'fasta': 'fasta',
            'bam': 'bam', 'sam': 'sam',
            'vcf': 'vcf', 'bed': 'bed'
        }
        file_format = format_map.get(ext, 'unknown')
        
        # Infer location type
        if uri.startswith('s3://'):
            location_type = "S3"
        elif uri.startswith('/') or uri.startswith('./'):
            location_type = "Local"
        else:
            location_type = "Unknown"
        
        data_inputs.append({
            "uri": uri,
            "size_bytes": None,  # Unknown for custom scenarios
            "location_type": location_type,
            "metadata": {"format": file_format}
        })
    
    workflow_plan_data = {
        "description": command,
        "operations": [
            {
                "operation_name": "custom_operation",
                    "tool_name": "unknown",
                "description": command
            }
        ],
        "data_inputs": data_inputs,
        "expected_outputs": []
    }
    
    workflow_plan = WorkflowPlan(**workflow_plan_data)
    
    # Run through orchestrator
    orchestrator = Orchestrator()
    trace = await orchestrator.run_pipeline(
        command=command,
        workflow_plan=workflow_plan,
        request_id="demo_custom"
    )
    
    if trace and trace.success:
        print("\n✅ Custom scenario completed successfully")
        # Display summary
        infra_inv = [inv for inv in trace.agent_invocations if inv.agent_name == "InfrastructureDecisionAgent"][0]
        impl_inv = [inv for inv in trace.agent_invocations if inv.agent_name == "ImplementationAgent"][0]
        
        print(f"\n📊 Summary:")
        print(f"  Infrastructure: {infra_inv.output.infrastructure}")
        print(f"  Confidence: {infra_inv.output.confidence_score:.2f}")
        print(f"  Tool: {impl_inv.output.tool_name}")
        print(f"  Duration: {trace.duration_ms:.2f} ms")
    else:
        print("\n❌ Custom scenario failed")
    
    return trace


# ============================================================================
# Main CLI
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="CLI Demo for Multi-Agent Refactored Architecture (Phases 1-4)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run predefined scenarios
  python backend/cli_demo.py --scenario small_local
  python backend/cli_demo.py --scenario large_s3
  python backend/cli_demo.py --scenario unknown_sizes
  python backend/cli_demo.py --all
  
  # Custom scenario
  python backend/cli_demo.py --custom --command "Run FastQC" --files "s3://bucket/file.fq"
  
Available scenarios: small_local, large_s3, unknown_sizes, mixed_locations, gpu_required
        """
    )
    
    parser.add_argument(
        "--scenario",
        choices=list(SCENARIOS.keys()),
        help="Run a predefined scenario"
    )
    
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all predefined scenarios"
    )
    
    parser.add_argument(
        "--custom",
        action="store_true",
        help="Run a custom scenario (requires --command and --files)"
    )
    
    parser.add_argument(
        "--command",
        help="Custom command description (for --custom)"
    )
    
    parser.add_argument(
        "--files",
        help="Comma-separated file URIs (for --custom)"
    )
    
    parser.add_argument(
        "--no-mock",
        action="store_true",
        help="Disable mock mode for Phase 2 tools (use real S3/file access)"
    )
    
    parser.add_argument(
        "--json",
        action="store_true",
        help="Output results as JSON"
    )
    
    args = parser.parse_args()
    
    # Validate arguments
    if not (args.scenario or args.all or args.custom):
        parser.error("Must specify --scenario, --all, or --custom")
    
    if args.custom and (not args.command or not args.files):
        parser.error("--custom requires both --command and --files")
    
    mock_mode = not args.no_mock
    
    # Run scenarios
    if args.all:
        print_header("Running All Predefined Scenarios")
        for scenario_name in SCENARIOS.keys():
            asyncio.run(run_scenario(scenario_name, SCENARIOS[scenario_name], mock_mode))
            print("\n" + "=" * 80 + "\n")
    
    elif args.scenario:
        asyncio.run(run_scenario(args.scenario, SCENARIOS[args.scenario], mock_mode))
    
    elif args.custom:
        files = [f.strip() for f in args.files.split(",")]
        asyncio.run(run_custom_scenario(args.command, files, mock_mode))


if __name__ == "__main__":
    main()
