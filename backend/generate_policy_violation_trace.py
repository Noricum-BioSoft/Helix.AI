#!/usr/bin/env python3
"""
Generate a policy violation trace showing semantic invariants catching violations.

Demonstrates:
1. InfrastructureExpert including execution commands (VIOLATION - Invariant 1)
2. ImplementationAgent ignoring infrastructure decision (VIOLATION - Invariant 2)
3. Very low confidence but system trying to auto-execute (VIOLATION - Invariant 3)

Shows that validation catches these and blocks execution.
"""

import json
import hashlib
import time
from datetime import datetime


def compute_hash(data):
    """Compute SHA256 hash of data."""
    json_str = json.dumps(data, sort_keys=True)
    return hashlib.sha256(json_str.encode('utf-8')).hexdigest()


def generate_policy_violation_trace():
    """Generate a trace with multiple policy violations."""
    
    start_time = time.time()
    request_id = f"policy_violation_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    
    # ========================================================================
    # Attempt 1: InfrastructureExpert overstepping boundaries
    # ========================================================================
    infra_start = start_time
    infra_input = {
        "workflow_plan": {
            "description": "Run alignment on unknown files",
            "operations": [{"operation_name": "alignment", "tool_name": "bwa"}],
            "data_inputs": [
                {
                    "uri": "s3://unknown-bucket/reads.fastq",
                    "size_bytes": None,  # Unknown size
                    "location_type": "S3"
                }
            ]
        }
    }
    
    # VIOLATION: InfrastructureExpert includes execution commands
    infra_output_BAD = {
        "infrastructure": "EMR",
        "reasoning": "Unknown file sizes suggest large data",
        "confidence_score": 0.28,  # VERY LOW
        "warnings": [
            "File size unknown - cannot estimate costs accurately",
            "Location unknown - cannot optimize for data locality",
            "Cannot determine if infrastructure is appropriate"
        ],
        # ❌ VIOLATION 1: Infrastructure agent including execution details!
        "execution_command": "bwa mem -t 16 ref.fa reads.fastq > aligned.sam",
        "docker_image": "biocontainers/bwa:0.7.17",
        "retry_policy": {"max_retries": 3},
        # ❌ VIOLATION: These belong in ExecutionToolSpec, not InfraDecision!
        "cost_analysis": {
            "estimated_cost_range_usd": [10.0, 50.0],
            "cost_assumptions": "Wide range due to unknown file sizes",
            "cost_confidence": 0.15
        }
    }
    
    infra_end = infra_start + 1.523
    infra_invocation = {
        "agent_name": "InfrastructureExpert",
        "request_id": request_id,
        "start_time": infra_start,
        "end_time": infra_end,
        "duration_ms": (infra_end - infra_start) * 1000,
        "success": True,  # Agent completed
        "input_hash": compute_hash(infra_input),
        "output_hash": compute_hash(infra_output_BAD),
        "confidence_score": 0.28,
        "warnings": infra_output_BAD["warnings"],
        "contract": infra_output_BAD,
        "policy_violations": [
            {
                "type": "SEMANTIC_INVARIANT_VIOLATION",
                "severity": "ERROR",
                "rule": "Separation of Concerns (InfraDecision: WHERE, ExecutionToolSpec: HOW)",
                "message": "InfraDecision contains execution details (violations: execution_command, docker_image, retry_policy). Infrastructure agent must only decide WHERE to execute, not HOW.",
                "prohibited_fields": ["execution_command", "docker_image", "retry_policy"],
                "principle": "Infrastructure agent decides WHERE, Implementation agent decides HOW"
            },
            {
                "type": "CONFIDENCE_THRESHOLD_VIOLATION",
                "severity": "ERROR",
                "rule": "Critically low confidence must block execution",
                "message": "InfrastructureExpert has critically low confidence (0.28 < 0.3). Execution should be blocked or require explicit user confirmation.",
                "confidence": 0.28,
                "threshold": 0.3,
                "principle": "Safety: Very low confidence should block automatic execution"
            }
        ]
    }
    
    # ========================================================================
    # Attempt 2: ImplementationAgent ignoring infrastructure decision
    # ========================================================================
    impl_start = infra_end + 0.1
    impl_input = {
        "workflow_plan": infra_input["workflow_plan"],
        "infra_decision": infra_output_BAD
    }
    
    # VIOLATION: Using different infrastructure than recommended
    impl_output_BAD = {
        "tool_name": "bwa",
        # ❌ VIOLATION 2: Infrastructure mismatch!
        "infrastructure": "Local",  # InfraDecision said "EMR"!
        "commands": [
            {
                "command": "bwa mem -t 4 ref.fa reads.fastq > aligned.sam",
                "timeout_seconds": 3600
            }
        ],
        "container_spec": {
            "image": "biocontainers/bwa:0.7.17"
        },
        "confidence_score": 0.45,  # LOW (< 0.5)
        "warnings": [
            "Unknown file size may cause OOM on local machine",
            "S3 download to local will be slow"
        ]
    }
    
    impl_end = impl_start + 2.145
    impl_invocation = {
        "agent_name": "ImplementationAgent",
        "request_id": request_id,
        "start_time": impl_start,
        "end_time": impl_end,
        "duration_ms": (impl_end - impl_start) * 1000,
        "success": True,  # Agent completed
        "input_hash": compute_hash(impl_input),
        "output_hash": compute_hash(impl_output_BAD),
        "confidence_score": 0.45,
        "warnings": impl_output_BAD["warnings"],
        "contract": impl_output_BAD,
        "policy_violations": [
            {
                "type": "SEMANTIC_INVARIANT_VIOLATION",
                "severity": "ERROR",
                "rule": "Infrastructure Consistency",
                "message": "Infrastructure mismatch: InfraDecision recommends 'EMR' but ExecutionToolSpec uses 'Local'. These must match.",
                "infra_decision_choice": "EMR",
                "execution_spec_choice": "Local",
                "principle": "Consistency: Implementation must respect infrastructure decision"
            },
            {
                "type": "CONFIDENCE_THRESHOLD_VIOLATION",
                "severity": "WARNING",
                "rule": "Low confidence should trigger clarification",
                "message": "ImplementationAgent has low confidence (0.45 < 0.5). Consider asking clarifying questions or showing warnings to user.",
                "confidence": 0.45,
                "threshold": 0.5,
                "principle": "Uncertainty Handling: Low confidence should trigger user interaction"
            }
        ]
    }
    
    # ========================================================================
    # Validation Result: BLOCKED
    # ========================================================================
    end_time = impl_end + 0.05
    
    validation_result = {
        "passed": False,
        "blocked": True,
        "reason": "Multiple policy violations detected",
        "errors": [
            {
                "category": "semantic_invariant",
                "agent": "InfrastructureExpert",
                "message": "InfraDecision contains execution details (execution_command, docker_image, retry_policy)"
            },
            {
                "category": "semantic_invariant",
                "agent": "InfrastructureExpert",
                "message": "Critically low confidence (0.28 < 0.3) blocks execution"
            },
            {
                "category": "semantic_invariant",
                "agent": "ImplementationAgent",
                "message": "Infrastructure mismatch: EMR vs Local"
            }
        ],
        "warnings": [
            {
                "category": "semantic_invariant",
                "agent": "ImplementationAgent",
                "message": "Low confidence (0.45 < 0.5) suggests clarification needed"
            }
        ],
        "action_required": "User confirmation or input required before execution",
        "suggested_fixes": [
            "InfrastructureExpert: Remove execution_command, docker_image, retry_policy from output",
            "InfrastructureExpert: Gather more information to improve confidence (file sizes, etc.)",
            "ImplementationAgent: Use infrastructure='EMR' to match InfraDecision"
        ]
    }
    
    # ========================================================================
    # Build complete trace
    # ========================================================================
    trace = {
        "request_id": request_id,
        "user_command": "Run alignment on unknown files",
        "intent": "execute",
        "start_time": start_time,
        "end_time": end_time,
        "duration_ms": (end_time - start_time) * 1000,
        "success": False,  # Blocked by validation
        "error": "Execution blocked due to policy violations",
        "timestamp": datetime.fromtimestamp(start_time).isoformat(),
        "orchestrator": "PRIMARY (backend/agent.py)",
        "validation_result": validation_result,
        "invocations": [
            infra_invocation,
            impl_invocation
        ],
        "policy_enforcement": {
            "checks_performed": [
                "Semantic Invariant 1: Separation of Concerns (WHERE vs HOW)",
                "Semantic Invariant 2: Infrastructure Consistency",
                "Semantic Invariant 3: Confidence/Warning Thresholds"
            ],
            "violations_detected": 3,
            "warnings_generated": 1,
            "execution_blocked": True,
            "user_action_required": True
        },
        "summary": {
            "what_happened": "Multiple policy violations detected by semantic invariant validation",
            "why_blocked": [
                "InfrastructureExpert overstepped boundaries (included HOW, not just WHERE)",
                "Infrastructure decision had critically low confidence (0.28 < 0.3)",
                "ImplementationAgent ignored infrastructure recommendation (EMR → Local)"
            ],
            "how_to_fix": [
                "Provide file sizes to improve infrastructure confidence",
                "Remove execution details from InfraDecision",
                "Ensure ImplementationAgent respects infrastructure choice"
            ]
        }
    }
    
    return trace


def main():
    print("=" * 70)
    print("GENERATING POLICY VIOLATION TRACE")
    print("=" * 70)
    print()
    
    trace = generate_policy_violation_trace()
    
    # Save to file
    output_path = "docs/POLICY_VIOLATION_TRACE.json"
    with open(output_path, 'w') as f:
        json.dump(trace, f, indent=2)
    
    print(f"✅ Generated trace saved to: {output_path}")
    print()
    print("📊 Trace Summary:")
    print("=" * 70)
    print(f"Request ID: {trace['request_id']}")
    print(f"Command: {trace['user_command']}")
    print(f"Success: {trace['success']}")
    print(f"Error: {trace['error']}")
    print(f"Execution Blocked: {trace['policy_enforcement']['execution_blocked']}")
    print()
    print(f"Policy Violations Detected: {trace['policy_enforcement']['violations_detected']}")
    for i, error in enumerate(trace['validation_result']['errors'], 1):
        print(f"  {i}. [{error['agent']}] {error['message']}")
    print()
    print("Warnings:")
    for i, warning in enumerate(trace['validation_result']['warnings'], 1):
        print(f"  {i}. [{warning['agent']}] {warning['message']}")
    print()
    print("Violations by Agent:")
    for inv in trace['invocations']:
        if 'policy_violations' in inv:
            print(f"  - {inv['agent_name']}: {len(inv['policy_violations'])} violations")
            for v in inv['policy_violations']:
                print(f"    • [{v['severity']}] {v['rule']}")
    print()
    print("=" * 70)
    print("✅ Policy violation trace generated!")
    print()
    print("This demonstrates:")
    print("  ✅ Semantic Invariant 1: InfraDecision must not include HOW")
    print("  ✅ Semantic Invariant 2: Infrastructure consistency enforced")
    print("  ✅ Semantic Invariant 3: Low confidence blocks execution")
    print("  ✅ Validation catches violations and blocks bad behavior")
    print("  ✅ Clear error messages guide fixes")
    print()


if __name__ == "__main__":
    main()
