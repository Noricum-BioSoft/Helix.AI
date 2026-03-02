"""
Scenario reporter - formats and displays results.

Provides multiple output formats:
- Human-readable console output
- JSON reports for CI/CD
- Diff views for baseline comparison
"""

from typing import List, Dict, Any, Optional
from pathlib import Path
import json
from datetime import datetime

from .executor import ExecutionTrace
from .validator import ValidationResult


class ScenarioReporter:
    """
    Reports on scenario execution results.
    
    Supports multiple output formats and diff views for regression detection.
    """
    
    def __init__(self, verbose: bool = False):
        self.verbose = verbose
    
    def report_single(self, trace: ExecutionTrace, show_contracts: bool = False) -> str:
        """
        Generate a report for a single scenario execution.
        
        Args:
            trace: Execution trace to report on
            show_contracts: Whether to include full contract details
        
        Returns:
            Formatted report string
        """
        lines = []
        lines.append(f"\n{'='*70}")
        lines.append(f"Scenario: {trace.scenario_id}")
        lines.append(f"{'='*70}")
        
        # Metadata
        lines.append(f"\nCategory: {trace.metadata.get('category', 'unknown')}")
        if trace.metadata.get('tags'):
            lines.append(f"Tags: {', '.join(trace.metadata['tags'])}")
        lines.append(f"Timestamp: {trace.timestamp}")
        lines.append(f"Duration: {trace.duration_ms}ms")
        
        # Agent sequence
        lines.append(f"\nAgent Sequence:")
        agent_flow = " → ".join(trace.agent_sequence)
        lines.append(f"  {agent_flow}")
        
        # Contracts
        lines.append(f"\nContracts Produced: {len(trace.contracts)}")
        if show_contracts:
            for i, contract in enumerate(trace.contracts, 1):
                lines.append(f"  {i}. {contract['agent']} → {contract['contract_type']}")
        
        # Validation results
        if trace.validation_result:
            validation = ValidationResult(**trace.validation_result)
            lines.append(f"\nValidation: {validation.summary}")
            
            if validation.get_errors():
                lines.append(f"\n  ❌ Errors ({len(validation.get_errors())}):")
                for error in validation.get_errors():
                    lines.append(f"    • [{error.category}] {error.message}")
            
            if validation.get_warnings():
                lines.append(f"\n  ⚠  Warnings ({len(validation.get_warnings())}):")
                for warning in validation.get_warnings():
                    lines.append(f"    • [{warning.category}] {warning.message}")
        
        # Errors
        if trace.errors:
            lines.append(f"\n❌ Execution Errors:")
            for error in trace.errors:
                lines.append(f"  • {error}")
        
        lines.append(f"{'='*70}\n")
        
        return "\n".join(lines)
    
    def report_batch(
        self,
        traces: List[ExecutionTrace],
        group_by: Optional[str] = None
    ) -> str:
        """
        Generate a summary report for multiple scenarios.
        
        Args:
            traces: List of execution traces
            group_by: Optional grouping key ('category', 'status')
        
        Returns:
            Formatted summary report
        """
        lines = []
        lines.append(f"\n{'='*70}")
        lines.append(f"Scenario Execution Summary")
        lines.append(f"{'='*70}")
        lines.append(f"\nTotal Scenarios: {len(traces)}")
        
        # Count pass/fail
        passed = sum(
            1 for t in traces
            if t.validation_result and ValidationResult(**t.validation_result).passed
        )
        failed = len(traces) - passed
        
        lines.append(f"Passed: {passed} ✓")
        lines.append(f"Failed: {failed} ✗")
        
        # Average duration
        avg_duration = sum(t.duration_ms for t in traces) / len(traces) if traces else 0
        lines.append(f"Average Duration: {avg_duration:.0f}ms")
        
        # Group by category if requested
        if group_by == "category":
            lines.append(f"\nBy Category:")
            by_category: Dict[str, List[ExecutionTrace]] = {}
            for trace in traces:
                cat = trace.metadata.get('category', 'unknown')
                by_category.setdefault(cat, []).append(trace)
            
            for category, cat_traces in sorted(by_category.items()):
                cat_passed = sum(
                    1 for t in cat_traces
                    if t.validation_result and ValidationResult(**t.validation_result).passed
                )
                lines.append(f"  {category}: {cat_passed}/{len(cat_traces)} passed")
        
        # List failures
        if failed > 0:
            lines.append(f"\nFailed Scenarios:")
            for trace in traces:
                if trace.validation_result:
                    validation = ValidationResult(**trace.validation_result)
                    if not validation.passed:
                        lines.append(f"  ✗ {trace.scenario_id}")
                        lines.append(f"    {validation.summary}")
        
        lines.append(f"{'='*70}\n")
        
        return "\n".join(lines)
    
    def report_diff(
        self,
        current: ExecutionTrace,
        baseline: ExecutionTrace
    ) -> str:
        """
        Generate a diff report comparing current and baseline traces.
        
        Args:
            current: Current execution trace
            baseline: Baseline execution trace
        
        Returns:
            Formatted diff report
        """
        lines = []
        lines.append(f"\n{'='*70}")
        lines.append(f"Diff Report: {current.scenario_id}")
        lines.append(f"{'='*70}")
        
        # Agent sequence diff
        if current.agent_sequence != baseline.agent_sequence:
            lines.append(f"\n❌ Agent Sequence Changed:")
            lines.append(f"  Baseline: {' → '.join(baseline.agent_sequence)}")
            lines.append(f"  Current:  {' → '.join(current.agent_sequence)}")
        else:
            lines.append(f"\n✓ Agent Sequence: Unchanged")
        
        # Contract count diff
        baseline_contract_count = len(baseline.contracts)
        current_contract_count = len(current.contracts)
        
        if current_contract_count != baseline_contract_count:
            lines.append(f"\n❌ Contract Count Changed:")
            lines.append(f"  Baseline: {baseline_contract_count}")
            lines.append(f"  Current:  {current_contract_count}")
        else:
            lines.append(f"\n✓ Contract Count: {current_contract_count} (unchanged)")
        
        # Duration diff
        duration_delta = current.duration_ms - baseline.duration_ms
        duration_pct = (duration_delta / baseline.duration_ms * 100) if baseline.duration_ms else 0
        
        lines.append(f"\nDuration:")
        lines.append(f"  Baseline: {baseline.duration_ms}ms")
        lines.append(f"  Current:  {current.duration_ms}ms")
        lines.append(f"  Delta:    {duration_delta:+d}ms ({duration_pct:+.1f}%)")
        
        # Validation status diff
        baseline_passed = (
            baseline.validation_result
            and ValidationResult(**baseline.validation_result).passed
        )
        current_passed = (
            current.validation_result
            and ValidationResult(**current.validation_result).passed
        )
        
        if baseline_passed != current_passed:
            lines.append(f"\n❌ Validation Status Changed:")
            lines.append(f"  Baseline: {'PASSED' if baseline_passed else 'FAILED'}")
            lines.append(f"  Current:  {'PASSED' if current_passed else 'FAILED'}")
        else:
            status = "PASSED" if current_passed else "FAILED"
            lines.append(f"\n✓ Validation Status: {status} (unchanged)")
        
        lines.append(f"{'='*70}\n")
        
        return "\n".join(lines)
    
    def export_json(
        self,
        traces: List[ExecutionTrace],
        output_path: Path
    ) -> None:
        """
        Export traces as JSON for CI/CD integration.
        
        Args:
            traces: Execution traces to export
            output_path: Path to write JSON report
        """
        report = {
            "timestamp": datetime.utcnow().isoformat(),
            "total_scenarios": len(traces),
            "passed": sum(
                1 for t in traces
                if t.validation_result and ValidationResult(**t.validation_result).passed
            ),
            "failed": sum(
                1 for t in traces
                if not t.validation_result or not ValidationResult(**t.validation_result).passed
            ),
            "traces": [t.model_dump() for t in traces]
        }
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2)
    
    def export_junit(
        self,
        traces: List[ExecutionTrace],
        output_path: Path
    ) -> None:
        """
        Export traces as JUnit XML for CI/CD integration.
        
        Args:
            traces: Execution traces to export
            output_path: Path to write JUnit XML
        """
        # Simple JUnit XML generation
        lines = []
        lines.append('<?xml version="1.0" encoding="UTF-8"?>')
        
        total = len(traces)
        failures = sum(
            1 for t in traces
            if not t.validation_result or not ValidationResult(**t.validation_result).passed
        )
        
        lines.append(f'<testsuite name="Helix.AI Scenarios" tests="{total}" failures="{failures}">')
        
        for trace in traces:
            passed = (
                trace.validation_result
                and ValidationResult(**trace.validation_result).passed
            )
            
            lines.append(f'  <testcase name="{trace.scenario_id}" time="{trace.duration_ms / 1000:.3f}">')
            
            if not passed:
                if trace.validation_result:
                    validation = ValidationResult(**trace.validation_result)
                    message = validation.summary
                    errors = "\n".join(e.message for e in validation.get_errors())
                else:
                    message = "Execution failed"
                    errors = "\n".join(trace.errors)
                
                lines.append(f'    <failure message="{message}">')
                lines.append(f'      {errors}')
                lines.append(f'    </failure>')
            
            lines.append(f'  </testcase>')
        
        lines.append('</testsuite>')
        
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write('\n'.join(lines))
