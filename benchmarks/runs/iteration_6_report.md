# Iteration Report (iteration_6_report)

- Generated: 2026-03-15T20:17:58.724964+00:00
- Session: 8e89c58c-d2b1-409d-ad55-536af858506e
- Scenario: bulk-rnaseq-e2e
- Raw Score: 29 / 40
- Overall Percentage: 72.5%
- Threshold Met: No
- Critical failures: 3

## Lessons Learned
- 1 critical failure(s): expectation=rerun_expected; tool_alignment=aligned; execution_vs_planning=uncertain; artifact_resolution=resolved_or_not_required; output_usefulness=specific; deductions=hard_execution_failure
- 1 critical failure(s): expectation=visualization_update_expected; tool_alignment=aligned; execution_vs_planning=uncertain; artifact_resolution=not_primary; output_usefulness=specific; deductions=hard_execution_failure
- 1 critical failure(s): expectation=historical_recreation_expected; tool_alignment=misaligned; execution_vs_planning=incorrect_planning_instead_of_execution; artifact_resolution=resolved_or_not_required; output_usefulness=weak; deductions=output_not_specific,plan_instead_of_execution,wrong_tool_class

## Remediation Targets
- Harden execution path robustness and fallback behavior.

## Critical Failures
- `turn_09` tool=`patch_and_rerun` status=`error`: expectation=rerun_expected; tool_alignment=aligned; execution_vs_planning=uncertain; artifact_resolution=resolved_or_not_required; output_usefulness=specific; deductions=hard_execution_failure
- `turn_10` tool=`patch_and_rerun` status=`error`: expectation=visualization_update_expected; tool_alignment=aligned; execution_vs_planning=uncertain; artifact_resolution=not_primary; output_usefulness=specific; deductions=hard_execution_failure
- `turn_19` tool=`__plan__` status=`workflow_planned`: expectation=historical_recreation_expected; tool_alignment=misaligned; execution_vs_planning=incorrect_planning_instead_of_execution; artifact_resolution=resolved_or_not_required; output_usefulness=weak; deductions=output_not_specific,plan_instead_of_execution,wrong_tool_class
