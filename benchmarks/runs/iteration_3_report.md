# Iteration Report (iteration_3_report)

- Generated: 2026-03-15T20:08:14.546151+00:00
- Session: 983bcde1-8621-4e1e-9ebc-1c09424347a8
- Scenario: bulk-rnaseq-e2e
- Raw Score: 31 / 40
- Overall Percentage: 77.5%
- Threshold Met: No
- Critical failures: 1

## Lessons Learned
- 1 critical failure(s): expectation=historical_recreation_expected; tool_alignment=misaligned; execution_vs_planning=incorrect_planning_instead_of_execution; artifact_resolution=resolved_or_not_required; output_usefulness=weak; deductions=output_not_specific,plan_instead_of_execution,wrong_tool_class

## Remediation Targets
- Harden execution path robustness and fallback behavior.

## Critical Failures
- `turn_19` tool=`__plan__` status=`workflow_planned`: expectation=historical_recreation_expected; tool_alignment=misaligned; execution_vs_planning=incorrect_planning_instead_of_execution; artifact_resolution=resolved_or_not_required; output_usefulness=weak; deductions=output_not_specific,plan_instead_of_execution,wrong_tool_class
