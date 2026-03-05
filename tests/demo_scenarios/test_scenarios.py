"""
Pytest integration for scenario-driven agent testing.

**SYSTEM TESTED:** PRIMARY Orchestrator (backend/agent.py) 🟢

This module validates the PRODUCTION multi-agent system used by the HTTP API.
Tests assert agent names from the PRIMARY system (BioinformaticsExecutor, 
InfrastructureExpert, etc.).

For tests of the EXPERIMENTAL orchestrator (backend/orchestrator.py), 
see tests/unit/backend/test_phase3_agents.py.

See docs/ORCHESTRATION_DUALITY.md for details on the two orchestration systems.

---

This module provides pytest-based regression tests for the multi-agent system.
Each scenario is executed and validated automatically.

Usage:
    # Run all scenarios
    pytest tests/demo_scenarios/test_scenarios.py -v
    
    # Run specific category
    pytest tests/demo_scenarios/test_scenarios.py -k "ask"
    
    # Run with detailed traces
    pytest tests/demo_scenarios/test_scenarios.py -v --show-traces
    
    # Update baselines
    pytest tests/demo_scenarios/test_scenarios.py --update-baseline
"""

import pytest
import asyncio
from pathlib import Path
from typing import List

from framework import (
    Scenario,
    ScenarioLoader,
    ScenarioExecutor,
    ScenarioReporter,
    ExecutionTrace,
)


# Configuration
SCENARIOS_DIR = Path(__file__).parent / "scenarios"
BASELINES_DIR = Path(__file__).parent / "baselines"


# Note: pytest options are now defined in conftest.py


@pytest.fixture(scope="session")
def scenario_loader():
    """Create a scenario loader."""
    return ScenarioLoader(SCENARIOS_DIR)


@pytest.fixture
def scenario_executor(request):
    """Create a scenario executor (per-test scope to allow dynamic mock_mode)."""
    show_traces = request.config.getoption("--show-traces")
    return ScenarioExecutor(
        mock_mode=True,  # Default, can be overridden per scenario
        baseline_dir=BASELINES_DIR,
        verbose=show_traces
    )


@pytest.fixture(scope="session")
def scenario_reporter(request):
    """Create a scenario reporter."""
    show_traces = request.config.getoption("--show-traces")
    return ScenarioReporter(verbose=show_traces)


@pytest.fixture(scope="session")
def all_scenarios(scenario_loader):
    """Load all scenarios."""
    return scenario_loader.load_all_scenarios()


# Parametrize tests by scenario
def pytest_generate_tests(metafunc):
    """Generate test cases for each scenario."""
    if "scenario" in metafunc.fixturenames:
        loader = ScenarioLoader(SCENARIOS_DIR)
        scenarios = loader.load_all_scenarios()
        
        if not scenarios:
            pytest.skip("No scenarios found")
        
        metafunc.parametrize(
            "scenario",
            scenarios,
            ids=[s.metadata.id for s in scenarios]
        )


@pytest.mark.asyncio
async def test_scenario(
    scenario: Scenario,
    scenario_reporter: ScenarioReporter,
    request
):
    """
    Test a single scenario.
    
    This is the main test function that:
    1. Executes the scenario
    2. Validates the results
    3. Optionally compares with baseline
    4. Optionally updates baseline
    """
    # Create executor with scenario-specific mock_mode
    show_traces = request.config.getoption("--show-traces")
    mock_mode = scenario.metadata.requires_mock if scenario.metadata.requires_mock is not None else True
    
    # Skip real execution in regular test runs (too slow, needs backend + LLM + network)
    if not mock_mode:
        pytest.skip(f"Real execution scenario - run separately with: python test_real_execution.py --scenario {scenario.metadata.id}")
    
    scenario_executor = ScenarioExecutor(
        mock_mode=mock_mode,
        baseline_dir=BASELINES_DIR,
        verbose=show_traces
    )
    
    # Execute the scenario
    trace = await scenario_executor.execute(scenario)
    
    # Handle baseline operations
    update_baseline = request.config.getoption("--update-baseline")
    compare_baseline = request.config.getoption("--compare-baseline")
    
    if update_baseline:
        scenario_executor.save_baseline(trace)
        pytest.skip(f"Updated baseline for {scenario.metadata.id}")
    
    if compare_baseline:
        baseline = scenario_executor.load_baseline(scenario.metadata.id)
        if baseline:
            diff = scenario_executor.compare_with_baseline(trace, baseline)
            if diff["agent_sequence_changed"] or diff["contract_count_changed"]:
                diff_report = scenario_reporter.report_diff(trace, baseline)
                pytest.fail(f"Scenario diverged from baseline:\n{diff_report}")
    
    # Check for execution errors
    if trace.errors:
        pytest.fail(f"Execution errors: {', '.join(trace.errors)}")
    
    # Validate the trace
    if not trace.validation_result:
        pytest.fail("No validation result produced")
    
    from framework.validator import ValidationResult
    validation = ValidationResult(**trace.validation_result)
    
    # Report errors
    if not validation.passed:
        error_details = "\n".join(
            f"  • [{e.category}] {e.message}"
            for e in validation.get_errors()
        )
        pytest.fail(
            f"Scenario validation failed:\n"
            f"{validation.summary}\n"
            f"Errors:\n{error_details}"
        )
    
    # Warnings don't fail the test, but we can log them
    if validation.get_warnings():
        for warning in validation.get_warnings():
            print(f"⚠ Warning: [{warning.category}] {warning.message}")


# Category-specific test collections
class TestAskScenarios:
    """Test collection for ask/Q&A scenarios."""
    
    @pytest.mark.asyncio
    async def test_ask_scenarios(
        self,
        scenario_loader: ScenarioLoader,
        scenario_executor: ScenarioExecutor
    ):
        """Test all ask scenarios."""
        from framework.scenario import ScenarioCategory
        scenarios = scenario_loader.load_by_category(ScenarioCategory.ASK)
        
        results = []
        for scenario in scenarios:
            trace = await scenario_executor.execute(scenario)
            results.append((scenario, trace))
        
        # Assert all passed
        for scenario, trace in results:
            if trace.validation_result:
                from framework.validator import ValidationResult
                validation = ValidationResult(**trace.validation_result)
                assert validation.passed, f"Scenario {scenario.metadata.id} failed"


class TestExecuteScenarios:
    """Test collection for execute/workflow scenarios."""
    
    @pytest.mark.asyncio
    async def test_execute_scenarios(
        self,
        scenario_loader: ScenarioLoader,
        request
    ):
        """Test all execute scenarios."""
        from framework.scenario import ScenarioCategory
        scenarios = scenario_loader.load_by_category(ScenarioCategory.EXECUTE)
        
        results = []
        for scenario in scenarios:
            # Create executor with scenario-specific mock_mode
            show_traces = request.config.getoption("--show-traces")
            mock_mode = scenario.metadata.requires_mock if scenario.metadata.requires_mock is not None else True
            executor = ScenarioExecutor(mock_mode=mock_mode, baseline_dir=BASELINES_DIR, verbose=show_traces)
            
            trace = await executor.execute(scenario)
            results.append((scenario, trace))
        
        # Assert all passed
        for scenario, trace in results:
            if trace.validation_result:
                from framework.validator import ValidationResult
                validation = ValidationResult(**trace.validation_result)
                assert validation.passed, f"Scenario {scenario.metadata.id} failed"


class TestEdgeCases:
    """Test collection for edge case scenarios."""
    
    @pytest.mark.asyncio
    async def test_edge_cases(
        self,
        scenario_loader: ScenarioLoader,
        scenario_executor: ScenarioExecutor
    ):
        """Test all edge case scenarios."""
        from framework.scenario import ScenarioCategory
        scenarios = scenario_loader.load_by_category(ScenarioCategory.EDGE)
        
        results = []
        for scenario in scenarios:
            trace = await scenario_executor.execute(scenario)
            results.append((scenario, trace))
        
        # Assert all passed
        for scenario, trace in results:
            if trace.validation_result:
                from framework.validator import ValidationResult
                validation = ValidationResult(**trace.validation_result)
                assert validation.passed, f"Scenario {scenario.metadata.id} failed"


# Policy enforcement tests
class TestPolicyEnforcement:
    """Test that policy rules are strictly enforced."""
    
    @pytest.mark.asyncio
    async def test_intent_detector_always_first(
        self,
        all_scenarios: List[Scenario],
        scenario_executor: ScenarioExecutor
    ):
        """Intent Detector must always be the first agent."""
        for scenario in all_scenarios:
            trace = await scenario_executor.execute(scenario)
            
            assert len(trace.agent_sequence) > 0, "No agents called"
            assert trace.agent_sequence[0] == "IntentDetector", (
                f"Scenario {scenario.metadata.id}: "
                f"First agent must be IntentDetector, got {trace.agent_sequence[0]}"
            )
    
    @pytest.mark.asyncio
    async def test_ask_intent_routes_to_guru(
        self,
        scenario_loader: ScenarioLoader,
        scenario_executor: ScenarioExecutor
    ):
        """
        Ask intent must route to question-answering agent.
        
        NOTE: This test uses OLD system name "BioinformaticsGuru".
        For NEW orchestrator system, this would be "WorkflowPlannerAgent" or similar.
        See docs/AGENT_NAMING_MIGRATION.md for details.
        """
        from framework.scenario import ScenarioCategory
        from framework.executor import AGENT_NAME_MAPPING
        
        scenarios = scenario_loader.load_by_category(ScenarioCategory.ASK)
        
        for scenario in scenarios:
            trace = await scenario_executor.execute(scenario)
            
            assert len(trace.agent_sequence) >= 2, "Expected at least 2 agents"
            
            # Accept both OLD and NEW canonical names
            expected_names = ["BioinformaticsGuru", "WorkflowPlannerAgent"]
            actual_agent = trace.agent_sequence[1]
            
            assert actual_agent in expected_names, (
                f"Scenario {scenario.metadata.id}: "
                f"Ask intent must route to {' or '.join(expected_names)}, got {actual_agent}"
            )
    
    @pytest.mark.asyncio
    async def test_execute_intent_routes_to_planner(
        self,
        scenario_loader: ScenarioLoader,
        scenario_executor: ScenarioExecutor
    ):
        """
        Execute intent must route to workflow planning/execution agent.
        
        NOTE: This test accepts both OLD ("BioinformaticsExecutor") and 
        NEW ("ImplementationAgent", "WorkflowPlannerAgent") system names.
        See docs/AGENT_NAMING_MIGRATION.md for details.
        """
        from framework.scenario import ScenarioCategory
        
        scenarios = scenario_loader.load_by_category(ScenarioCategory.EXECUTE)
        
        for scenario in scenarios:
            trace = await scenario_executor.execute(scenario)
            
            assert len(trace.agent_sequence) >= 2, "Expected at least 2 agents"
            
            # Accept OLD or NEW canonical names for execution agents
            valid_execution_agents = [
                "BioinformaticsExecutor",  # OLD system
                "WorkflowPlannerAgent",     # NEW system (first step)
                "ImplementationAgent",      # NEW system (could be second agent if no planner)
            ]
            actual_agent = trace.agent_sequence[1]
            
            assert actual_agent in valid_execution_agents, (
                f"Scenario {scenario.metadata.id}: "
                f"Execute intent must route to {' or '.join(valid_execution_agents)}, got {actual_agent}"
            )
    
    @pytest.mark.asyncio
    async def test_guru_never_executes(
        self,
        scenario_loader: ScenarioLoader,
        scenario_executor: ScenarioExecutor
    ):
        """
        Question-answering agents should never trigger job execution without escalation.
        
        NOTE: This test checks for OLD system names. For NEW orchestrator system,
        execution happens via ImplementationAgent → external runner (not ExecutionBroker).
        See docs/AGENT_NAMING_MIGRATION.md for details.
        """
        from framework.scenario import ScenarioCategory
        scenarios = scenario_loader.load_by_category(ScenarioCategory.ASK)
        
        for scenario in scenarios:
            # Skip scenarios that explicitly test Guru escalation
            if "escalation" in scenario.metadata.tags:
                continue
            
            trace = await scenario_executor.execute(scenario)
            
            # If QA agent is in sequence, execution components should not be
            qa_agents = ["BioinformaticsGuru", "WorkflowPlannerAgent"]
            execution_components = ["ExecutionBroker", "ImplementationAgent"]
            
            has_qa_agent = any(agent in trace.agent_sequence for agent in qa_agents)
            has_execution = any(agent in trace.agent_sequence for agent in execution_components)
            
            if has_qa_agent:
                assert not has_execution, (
                    f"Scenario {scenario.metadata.id}: "
                    f"QA agents should not trigger execution components without escalation. "
                    f"Sequence: {' → '.join(trace.agent_sequence)}"
                )


# Baseline regression tests
class TestBaselineRegression:
    """Test that scenarios don't regress from baselines."""
    
    @pytest.mark.asyncio
    async def test_no_agent_sequence_regression(
        self,
        all_scenarios: List[Scenario],
        scenario_executor: ScenarioExecutor
    ):
        """Agent sequences should not change from baseline."""
        regressions = []
        
        for scenario in all_scenarios:
            baseline = scenario_executor.load_baseline(scenario.metadata.id)
            if not baseline:
                continue  # No baseline to compare
            
            trace = await scenario_executor.execute(scenario)
            
            if trace.agent_sequence != baseline.agent_sequence:
                regressions.append({
                    "scenario": scenario.metadata.id,
                    "baseline": baseline.agent_sequence,
                    "current": trace.agent_sequence
                })
        
        if regressions:
            details = "\n".join(
                f"  • {r['scenario']}: {' → '.join(r['baseline'])} became {' → '.join(r['current'])}"
                for r in regressions
            )
            pytest.fail(
                f"Agent sequence regressions detected:\n{details}\n"
                f"Run with --update-baseline to accept these changes."
            )


# Performance tests
class TestPerformance:
    """Test performance characteristics."""
    
    @pytest.mark.asyncio
    async def test_scenario_duration_within_bounds(
        self,
        all_scenarios: List[Scenario],
        scenario_executor: ScenarioExecutor
    ):
        """Scenarios should complete within specified time bounds."""
        slow_scenarios = []
        
        for scenario in all_scenarios:
            trace = await scenario_executor.execute(scenario)
            
            if scenario.expected_behavior.max_duration_ms:
                if trace.duration_ms > scenario.expected_behavior.max_duration_ms:
                    slow_scenarios.append({
                        "scenario": scenario.metadata.id,
                        "expected_max": scenario.expected_behavior.max_duration_ms,
                        "actual": trace.duration_ms
                    })
        
        if slow_scenarios:
            details = "\n".join(
                f"  • {s['scenario']}: {s['actual']}ms (max: {s['expected_max']}ms)"
                for s in slow_scenarios
            )
            pytest.fail(f"Slow scenarios detected:\n{details}")
