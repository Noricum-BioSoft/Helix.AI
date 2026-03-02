"""
Scenario definition models and loader.

Scenarios are the source of truth for expected multi-agent behavior.
They are defined in YAML and loaded into structured Pydantic models.
"""

import yaml
from pathlib import Path
from typing import Any, Dict, List, Optional, Literal
from pydantic import BaseModel, Field
from enum import Enum


class ScenarioCategory(str, Enum):
    """Scenario categories for organization."""
    ASK = "ask"
    EXECUTE = "execute"
    EDGE = "edge"
    MULTI_TURN = "multi_turn"


class ValidationOperator(str, Enum):
    """Operators for field validation."""
    EQUALS = "equals"
    CONTAINS = "contains"
    GREATER_THAN = "greater_than"
    LESS_THAN = "less_than"
    MIN_LENGTH = "min_length"
    MAX_LENGTH = "max_length"
    REGEX_MATCH = "regex_match"
    NOT_EMPTY = "not_empty"


class FieldValidation(BaseModel):
    """Validation rule for a contract field."""
    field: str
    operator: Optional[ValidationOperator] = None
    # Legacy support for inline operators
    equals: Optional[Any] = None
    contains: Optional[List[str]] = None
    greater_than: Optional[float] = None
    less_than: Optional[float] = None
    min_length: Optional[int] = None
    max_length: Optional[int] = None
    regex_match: Optional[str] = None
    not_empty: Optional[bool] = None
    
    def get_operator_and_value(self) -> tuple[ValidationOperator, Any]:
        """Extract operator and value from the validation."""
        if self.equals is not None:
            return ValidationOperator.EQUALS, self.equals
        elif self.contains is not None:
            return ValidationOperator.CONTAINS, self.contains
        elif self.greater_than is not None:
            return ValidationOperator.GREATER_THAN, self.greater_than
        elif self.less_than is not None:
            return ValidationOperator.LESS_THAN, self.less_than
        elif self.min_length is not None:
            return ValidationOperator.MIN_LENGTH, self.min_length
        elif self.max_length is not None:
            return ValidationOperator.MAX_LENGTH, self.max_length
        elif self.regex_match is not None:
            return ValidationOperator.REGEX_MATCH, self.regex_match
        elif self.not_empty:
            return ValidationOperator.NOT_EMPTY, True
        else:
            raise ValueError(f"No validation operator specified for field '{self.field}'")


class ContractExpectation(BaseModel):
    """Expected contract output from an agent."""
    agent: str
    output_type: str
    validations: List[FieldValidation] = Field(default_factory=list)


class PolicyCheck(str, Enum):
    """Policy checks to enforce."""
    NO_EXECUTION_SIDE_EFFECTS = "no_execution_side_effects"
    NO_INFRASTRUCTURE_DECISIONS = "no_infrastructure_decisions"
    STRICT_ROLE_SEPARATION = "strict_role_separation"
    NO_CONTRACT_MUTATION = "no_contract_mutation"
    CONFIDENCE_THRESHOLDS = "confidence_thresholds"


class IntentExpectation(BaseModel):
    """Expected intent classification."""
    type: Literal["ask", "execute"]
    confidence_min: float = Field(ge=0.0, le=1.0)


class ExpectedBehavior(BaseModel):
    """Expected behavior for a scenario."""
    agent_sequence: List[str]
    intent: Optional[IntentExpectation] = None
    contracts: List[ContractExpectation] = Field(default_factory=list)
    policy_checks: List[PolicyCheck] = Field(default_factory=list)
    min_duration_ms: Optional[int] = None
    max_duration_ms: Optional[int] = None


class ScenarioInput(BaseModel):
    """Input for a scenario."""
    user_prompt: str
    session_context: Dict[str, Any] = Field(default_factory=dict)
    uploaded_files: List[str] = Field(default_factory=list)
    environment_vars: Dict[str, str] = Field(default_factory=dict)


class ScenarioMetadata(BaseModel):
    """Metadata about a scenario."""
    id: str
    category: ScenarioCategory
    description: str
    tags: List[str] = Field(default_factory=list)
    author: Optional[str] = None
    created_at: Optional[str] = None
    requires_mock: bool = Field(
        default=True,
        description="If True, runs in mock mode without real LLM/execution"
    )


class Scenario(BaseModel):
    """
    A complete scenario definition.
    
    Scenarios are loaded from YAML files and represent expected
    system behavior for a given input.
    """
    metadata: ScenarioMetadata
    input: ScenarioInput
    expected_behavior: ExpectedBehavior
    
    @classmethod
    def from_yaml(cls, path: Path) -> "Scenario":
        """Load a scenario from a YAML file."""
        with open(path, 'r') as f:
            data = yaml.safe_load(f)
        return cls(**data)
    
    @classmethod
    def from_yaml_string(cls, yaml_str: str) -> "Scenario":
        """Load a scenario from a YAML string."""
        data = yaml.safe_load(yaml_str)
        return cls(**data)
    
    def to_yaml(self, path: Path) -> None:
        """Save a scenario to a YAML file."""
        with open(path, 'w') as f:
            yaml.safe_dump(
                self.model_dump(mode='json', exclude_none=True),
                f,
                default_flow_style=False,
                sort_keys=False
            )


class ScenarioLoader:
    """Loads scenarios from a directory."""
    
    def __init__(self, scenarios_dir: Path):
        self.scenarios_dir = scenarios_dir
    
    def load_scenario(self, scenario_id: str) -> Scenario:
        """Load a single scenario by ID."""
        # Try various naming patterns
        patterns = [
            f"{scenario_id}.yaml",
            f"{scenario_id}.yml",
            f"**/{scenario_id}.yaml",
            f"**/{scenario_id}.yml",
        ]
        
        for pattern in patterns:
            matches = list(self.scenarios_dir.glob(pattern))
            if matches:
                return Scenario.from_yaml(matches[0])
        
        raise FileNotFoundError(
            f"Scenario '{scenario_id}' not found in {self.scenarios_dir}"
        )
    
    def load_all_scenarios(self) -> List[Scenario]:
        """Load all scenarios from the directory."""
        scenarios = []
        for yaml_file in self.scenarios_dir.glob("**/*.yaml"):
            try:
                scenario = Scenario.from_yaml(yaml_file)
                scenarios.append(scenario)
            except Exception as e:
                print(f"Warning: Failed to load {yaml_file}: {e}")
        
        for yaml_file in self.scenarios_dir.glob("**/*.yml"):
            try:
                scenario = Scenario.from_yaml(yaml_file)
                scenarios.append(scenario)
            except Exception as e:
                print(f"Warning: Failed to load {yaml_file}: {e}")
        
        return scenarios
    
    def load_by_category(self, category: ScenarioCategory) -> List[Scenario]:
        """Load all scenarios in a category."""
        all_scenarios = self.load_all_scenarios()
        return [s for s in all_scenarios if s.metadata.category == category]
    
    def load_by_tag(self, tag: str) -> List[Scenario]:
        """Load all scenarios with a specific tag."""
        all_scenarios = self.load_all_scenarios()
        return [s for s in all_scenarios if tag in s.metadata.tags]
