"""
Workflow Checkpoint Manager
===========================
Provides an explicit, persisted state machine for every session.

States
------
IDLE                    – No workflow in progress.
WAITING_FOR_CLARIFICATION – Planner asked a question; waiting for the user's answer.
WAITING_FOR_INPUTS      – Plan is ready but required data inputs are missing.
WAITING_FOR_APPROVAL    – Plan is staged and requires explicit user approval before execution.
PLANNING                – Planner is actively constructing a plan (transient; cleared on save).
READY_TO_EXECUTE        – All inputs bound, no approval required; will execute on next turn.
EXECUTING               – Pipeline is running.
FAILED                  – Last step failed; recovery options available.
FAILED_WAITING_FOR_USER – Failure recovery needs user guidance before retrying.
COMPLETED               – Workflow finished successfully.
"""

from __future__ import annotations

import time
import uuid
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional


class WorkflowState(str, Enum):
    IDLE = "IDLE"
    WAITING_FOR_CLARIFICATION = "WAITING_FOR_CLARIFICATION"
    WAITING_FOR_INPUTS = "WAITING_FOR_INPUTS"
    WAITING_FOR_APPROVAL = "WAITING_FOR_APPROVAL"
    PLANNING = "PLANNING"
    READY_TO_EXECUTE = "READY_TO_EXECUTE"
    EXECUTING = "EXECUTING"
    FAILED = "FAILED"
    FAILED_WAITING_FOR_USER = "FAILED_WAITING_FOR_USER"
    COMPLETED = "COMPLETED"


# States where the system is actively waiting for the user to act.
WAITING_STATES = {
    WorkflowState.WAITING_FOR_CLARIFICATION,
    WorkflowState.WAITING_FOR_INPUTS,
    WorkflowState.WAITING_FOR_APPROVAL,
    WorkflowState.FAILED_WAITING_FOR_USER,
}

# States that block a new independent workflow from starting.
BLOCKING_STATES = WAITING_STATES | {WorkflowState.EXECUTING}


@dataclass
class WorkflowCheckpoint:
    """Persisted workflow state for a session.

    Stored under session["__checkpoint__"] in history_manager.
    """

    workflow_id: str = field(default_factory=lambda: str(uuid.uuid4()))
    state: WorkflowState = WorkflowState.IDLE

    # What the system is waiting for
    pending_plan: Optional[Dict[str, Any]] = None        # staged plan
    pending_question: Optional[str] = None               # clarification question asked
    missing_inputs: Optional[List[str]] = None           # list of missing param names
    failed_step: Optional[str] = None                    # step id that failed
    failure_message: Optional[str] = None                # error text from failed step

    # Execution tracking
    run_id: Optional[str] = None
    current_step: Optional[str] = None
    resume_node: Optional[str] = None                    # logical label for where to resume

    # Timestamps (unix seconds)
    created_at: float = field(default_factory=time.time)
    updated_at: float = field(default_factory=time.time)

    def transition(self, new_state: WorkflowState) -> "WorkflowCheckpoint":
        """Return a new checkpoint with state updated and timestamp bumped."""
        import copy
        cp = copy.copy(self)
        cp.state = new_state
        cp.updated_at = time.time()
        return cp

    def to_dict(self) -> Dict[str, Any]:
        return {
            "workflow_id": self.workflow_id,
            "state": self.state.value,
            "pending_plan": self.pending_plan,
            "pending_question": self.pending_question,
            "missing_inputs": self.missing_inputs,
            "failed_step": self.failed_step,
            "failure_message": self.failure_message,
            "run_id": self.run_id,
            "current_step": self.current_step,
            "resume_node": self.resume_node,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "WorkflowCheckpoint":
        return cls(
            workflow_id=d.get("workflow_id") or str(uuid.uuid4()),
            state=WorkflowState(d.get("state", WorkflowState.IDLE.value)),
            pending_plan=d.get("pending_plan"),
            pending_question=d.get("pending_question"),
            missing_inputs=d.get("missing_inputs"),
            failed_step=d.get("failed_step"),
            failure_message=d.get("failure_message"),
            run_id=d.get("run_id"),
            current_step=d.get("current_step"),
            resume_node=d.get("resume_node"),
            created_at=d.get("created_at", time.time()),
            updated_at=d.get("updated_at", time.time()),
        )

    @classmethod
    def idle(cls) -> "WorkflowCheckpoint":
        return cls(state=WorkflowState.IDLE)

    # ── convenience factories ────────────────────────────────────────────────

    @classmethod
    def waiting_for_approval(
        cls,
        pending_plan: Dict[str, Any],
        *,
        resume_node: str = "execute_pipeline",
    ) -> "WorkflowCheckpoint":
        return cls(
            state=WorkflowState.WAITING_FOR_APPROVAL,
            pending_plan=pending_plan,
            resume_node=resume_node,
        )

    @classmethod
    def waiting_for_inputs(
        cls,
        pending_plan: Dict[str, Any],
        missing_inputs: List[str],
        *,
        resume_node: str = "execute_pipeline",
    ) -> "WorkflowCheckpoint":
        return cls(
            state=WorkflowState.WAITING_FOR_INPUTS,
            pending_plan=pending_plan,
            missing_inputs=missing_inputs,
            resume_node=resume_node,
        )

    @classmethod
    def waiting_for_clarification(
        cls,
        question: str,
        *,
        resume_node: str = "planner",
    ) -> "WorkflowCheckpoint":
        return cls(
            state=WorkflowState.WAITING_FOR_CLARIFICATION,
            pending_question=question,
            resume_node=resume_node,
        )

    @classmethod
    def failed(
        cls,
        failed_step: str,
        failure_message: str,
        *,
        run_id: Optional[str] = None,
        pending_plan: Optional[Dict[str, Any]] = None,
    ) -> "WorkflowCheckpoint":
        return cls(
            state=WorkflowState.FAILED_WAITING_FOR_USER,
            failed_step=failed_step,
            failure_message=failure_message,
            run_id=run_id,
            pending_plan=pending_plan,
            resume_node="failure_recovery",
        )

    @classmethod
    def executing(cls, run_id: str, current_step: str = "step1") -> "WorkflowCheckpoint":
        return cls(
            state=WorkflowState.EXECUTING,
            run_id=run_id,
            current_step=current_step,
        )

    @classmethod
    def completed(cls, run_id: Optional[str] = None) -> "WorkflowCheckpoint":
        return cls(state=WorkflowState.COMPLETED, run_id=run_id)
