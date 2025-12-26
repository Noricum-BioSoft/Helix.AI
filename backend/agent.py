# /backend/agent.py

import asyncio
import os
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional

from langgraph.prebuilt import create_react_agent
from langgraph.checkpoint.memory import MemorySaver

from langchain_core.tools import tool

from pydantic import BaseModel, Field

# load the environment
from dotenv import load_dotenv

logger = logging.getLogger(__name__)

try:
    # Best-effort: in some sandbox/CI environments the .env file may be unreadable
    # (e.g., ignored files are blocked). This should never be fatal.
    load_dotenv()
except PermissionError as e:
    logger.warning(f"Could not read .env (permission denied); continuing without it: {e}")
except Exception as e:
    logger.warning(f"Could not load .env; continuing without it: {e}")

# Try to import langchain globals (may not exist in all versions)
try:
    from langchain.globals import set_verbose, set_debug
    _has_langchain_globals = True
except ImportError:
    # langchain.globals doesn't exist in this version - use environment variables instead
    _has_langchain_globals = False
    def set_verbose(value: bool):
        """No-op if langchain.globals not available"""
        pass
    def set_debug(value: bool):
        """No-op if langchain.globals not available"""
        pass

from backend.prompts.templates import build_react_prompt
from backend.context_builder import build_context_snippet, build_session_brief

# Disable verbose/debug logging by default to prevent large tool results (like sequences) from being printed
# Can be enabled via LANGCHAIN_VERBOSE=true or LANGCHAIN_DEBUG=true environment variables for debugging
verbose_enabled = os.getenv("LANGCHAIN_VERBOSE", "false").lower() == "true"
debug_enabled = os.getenv("LANGCHAIN_DEBUG", "false").lower() == "true"
set_verbose(verbose_enabled)
set_debug(debug_enabled)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
AGENT_PROMPT_PATH = PROJECT_ROOT / "agent.md"

# Load system prompt (static, no placeholders)
try:
    BIOAGENT_SYSTEM_PROMPT = AGENT_PROMPT_PATH.read_text()
except Exception:
    BIOAGENT_SYSTEM_PROMPT = (
        "You are BioAgent, an autonomous bioinformatics assistant. "
        "Classify prompts, plan, execute with real tools, and return structured JSON "
        "for browser rendering. Decline non-bioinformatics or infeasible requests."
    )

class OutputFormatter(BaseModel):
    """Always use this tool to structure your response to the user."""

    input: str = Field(description="The input of the tool/function call")
    output: str = Field(description="The output of the tool/function call")
    plot: str = Field(description="An optional plot visualizing the output")


# Import all agent tools from dedicated module
from backend.agent_tools import (
    toolbox_inventory,
    sequence_alignment,
    mutate_sequence,
    dna_vendor_research,
    phylogenetic_tree,
    sequence_selection,
    synthesis_submission,
    create_session,
    plasmid_visualization,
    plasmid_for_representatives,
    single_cell_analysis,
    fetch_ncbi_sequence,
    query_uniprot,
    lookup_go_term,
    bulk_rnaseq_analysis,
    fastqc_quality_analysis,
)


_llm = None
_agent = None
_memory = None


def _get_llm():
    """
    Lazily initialize and return the LLM for BioAgent.

    Important:
    - In tests/CI we often run with `HELIX_MOCK_MODE=1` (no LLM, no network).
    - Imports for OpenAI/DeepSeek can trigger SSL/cert loading that may be blocked
      in sandboxed environments; keep them local to this function.
    """
    if os.getenv("HELIX_MOCK_MODE") == "1":
        raise RuntimeError("LLM is disabled in HELIX_MOCK_MODE")

    openai_key = os.getenv("OPENAI_API_KEY", "").strip()
    deepseek_key = os.getenv("DEEPSEEK_API_KEY", "").strip()

    openai_enabled = openai_key and openai_key not in ["", "disabled", "your_openai_api_key_here", "none"]
    deepseek_enabled = deepseek_key and deepseek_key not in ["", "disabled", "your_deepseek_api_key_here", "none"]

    if openai_enabled:
        from langchain.chat_models import init_chat_model
        return init_chat_model("openai:gpt-4o", temperature=0)

    if deepseek_enabled:
        from langchain_deepseek import ChatDeepSeek
        return ChatDeepSeek(
            model="deepseek-chat",
            temperature=0,
            max_tokens=None,
            timeout=60.0,
            max_retries=1,
        )

    raise ValueError(
        "No API keys found. Please set either OPENAI_API_KEY or DEEPSEEK_API_KEY "
        "in your environment variables."
    )


def _get_agent():
    """Lazily create and return the LangGraph agent instance."""
    global _llm, _agent, _memory
    if _agent is not None:
        return _agent

    _llm = _get_llm()
    _memory = MemorySaver()
    _agent = create_react_agent(
        model=_llm,
        tools=[
            toolbox_inventory,
            sequence_alignment,
            mutate_sequence,
            dna_vendor_research,
            phylogenetic_tree,
            sequence_selection,
            synthesis_submission,
            plasmid_visualization,
            plasmid_for_representatives,
            single_cell_analysis,
            fetch_ncbi_sequence,
            query_uniprot,
            lookup_go_term,
            bulk_rnaseq_analysis,
            fastqc_quality_analysis,
        ],
        checkpointer=_memory,
    )
    return _agent


agent = None


class CommandProcessor:
    """
    Processes user commands by mapping them to tools or generating responses.
    
    This class encapsulates the command processing logic, making it easier to test
    and maintain. Dependencies are injected to allow for mocking in tests.
    """
    
    def __init__(
        self,
        agent=None,
        intent_classifier=None,
        router=None,
        tool_generator=None,
        memory=None,
        system_prompt: str = None,
        is_mock_mode: bool = False
    ):
        """
        Initialize CommandProcessor with dependencies.
        
        Args:
            agent: LangGraph agent instance (can be None, will be lazy-loaded)
            intent_classifier: Function to classify user intent (default: classify_intent)
            router: CommandRouter instance (default: CommandRouter())
            tool_generator: Tool generator function (default: generate_and_execute_tool)
            memory: MemorySaver instance for checkpointing (default: _memory)
            system_prompt: System prompt for agent (default: BIOAGENT_SYSTEM_PROMPT)
            is_mock_mode: Whether running in mock mode (default: from HELIX_MOCK_MODE env)
        """
        self.agent = agent
        self._intent_classifier = intent_classifier
        self._router = router
        self._tool_generator = tool_generator
        self.memory = memory
        self.system_prompt = system_prompt or BIOAGENT_SYSTEM_PROMPT
        self.is_mock_mode = is_mock_mode or (os.getenv("HELIX_MOCK_MODE") == "1")
    
    def _get_intent_classifier(self):
        """Lazy load intent classifier."""
        if self._intent_classifier is None:
            from backend.intent_classifier import classify_intent
            self._intent_classifier = classify_intent
        return self._intent_classifier
    
    def _get_router(self):
        """Lazy load command router."""
        if self._router is None:
            from backend.command_router import CommandRouter
            self._router = CommandRouter()
        return self._router
    
    def _prepare_messages(self, command: str, session_context: Dict) -> tuple:
        """
        Prepare user and system messages for agent.
        
        Returns:
            Tuple of (HumanMessage, SystemMessage, session_brief)
        """
        from langchain_core.messages import HumanMessage, SystemMessage
        
        # Build session brief
        session_brief = ""
        if session_context:
            session_brief = build_session_brief(session_context, max_tokens=800)
            brief_tokens = len(session_brief) // 4
            print(f"[CommandProcessor] Session Brief: {len(session_brief):,} chars (~{brief_tokens} tokens)")
        
        # Build user content
        user_content = command
        if session_brief:
            user_content = f"## SESSION BRIEF (injected, ≤800 tokens)\n```json\n{session_brief}\n```\n\n## USER MESSAGE\n{command}"
            print(f"[CommandProcessor] Prepended Session Brief to user message")
        
        # Log message sizes
        total_message_size = len(user_content) + len(self.system_prompt)
        total_tokens = total_message_size // 4
        print(f"[CommandProcessor] Total message size to LLM: {total_message_size:,} chars (~{total_tokens:,} tokens estimated)")
        print(f"[CommandProcessor] System prompt size: {len(self.system_prompt):,} chars (static)")
        
        input_message = HumanMessage(content=user_content)
        system_message = SystemMessage(content=self.system_prompt)
        
        return input_message, system_message, session_brief
    
    def _create_session_config(self, session_id: str) -> Dict:
        """Create a temporary session config for tool mapping."""
        import uuid
        temp_session_id = f"{session_id}_mapping_{uuid.uuid4().hex[:8]}"
        return {"configurable": {"thread_id": temp_session_id}}
    
    def _deduplicate_messages(self, messages: List) -> List:
        """Remove duplicate messages from a list."""
        seen_ids = set()
        seen_content_hashes = set()
        deduplicated_messages = []
        
        for msg in messages:
            # Try to use message ID first (most reliable)
            msg_id = getattr(msg, 'id', None)
            if msg_id and msg_id in seen_ids:
                continue
            if msg_id:
                seen_ids.add(msg_id)
            
            # Fallback: use content hash for messages without IDs
            if not msg_id:
                import hashlib
                msg_type = type(msg).__name__
                msg_content = getattr(msg, 'content', str(msg))
                content_hash = hashlib.md5(f"{msg_type}:{msg_content}".encode()).hexdigest()
                if content_hash in seen_content_hashes:
                    continue
                seen_content_hashes.add(content_hash)
            
            deduplicated_messages.append(msg)
        
        if len(messages) != len(deduplicated_messages):
            print(f"[CommandProcessor] Deduplicated messages: {len(messages)} -> {len(deduplicated_messages)}")
        
        return deduplicated_messages
    
    async def _handle_qa_intent(
        self, input_message, system_message, session_config: Dict
    ) -> Dict:
        """Handle Q&A intent by letting agent respond directly."""
        print("[CommandProcessor] Skipping tool mapping for Q&A intent, letting agent respond directly...")
        
        if self.agent is None:
            self.agent = _get_agent()
        
        result = await asyncio.wait_for(
            asyncio.to_thread(self.agent.invoke, {"messages": [system_message, input_message]}, session_config),
            timeout=60.0  # Longer timeout for Q&A responses
        )
        
        # Deduplicate messages
        if isinstance(result, dict) and "messages" in result:
            result["messages"] = self._deduplicate_messages(result.get("messages", []))
        
        return result
    
    def _extract_tool_from_stream_event(self, event: Dict) -> Optional[Dict]:
        """Extract tool mapping from a single stream event."""
        if not isinstance(event, dict):
            return None
        
        for node_name, node_output in event.items():
            # Look for "tools" node
            if node_name == "tools" or "tools" in node_name.lower():
                if isinstance(node_output, dict) and "messages" in node_output:
                    messages = node_output["messages"]
                    for msg in messages:
                        if hasattr(msg, 'name'):
                            tool_name = msg.name
                            tool_args = {}
                            if hasattr(msg, 'content'):
                                content = msg.content
                                if isinstance(content, str):
                                    import json
                                    try:
                                        tool_args = json.loads(content)
                                    except:
                                        pass
                            elif hasattr(msg, 'input'):
                                tool_args = msg.input
                            
                            return {
                                "tool_name": tool_name,
                                "parameters": tool_args if isinstance(tool_args, dict) else {}
                            }
            
            # Check for AIMessage with tool_calls
            if isinstance(node_output, dict) and "messages" in node_output:
                messages = node_output["messages"]
                for msg in messages:
                    if hasattr(msg, 'tool_calls') and msg.tool_calls:
                        for tool_call in msg.tool_calls:
                            tool_name = None
                            tool_args = {}
                            
                            if hasattr(tool_call, 'name'):
                                tool_name = tool_call.name
                            elif isinstance(tool_call, dict):
                                tool_name = tool_call.get("name")
                            
                            if hasattr(tool_call, 'args'):
                                tool_args = tool_call.args
                            elif isinstance(tool_call, dict):
                                tool_args = tool_call.get("args", {})
                            
                            if tool_name:
                                return {
                                    "tool_name": tool_name,
                                    "parameters": tool_args if isinstance(tool_args, dict) else {}
                                }
        return None
    
    async def _extract_tool_mapping_from_stream(
        self, input_message, system_message, session_config: Dict
    ) -> Optional[Dict]:
        """Extract tool mapping by streaming agent execution."""
        if self.agent is None:
            self.agent = _get_agent()
        
        def capture_tool_call_sync():
            """Synchronous function to stream agent and capture first tool call."""
            try:
                for event in self.agent.stream(
                    {"messages": [system_message, input_message]},
                    session_config,
                    stream_mode="updates"
                ):
                    tool_mapping = self._extract_tool_from_stream_event(event)
                    if tool_mapping:
                        print(f"[CommandProcessor] Agent mapped tool from stream: {tool_mapping['tool_name']} with params: {tool_mapping['parameters']}")
                        return tool_mapping
            except Exception as e:
                print(f"[CommandProcessor] Error during tool call capture: {e}")
                import traceback
                traceback.print_exc()
            return None
        
        try:
            tool_mapping = await asyncio.wait_for(
                asyncio.to_thread(capture_tool_call_sync),
                timeout=30.0
            )
        except asyncio.TimeoutError:
            print("⚠️  Agent tool mapping timed out after 30s")
            # Try to extract from state
            if self.memory:
                extracted = _extract_tool_call_from_state(self.memory, session_config)
                if extracted:
                    tool_mapping = extracted
        
        return tool_mapping
    
    def _extract_tool_from_result_messages(self, result: Dict) -> Optional[Dict]:
        """Extract tool mapping from agent result messages."""
        if not isinstance(result, dict) or "messages" not in result:
            return None
        
        messages = result.get("messages", [])
        for msg in reversed(messages):
            if hasattr(msg, 'tool_calls') and msg.tool_calls:
                for tool_call in msg.tool_calls:
                    tool_name = getattr(tool_call, 'name', None)
                    tool_args = getattr(tool_call, 'args', None) or {}
                    if tool_name:
                        print(f"[CommandProcessor] Extracted tool mapping from result: {tool_name}")
                        return {
                            "tool_name": tool_name,
                            "parameters": tool_args
                        }
        return None
    
    def _build_tool_mapped_response(self, tool_mapping: Dict) -> Dict:
        """Build a standardized tool_mapped response."""
        return {
            "tool_mapping": tool_mapping,
            "tool_name": tool_mapping["tool_name"],
            "parameters": tool_mapping["parameters"],
            "status": "tool_mapped",
            "message": f"Agent identified tool: {tool_mapping['tool_name']}. Execution will be handled by router."
        }
    
    async def _try_router_fallback(self, command: str, session_context: Dict) -> Optional[Dict]:
        """Try deterministic router fallback for safe tools."""
        try:
            router = self._get_router()
            tool_name, parameters = router.route_command(command, session_context or {})
            
            safe_tools = {
                "toolbox_inventory",
                "read_merging",
                "read_trimming",
                "fastqc_quality_analysis",
                "sequence_alignment",
                "mutate_sequence",
                "plasmid_visualization",
                "phylogenetic_tree",
                "clustering_analysis",
                "variant_selection",
                "fetch_ncbi_sequence",
                "query_uniprot",
                "lookup_go_term",
                "bulk_rnaseq_analysis",
                "single_cell_analysis",
            }
            
            if tool_name in safe_tools:
                print(f"[CommandProcessor] Deterministic router fallback -> {tool_name}")
                return {
                    "tool_mapping": {"tool_name": tool_name, "parameters": parameters or {}},
                    "tool_name": tool_name,
                    "parameters": parameters or {},
                    "status": "tool_mapped",
                    "message": f"Deterministic router fallback identified tool: {tool_name}. Execution will be handled by router.",
                }
        except Exception:
            pass
        return None
    
    async def _try_tool_generator_fallback(
        self, command: str, session_id: str, session_context: Dict, intent
    ) -> Optional[Dict]:
        """Try tool generator agent as final fallback."""
        if intent.intent != "execute":
            print(f"[CommandProcessor] Skipping tool-generator-agent (intent={intent.intent}, reason={intent.reason})")
            return None
        
        try:
            if self._tool_generator is None:
                from backend.tool_generator_agent import generate_and_execute_tool, _discover_inputs_from_args, _discover_outputs_from_args
            else:
                # If injected, assume it's the function itself
                generate_and_execute_tool = self._tool_generator
                from backend.tool_generator_agent import _discover_inputs_from_args, _discover_outputs_from_args
            
            command_args = {"command": command}
            discovered_inputs = _discover_inputs_from_args(command_args, session_context)
            if discovered_inputs:
                print(f"[CommandProcessor] 🔧 Discovered {len(discovered_inputs)} input files for infrastructure decision")
            
            discovered_outputs = _discover_outputs_from_args(command_args, command)
            if discovered_outputs:
                print(f"[CommandProcessor] 🔧 Discovered {len(discovered_outputs)} output paths")
            
            tool_result = await generate_and_execute_tool(
                command=command,
                user_request=command,
                session_id=session_id,
                inputs=discovered_inputs,
                outputs=discovered_outputs
            )
            
            if tool_result.get("status") == "success":
                print(f"[CommandProcessor] ✅ Tool-generator-agent successfully generated tool")
                return {
                    "status": "success",
                    "tool_generated": True,
                    "tool_name": "generated_tool",
                    "result": tool_result,
                    "message": tool_result.get("explanation", "Tool generated and executed successfully")
                }
            else:
                error_msg = tool_result.get('error')
                if not error_msg:
                    execution_result = tool_result.get("execution_result", {})
                    if isinstance(execution_result, dict):
                        error_msg = execution_result.get("error")
                        if not error_msg and execution_result.get("stderr"):
                            error_msg = f"Execution failed: {execution_result.get('stderr', '')[:200]}"
                    error_msg = error_msg or "Unknown error"
                print(f"[CommandProcessor] ⚠️  Tool-generator-agent failed: {error_msg}")
        except Exception as e:
            print(f"[CommandProcessor] ❌ Tool-generator-agent exception: {e}")
            import traceback
            traceback.print_exc()
        
        return None
    
    async def process(
        self, command: str, session_id: str = "default", session_context: Dict = None
    ) -> Dict:
        """
        Main processing method that orchestrates command handling.
        
        Returns:
            Dict with tool mapping or agent response
        """
        print(f"[CommandProcessor] Processing command: {command}")
        print(f"[CommandProcessor] session_id: {session_id}")
        
        # Classify intent (using injected classifier if provided, otherwise default)
        intent = self._get_intent_classifier()(command)
        
        # Prepare messages
        input_message, system_message, session_brief = self._prepare_messages(command, session_context or {})
        
        # Create session config
        temp_session_config = self._create_session_config(session_id)
        
        # Handle mock mode
        if self.is_mock_mode:
            if self.agent is None:
                return {
                    "tool_mapping": {"tool_name": "toolbox_inventory", "parameters": {}},
                    "tool_name": "toolbox_inventory",
                    "parameters": {},
                    "status": "tool_mapped",
                    "message": "HELIX_MOCK_MODE=1: agent disabled; returning toolbox_inventory mapping.",
                }
            if hasattr(self.agent, "invoke"):
                return self.agent.invoke({"messages": [system_message, input_message]}, config=None)
        
        # Ensure agent exists
        if self.agent is None:
            self.agent = _get_agent()
        
        # Handle Q&A intent
        if intent.intent == "qa":
            return await self._handle_qa_intent(input_message, system_message, temp_session_config)
        
        # Try to extract tool mapping from stream
        try:
            tool_mapping = await self._extract_tool_mapping_from_stream(
                input_message, system_message, temp_session_config
            )
            
            if tool_mapping and tool_mapping.get("tool_name"):
                print(f"[CommandProcessor] Returning tool mapping (no execution): {tool_mapping['tool_name']}")
                return self._build_tool_mapped_response(tool_mapping)
            
            # Fallback: try regular invoke
            print("[CommandProcessor] Streaming didn't capture tool call, using regular invoke as fallback...")
            import uuid
            fallback_session_id = f"{session_id}_mapping_{uuid.uuid4().hex[:8]}"
            fallback_session_config = {"configurable": {"thread_id": fallback_session_id}}
            
            result = await asyncio.wait_for(
                asyncio.to_thread(self.agent.invoke, {"messages": [system_message, input_message]}, fallback_session_config),
                timeout=30.0
            )
            
            # Deduplicate messages
            if isinstance(result, dict) and "messages" in result:
                result["messages"] = self._deduplicate_messages(result.get("messages", []))
            
            # Try to extract tool from result
            tool_mapping = self._extract_tool_from_result_messages(result)
            if tool_mapping:
                return self._build_tool_mapped_response(tool_mapping)
            
            # Try router fallback
            router_result = await self._try_router_fallback(command, session_context or {})
            if router_result:
                return router_result
            
            # Try tool generator fallback
            toolgen_result = await self._try_tool_generator_fallback(command, session_id, session_context or {}, intent)
            if toolgen_result:
                return toolgen_result
            
            # Return original result
            msg_count = len(result.get("messages", [])) if isinstance(result, dict) else 0
            print(f"[CommandProcessor] result: {msg_count} messages in response (no tool mapping extracted)")
            return result
            
        except asyncio.TimeoutError:
            print("⚠️  Agent tool mapping timed out after 30s")
            print("🔍 Attempting to extract tool call from agent state...")
            if self.memory:
                extracted = _extract_tool_call_from_state(self.memory, temp_session_config)
                if extracted:
                    return self._build_tool_mapped_response(extracted)
            raise TimeoutError("Agent tool mapping timed out and could not extract tool call")
        
        except Exception as e:
            error_str = str(e)
            if any(keyword in error_str for keyword in ["Connection error", "APIConnectionError", "timeout", "INVALID_CHAT_HISTORY"]):
                print(f"⚠️  Agent tool mapping failed: {e}")
                print("🔍 Attempting to extract tool call from agent state...")
                if self.memory:
                    extracted = _extract_tool_call_from_state(self.memory, temp_session_config)
                    if extracted:
                        return self._build_tool_mapped_response(extracted)
            raise


async def handle_command(command: str, session_id: str = "default", session_context: Dict = None):
    """
    Use agent for tool mapping only - identify which tool to use and parameters.
    Tool execution is handled by the router/tool path in main_with_mcp.py.
    
    This is a thin wrapper around CommandProcessor for backward compatibility.
    """
    global agent, _memory
    print(f"[handle_command] command: {command}")
    print(f"[handle_command] session_id: {session_id}")
    print(f"[handle_command] Agent mode: TOOL MAPPING ONLY (no execution)")
    
    # Create processor with current global state
    processor = CommandProcessor(
        agent=agent,
        memory=_memory,
        is_mock_mode=(os.getenv("HELIX_MOCK_MODE") == "1")
    )
    
    # Process command
    result = await processor.process(command, session_id, session_context)
    
    # Update global agent reference if it was lazy-loaded
    if processor.agent is not None and agent is None:
        agent = processor.agent
    
    return result


def _extract_tool_call_from_state(memory, session_config):
    """Extract tool call (name + params) from agent state before execution."""
    try:
        thread_id = session_config.get("configurable", {}).get("thread_id", "default")
        
        # Access memory storage
        if hasattr(memory, 'storage'):
            storage = memory.storage
        elif hasattr(memory, '_storage'):
            storage = memory._storage
        else:
            return None
        
        if thread_id in storage:
            checkpoints = storage.get(thread_id, [])
            if checkpoints:
                latest_checkpoint = checkpoints[-1] if isinstance(checkpoints, list) else checkpoints
                if isinstance(latest_checkpoint, tuple) and len(latest_checkpoint) >= 2:
                    checkpoint = latest_checkpoint[1]
                else:
                    checkpoint = latest_checkpoint
                
                if checkpoint and isinstance(checkpoint, dict):
                    if "channel_values" in checkpoint:
                        messages = checkpoint["channel_values"].get("messages", [])
                        # Look for AIMessage with tool_calls
                        for msg in reversed(messages):
                            if hasattr(msg, 'tool_calls') and msg.tool_calls:
                                for tool_call in msg.tool_calls:
                                    tool_name = getattr(tool_call, 'name', None)
                                    tool_args = getattr(tool_call, 'args', None) or {}
                                    if tool_name:
                                        return {
                                            "tool_name": tool_name,
                                            "parameters": tool_args
                                        }
    except Exception as e:
        print(f"[handle_command] Error extracting tool call from state: {e}")
    return None


def _extract_tool_results_from_state(memory, session_config, error_msg):
    """Extract tool results from agent state when LLM call fails."""
    try:
        thread_id = session_config.get("configurable", {}).get("thread_id", "default")
        
        # MemorySaver stores checkpoints in a dict keyed by thread_id
        # Access the internal storage (this is implementation-specific but works for MemorySaver)
        if hasattr(memory, 'storage'):
            storage = memory.storage
        elif hasattr(memory, '_storage'):
            storage = memory._storage
        else:
            # Try to get via the public API if available
            try:
                # MemorySaver might have a list() or get() method
                checkpoints = list(memory.list(session_config)) if hasattr(memory, 'list') else []
                if checkpoints:
                    # Get the latest checkpoint
                    latest = checkpoints[-1] if checkpoints else None
                    if latest:
                        state = memory.get(session_config, latest) if hasattr(memory, 'get') else None
                        if state and "channel_values" in state:
                            messages = state["channel_values"].get("messages", [])
                            return _extract_from_messages(messages, error_msg)
            except Exception:
                pass
            print("⚠️  Could not access memory storage")
            return None
        
        if thread_id in storage:
            # MemorySaver stores checkpoints as a list of (config, checkpoint) tuples
            # Get the latest checkpoint for this thread
            checkpoints = storage.get(thread_id, [])
            if checkpoints:
                latest_checkpoint = checkpoints[-1] if isinstance(checkpoints, list) else checkpoints
                if isinstance(latest_checkpoint, tuple) and len(latest_checkpoint) >= 2:
                    checkpoint = latest_checkpoint[1]
                else:
                    checkpoint = latest_checkpoint
                
                if checkpoint and isinstance(checkpoint, dict):
                    # Checkpoint structure: {"channel_values": {"messages": [...]}}
                    if "channel_values" in checkpoint:
                        messages = checkpoint["channel_values"].get("messages", [])
                        return _extract_from_messages(messages, error_msg)
        
        print("⚠️  Could not extract tool results from state")
        return None
    except Exception as e:
        print(f"[handle_command] Error extracting tool results from state: {e}")
        return None


def _extract_from_messages(messages, error_msg):
    """Helper to extract tool results from messages."""
    # Look for ToolMessage with results
    for msg in reversed(messages):
        if hasattr(msg, 'name') and hasattr(msg, 'content'):
            # This is a ToolMessage with results
            return {
                "tool_name": msg.name,
                "result": msg.content,
                "error": error_msg
            }
    return None


def _extract_tool_call_from_state(memory, session_config):
    """Extract tool call (name + params) from agent state before execution."""
    try:
        thread_id = session_config.get("configurable", {}).get("thread_id", "default")
        
        # Access memory storage
        if hasattr(memory, 'storage'):
            storage = memory.storage
        elif hasattr(memory, '_storage'):
            storage = memory._storage
        else:
            return None
        
        if thread_id in storage:
            checkpoints = storage.get(thread_id, [])
            if checkpoints:
                latest_checkpoint = checkpoints[-1] if isinstance(checkpoints, list) else checkpoints
                if isinstance(latest_checkpoint, tuple) and len(latest_checkpoint) >= 2:
                    checkpoint = latest_checkpoint[1]
                else:
                    checkpoint = latest_checkpoint
                
                if checkpoint and isinstance(checkpoint, dict):
                    if "channel_values" in checkpoint:
                        messages = checkpoint["channel_values"].get("messages", [])
                        # Look for AIMessage with tool_calls
                        for msg in reversed(messages):
                            if hasattr(msg, 'tool_calls') and msg.tool_calls:
                                for tool_call in msg.tool_calls:
                                    tool_name = getattr(tool_call, 'name', None)
                                    tool_args = getattr(tool_call, 'args', None) or {}
                                    if tool_name:
                                        return {
                                            "tool_name": tool_name,
                                            "parameters": tool_args
                                        }
    except Exception as e:
        print(f"[handle_command] Error extracting tool call from state: {e}")
    return None


def _extract_tool_results_from_state(memory, session_config, error_msg):
    """Extract tool results from agent state when LLM call fails."""
    try:
        thread_id = session_config.get("configurable", {}).get("thread_id", "default")
        
        # MemorySaver stores checkpoints in a dict keyed by thread_id
        # Access the internal storage (this is implementation-specific but works for MemorySaver)
        if hasattr(memory, 'storage'):
            storage = memory.storage
        elif hasattr(memory, '_storage'):
            storage = memory._storage
        else:
            # Try to get via the public API if available
            try:
                # MemorySaver might have a list() or get() method
                checkpoints = list(memory.list(session_config)) if hasattr(memory, 'list') else []
                if checkpoints:
                    # Get the latest checkpoint
                    latest = checkpoints[-1] if checkpoints else None
                    if latest:
                        state = memory.get(session_config, latest) if hasattr(memory, 'get') else None
                        if state and "channel_values" in state:
                            messages = state["channel_values"].get("messages", [])
                            return _extract_from_messages(messages, error_msg)
            except Exception:
                pass
            print("⚠️  Could not access memory storage")
            return None
        
        if thread_id in storage:
            # MemorySaver stores checkpoints as a list of (config, checkpoint) tuples
            # Get the latest checkpoint for this thread
            checkpoints = storage.get(thread_id, [])
            if checkpoints:
                latest_checkpoint = checkpoints[-1] if isinstance(checkpoints, list) else checkpoints
                if isinstance(latest_checkpoint, tuple) and len(latest_checkpoint) >= 2:
                    checkpoint = latest_checkpoint[1]
                else:
                    checkpoint = latest_checkpoint
                
                if checkpoint and isinstance(checkpoint, dict):
                    # Checkpoint structure: {"channel_values": {"messages": [...]}}
                    if "channel_values" in checkpoint:
                        messages = checkpoint["channel_values"].get("messages", [])
                        return _extract_from_messages(messages, error_msg)
        
        print("⚠️  Could not extract tool results from state")
        return None
    except Exception as extract_err:
        print(f"❌ Could not extract tool results: {extract_err}")
        import traceback
        traceback.print_exc()
        return None


def _extract_from_messages(messages, error_msg):
    """Extract tool results from message list."""
    # Look for the last tool message (which contains the actual results)
    for msg in reversed(messages):
        if hasattr(msg, 'type') and msg.type == 'tool':
            tool_name = getattr(msg, 'name', 'unknown')
            print(f"✅ Found tool result in agent state: {tool_name}")
            # Extract the tool result
            tool_result = msg.content if hasattr(msg, 'content') else str(msg)
            return {
                "messages": messages,
                "tool_result": tool_result,
                "tool_name": tool_name,
                "status": "partial_success",
                "error": error_msg
            }
    return None


if __name__ == "__main__":
    seq_align_results = asyncio.run(handle_command("align the given sequences"))
    print(seq_align_results)

    mut_results = asyncio.run(handle_command("mutate the given sequence."))
    print(mut_results)


def _maybe_handle_vendor_request(command: str):
    """Detect vendor-focused requests and handle them directly."""
    command_lower = command.lower()
    
    # Exclude "test" if it's part of a file path (e.g., s3://bucket/test/file.fq)
    # Check if "test" appears in an S3 path or file path context
    is_test_in_path = bool(re.search(r's3://[^/]+/[^/]*test[^/]*/', command_lower) or 
                           re.search(r'[/\\][^/\\]*test[^/\\]*[/\\]', command_lower))
    
    # FastQC detection (HIGH PRIORITY - before vendor check)
    if any(phrase in command_lower for phrase in ['fastqc', 'fastqc analysis', 'quality control', 'quality analysis', 'perform fastqc']):
        return None  # Let it be handled by the agent's tool selection
    
    vendor_keywords = ['order', 'vendor', 'synthesis', 'assay', 'expression', 'function', 'binding']
    # Only include "test" if it's not in a file path
    if not is_test_in_path:
        vendor_keywords.append('test')
    
    if not any(keyword in command_lower for keyword in vendor_keywords):
        return None

    print("[handle_command] Detected vendor research command, using dna_vendor_research tool")

    sequence_length = None
    if '96' in command or 'variants' in command:
        sequence_length = 1000

    from dna_vendor_research import run_dna_vendor_research_raw

    result = run_dna_vendor_research_raw(command, sequence_length, 'large')

    return {
        "text": f"DNA vendor research completed: {result.get('message', 'Research done')}",
        "input": command,
        "output": result,
        "plot": {
            "data": [
                {
                    "x": ["Vendors", "Testing Options"],
                    "y": [result.get('total_vendors', 0), result.get('total_testing_options', 0)],
                    "type": "bar",
                }
            ],
            "layout": {"title": "DNA Vendor Research Results"},
        },
    }
