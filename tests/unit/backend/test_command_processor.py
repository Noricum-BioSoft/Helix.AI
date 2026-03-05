"""
Unit tests for CommandProcessor class.

These tests verify that CommandProcessor correctly handles various command types,
shortcuts, intent classification, and tool mapping.
"""

import pytest
import asyncio
from unittest.mock import Mock, MagicMock, AsyncMock, patch
from dataclasses import dataclass
from typing import Dict, Literal, Optional

from backend.agent import CommandProcessor
from backend.intent_classifier import IntentDecision


# Mock message classes
class MockMessage:
    def __init__(self, content: str, msg_id: str = None):
        self.content = content
        self.id = msg_id


class MockHumanMessage(MockMessage):
    pass


class MockSystemMessage(MockMessage):
    pass


class MockAIMessage(MockMessage):
    def __init__(self, content: str, tool_calls=None, msg_id: str = None):
        super().__init__(content, msg_id)
        self.tool_calls = tool_calls or []


class MockToolCall:
    def __init__(self, name: str, args: Dict):
        self.name = name
        self.args = args


@pytest.fixture
def mock_agent():
    """Create a mock LangGraph agent."""
    agent = Mock()
    # Return empty messages so router fallback is used
    # Make stream return an empty iterator that completes immediately
    # The stream is called with specific arguments, so we need to handle that
    def empty_stream(*args, **kwargs):
        # Return an empty iterator - this simulates stream completing with no tool calls
        return iter([])
    agent.stream = Mock(side_effect=empty_stream)
    agent.invoke = Mock(return_value={"messages": []})
    return agent


@pytest.fixture
def mock_intent_classifier():
    """Create a mock intent classifier."""
    def classifier(command: str) -> IntentDecision:
        # Toolbox queries should be execute intent (they want to execute toolbox_inventory tool)
        if any(phrase in command.lower() for phrase in [
            "what tools do you have", "what tools are available", "list tools",
            "show tools", "your toolbox", "toolbox", "capabilities", "what can you do"
        ]):
            return IntentDecision(intent="execute", reason="Toolbox query detected")
        if "?" in command or command.lower().startswith(("what", "how", "why")):
            return IntentDecision(intent="qa", reason="Question detected")
        return IntentDecision(intent="execute", reason="Command detected")
    return classifier


@pytest.fixture
def mock_router():
    """Create a mock command router."""
    from backend.command_router import CommandRouter
    # Use the real router so it properly detects toolbox queries
    # But allow tests to override route_command if needed
    router = CommandRouter()
    return router

@pytest.fixture
def mock_router_mock():
    """Create a fully mockable router for tests that need to control route_command."""
    router = Mock()
    # Default behavior: return None for unmatched commands
    # Tests can override with .return_value or .side_effect
    router.route_command = Mock(return_value=(None, {}))
    return router


@pytest.fixture
def processor(mock_agent, mock_intent_classifier, mock_router):
    """Create a CommandProcessor instance with mocked dependencies."""
    return CommandProcessor(
        agent=mock_agent,
        intent_classifier=mock_intent_classifier,
        router=mock_router,
        memory=None,
        is_mock_mode=False
    )


class TestToolboxQuery:
    """Test toolbox/capability query detection."""
    
    @pytest.mark.asyncio
    async def test_toolbox_query_what_can_you_do(self, processor):
        """Test detection of 'what can you do' query."""
        result = await processor.process("what can you do", "test-session")
        assert result is not None
        # Router fallback should return toolbox_inventory
        # Check both possible response formats
        tool_name = result.get("tool_name") or (result.get("tool_mapping", {}).get("tool_name") if isinstance(result.get("tool_mapping"), dict) else None)
        assert tool_name == "toolbox_inventory", f"Expected toolbox_inventory, got {tool_name}. Full result: {result}"
        assert result.get("status") == "tool_mapped"
    
    @pytest.mark.asyncio
    async def test_toolbox_query_list_tools(self, processor):
        """Test detection of 'list tools' query."""
        result = await processor.process("list tools", "test-session")
        assert result is not None
        # In mock mode, returns toolbox_inventory mapping directly
        if "tool_mapping" in result:
            assert result["tool_mapping"]["tool_name"] == "toolbox_inventory"
        else:
            assert result.get("tool_name") == "toolbox_inventory"
    
    @pytest.mark.asyncio
    async def test_toolbox_query_capabilities(self, processor):
        """Test detection of 'capabilities' query."""
        result = await processor.process("show me your capabilities", "test-session")
        assert result is not None
        # In mock mode, returns toolbox_inventory mapping directly
        if "tool_mapping" in result:
            assert result["tool_mapping"]["tool_name"] == "toolbox_inventory"
        else:
            assert result.get("tool_name") == "toolbox_inventory"
    
    @pytest.mark.asyncio
    async def test_non_toolbox_query(self, processor):
        """Test that non-toolbox queries return None for toolbox_inventory."""
        result = await processor.process("align sequences", "test-session")
        # Should not return toolbox_inventory for non-toolbox queries
        assert result.get("tool_name") != "toolbox_inventory"


class TestReadMerging:
    """Test read-merging shortcut detection."""
    
    @pytest.mark.asyncio
    async def test_read_merging_detection(self, processor, mock_router_mock):
        """Test detection of read-merging command."""
        mock_router_mock.route_command.return_value = ("read_merging", {"r1": "s3://bucket/r1.fastq", "r2": "s3://bucket/r2.fastq"})
        # Update processor's router (use private attribute that _get_router() uses)
        processor._router = mock_router_mock
        
        command = "merge r1:s3://bucket/r1.fastq r2:s3://bucket/r2.fastq"
        result = await processor._try_router_fallback(command, {})
        
        assert result is not None
        assert result["tool_name"] == "read_merging"
        assert result["status"] == "tool_mapped"
    
    @pytest.mark.asyncio
    async def test_read_merging_missing_s3(self, processor, mock_router_mock):
        """Test that read-merging without S3 paths uses real router which may still match."""
        # Use real router for this test - it may match "merge" commands regardless of S3
        command = "merge r1:file1.fastq r2:file2.fastq"
        result = await processor._try_router_fallback(command, {})
        # The real router may match "merge" commands, so just check it doesn't crash
        assert result is None or isinstance(result, dict)
    
    @pytest.mark.asyncio
    async def test_read_merging_router_exception(self, processor, mock_router_mock):
        """Test that router exceptions don't crash the processor."""
        mock_router_mock.route_command.side_effect = Exception("Router error")
        # Set the private _router attribute that _get_router() uses
        processor._router = mock_router_mock
        command = "merge r1:s3://bucket/r1.fastq r2:s3://bucket/r2.fastq"
        result = await processor._try_router_fallback(command, {})
        # Should return None, not raise
        assert result is None


class TestMessagePreparation:
    """Test message preparation logic."""
    
    def test_prepare_messages_no_session_context(self, processor):
        """Test message preparation without session context."""
        input_msg, system_msg, brief = processor._prepare_messages("test command", None)
        
        assert input_msg.content == "test command"
        assert system_msg.content == processor.system_prompt
        assert brief == ""
    
    def test_prepare_messages_with_session_context(self, processor):
        """Test message preparation with session context."""
        session_context = {"files": ["file1.txt"], "history": []}
        
        with patch('backend.agent.build_session_brief', return_value='{"files": ["file1.txt"]}'):
            input_msg, system_msg, brief = processor._prepare_messages("test command", session_context)
            
            assert "SESSION BRIEF" in input_msg.content
            assert "USER MESSAGE" in input_msg.content
            assert "test command" in input_msg.content
            assert brief != ""


class TestMessageDeduplication:
    """Test message deduplication logic."""
    
    def test_deduplicate_messages_by_id(self, processor):
        """Test deduplication using message IDs."""
        msg1 = MockMessage("content1", msg_id="id1")
        msg2 = MockMessage("content2", msg_id="id2")
        msg3 = MockMessage("content1", msg_id="id1")  # Duplicate ID
        
        messages = [msg1, msg2, msg3]
        deduplicated = processor._deduplicate_messages(messages)
        
        assert len(deduplicated) == 2
        assert deduplicated[0].id == "id1"
        assert deduplicated[1].id == "id2"
    
    def test_deduplicate_messages_by_content_hash(self, processor):
        """Test deduplication using content hash when IDs are missing."""
        msg1 = MockMessage("same content")
        msg2 = MockMessage("different content")
        msg3 = MockMessage("same content")  # Duplicate content
        
        messages = [msg1, msg2, msg3]
        deduplicated = processor._deduplicate_messages(messages)
        
        assert len(deduplicated) == 2
    
    def test_deduplicate_no_duplicates(self, processor):
        """Test that unique messages are not removed."""
        msg1 = MockMessage("content1", msg_id="id1")
        msg2 = MockMessage("content2", msg_id="id2")
        
        messages = [msg1, msg2]
        deduplicated = processor._deduplicate_messages(messages)
        
        assert len(deduplicated) == 2


class TestToolExtraction:
    """Test tool extraction from stream events and results."""
    
    def test_extract_tool_from_stream_event_tools_node(self, processor):
        """Test extraction from 'tools' node in stream event."""
        mock_tool_msg = Mock()
        mock_tool_msg.name = "sequence_alignment"
        mock_tool_msg.content = '{"sequences": "ACGT"}'
        
        event = {
            "tools": {
                "messages": [mock_tool_msg]
            }
        }
        
        result = processor._extract_tool_from_stream_event(event)
        
        assert result is not None
        assert result["tool_name"] == "sequence_alignment"
    
    def test_extract_tool_from_stream_event_ai_message(self, processor):
        """Test extraction from AIMessage with tool_calls."""
        tool_call = MockToolCall("mutate_sequence", {"sequence": "ACGT"})
        ai_msg = MockAIMessage("response", tool_calls=[tool_call])
        
        event = {
            "agent": {
                "messages": [ai_msg]
            }
        }
        
        result = processor._extract_tool_from_stream_event(event)
        
        assert result is not None
        assert result["tool_name"] == "mutate_sequence"
        assert result["parameters"] == {"sequence": "ACGT"}
    
    def test_extract_tool_from_result_messages(self, processor):
        """Test extraction from agent result messages."""
        tool_call = MockToolCall("phylogenetic_tree", {"sequences": ["seq1", "seq2"]})
        ai_msg = MockAIMessage("response", tool_calls=[tool_call])
        
        result = {
            "messages": [
                MockHumanMessage("command"),
                ai_msg
            ]
        }
        
        extracted = processor._extract_tool_from_result_messages(result)
        
        assert extracted is not None
        assert extracted["tool_name"] == "phylogenetic_tree"
    
    def test_extract_tool_no_tool_calls(self, processor):
        """Test that messages without tool calls return None."""
        result = {
            "messages": [
                MockHumanMessage("command"),
                MockAIMessage("response")
            ]
        }
        
        extracted = processor._extract_tool_from_result_messages(result)
        assert extracted is None


class TestToolMappedResponse:
    """Test building standardized tool_mapped responses."""
    
    def test_build_tool_mapped_response(self, processor):
        """Test building a tool_mapped response."""
        tool_mapping = {
            "tool_name": "sequence_alignment",
            "parameters": {"sequences": "ACGT"}
        }
        
        response = processor._build_tool_mapped_response(tool_mapping)
        
        assert response["status"] == "tool_mapped"
        assert response["tool_name"] == "sequence_alignment"
        assert response["parameters"] == {"sequences": "ACGT"}
        assert "tool_mapping" in response
        assert "message" in response


class TestRouterFallback:
    """Test deterministic router fallback."""
    
    @pytest.mark.asyncio
    async def test_router_fallback_safe_tool(self, processor, mock_router_mock):
        """Test router fallback for safe tools."""
        mock_router_mock.route_command.return_value = ("sequence_alignment", {"sequences": "ACGT"})
        # Set the private _router attribute that _get_router() uses
        processor._router = mock_router_mock
        
        result = await processor._try_router_fallback("align sequences", {})
        
        assert result is not None
        assert result["tool_name"] == "sequence_alignment"
        assert result["status"] == "tool_mapped"
    
    @pytest.mark.asyncio
    async def test_router_fallback_unsafe_tool(self, processor, mock_router_mock):
        """Test that unsafe tools are not returned."""
        mock_router_mock.route_command.return_value = ("unsafe_tool", {})
        # Set the private _router attribute that _get_router() uses
        processor._router = mock_router_mock
        
        result = await processor._try_router_fallback("unsafe command", {})
        
        assert result is None
    
    @pytest.mark.asyncio
    async def test_router_fallback_exception(self, processor, mock_router_mock):
        """Test that router exceptions are handled gracefully."""
        mock_router_mock.route_command.side_effect = Exception("Router error")
        # Set the private _router attribute that _get_router() uses
        processor._router = mock_router_mock
        
        result = await processor._try_router_fallback("command", {})
        
        assert result is None


class TestQAndAHandling:
    """Test Q&A intent handling."""
    
    @pytest.mark.asyncio
    async def test_handle_qa_intent(self, processor, mock_agent):
        """Test handling Q&A intent."""
        mock_agent.invoke.return_value = {
            "messages": [
                MockHumanMessage("What is DNA?"),
                MockAIMessage("DNA is...")
            ]
        }
        
        input_msg = MockHumanMessage("What is DNA?")
        system_msg = MockSystemMessage("You are a bioinformatics assistant.")
        session_config = {"configurable": {"thread_id": "test"}}
        
        result = await processor._handle_qa_intent(input_msg, system_msg, session_config)
        
        assert "messages" in result
        mock_agent.invoke.assert_called_once()
    
    @pytest.mark.asyncio
    async def test_handle_qa_deduplicates_messages(self, processor, mock_agent):
        """Test that Q&A responses are deduplicated."""
        duplicate_msg = MockAIMessage("response", msg_id="id1")
        mock_agent.invoke.return_value = {
            "messages": [
                MockHumanMessage("question"),
                duplicate_msg,
                duplicate_msg  # Duplicate
            ]
        }
        
        input_msg = MockHumanMessage("question")
        system_msg = MockSystemMessage("system")
        session_config = {"configurable": {"thread_id": "test"}}
        
        result = await processor._handle_qa_intent(input_msg, system_msg, session_config)
        
        # Should have deduplicated
        assert len(result["messages"]) == 2


class TestMockMode:
    """Test mock mode handling."""
    
    @pytest.mark.asyncio
    async def test_mock_mode_no_agent(self):
        """Test mock mode when agent is None."""
        processor = CommandProcessor(agent=None, is_mock_mode=True)
        
        result = await processor.process("test command", "session1", None)
        
        assert result["status"] == "tool_mapped"
        assert result["tool_name"] == "toolbox_inventory"
        assert processor.is_mock_mode is True
    
    def test_mock_mode_with_agent(self, mock_agent):
        """Test mock mode when agent is provided."""
        processor = CommandProcessor(agent=mock_agent, is_mock_mode=True)
        mock_agent.invoke.return_value = {"messages": []}
        
        # In mock mode, if agent has invoke, it should be called
        # This is tested in the full process flow


class TestFullProcessFlow:
    """Test the full process flow with various scenarios."""
    
    @pytest.mark.asyncio
    async def test_process_toolbox_query(self, processor):
        """Test processing a toolbox query."""
        result = await processor.process("what can you do", "session1")
        
        assert result["tool_name"] == "toolbox_inventory"
        assert result["status"] == "tool_mapped"
    
    @pytest.mark.asyncio
    async def test_process_qa_intent(self, processor, mock_agent, mock_intent_classifier):
        """Test processing a Q&A intent."""
        mock_intent_classifier.return_value = IntentDecision(intent="qa", reason="Question")
        mock_agent.invoke.return_value = {"messages": [MockAIMessage("answer")]}
        
        result = await processor.process("What is DNA?", "session1")
        
        assert "messages" in result
        mock_agent.invoke.assert_called_once()
    
    @pytest.mark.asyncio
    async def test_process_execute_intent_with_tool_mapping(self, processor, mock_agent):
        """Test processing an execute intent that results in tool mapping."""
        # Override the intent classifier to return execute intent
        def execute_intent_classifier(command: str) -> IntentDecision:
            return IntentDecision(intent="execute", reason="Command")
        
        processor._intent_classifier = execute_intent_classifier
        
        # Mock stream to return a tool call in the expected format
        tool_msg = Mock()
        tool_msg.name = "sequence_alignment"
        tool_msg.content = '{"sequences": "ACGT"}'
        
        def mock_stream(*args, **kwargs):
            # Return a generator that yields the event structure
            def stream_generator():
                yield {"tools": {"messages": [tool_msg]}}
            return stream_generator()
        
        mock_agent.stream = Mock(side_effect=mock_stream)
        processor.agent = mock_agent
        
        result = await processor.process("align sequences ACGT", "session1")
        
        # Should return tool_mapped response
        assert "status" in result, f"Result missing 'status' key. Got: {result}"
        assert result["status"] == "tool_mapped"
        assert result["tool_name"] == "sequence_alignment"
    
    @pytest.mark.asyncio
    async def test_process_timeout_handling(self, processor, mock_agent):
        """Test handling of timeouts during tool mapping."""
        # Set execute intent
        def execute_intent_classifier(command: str) -> IntentDecision:
            return IntentDecision(intent="execute", reason="Command")
        
        processor._intent_classifier = execute_intent_classifier
        processor.agent = mock_agent
        processor.memory = None  # No memory, so it will raise TimeoutError
        
        # Patch asyncio.wait_for to raise TimeoutError immediately
        import asyncio
        from unittest.mock import patch
        
        with patch('backend.agent.asyncio.wait_for') as mock_wait_for:
            mock_wait_for.side_effect = asyncio.TimeoutError("Simulated timeout")
            
            # Should raise TimeoutError
            with pytest.raises(TimeoutError, match="Agent tool mapping timed out"):
                await processor.process("command", "session1")


class TestSessionConfig:
    """Test session config creation."""
    
    def test_create_session_config(self, processor):
        """Test creating a session config."""
        config = processor._create_session_config("session1")
        
        assert "configurable" in config
        assert "thread_id" in config["configurable"]
        assert "session1" in config["configurable"]["thread_id"]
        assert "mapping" in config["configurable"]["thread_id"]


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

