import json
import uuid
from datetime import datetime
from typing import Dict, Any, List, Optional
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def serialize_langchain_messages(result: Any) -> Any:
    """Convert LangChain message objects to serializable format."""
    if isinstance(result, dict):
        if "messages" in result:
            # Convert LangChain messages to simple dict format
            serialized_messages = []
            for msg in result["messages"]:
                if hasattr(msg, 'content'):
                    serialized_msg = {
                        "type": getattr(msg, 'type', 'unknown'),
                        "content": getattr(msg, 'content', str(msg)),
                        "name": getattr(msg, 'name', None),
                        "id": getattr(msg, 'id', None)
                    }
                    # Add additional kwargs if they exist
                    if hasattr(msg, 'additional_kwargs'):
                        serialized_msg["additional_kwargs"] = msg.additional_kwargs
                    if hasattr(msg, 'response_metadata'):
                        serialized_msg["response_metadata"] = msg.response_metadata
                    serialized_messages.append(serialized_msg)
                else:
                    # Fallback for non-message objects
                    serialized_messages.append(str(msg))
            result["messages"] = serialized_messages
        return result
    elif hasattr(result, 'content'):
        # Single message object
        return {
            "type": getattr(result, 'type', 'unknown'),
            "content": getattr(result, 'content', str(result)),
            "name": getattr(result, 'name', None),
            "id": getattr(result, 'id', None)
        }
    else:
        return result

class HistoryManager:
    """Manages user session history and results for bioinformatics operations."""
    
    def __init__(self, storage_dir: str = "sessions"):
        self.storage_dir = Path(storage_dir)
        self.storage_dir.mkdir(exist_ok=True)
        self.sessions: Dict[str, Dict[str, Any]] = {}
        self._load_existing_sessions()
    
    def _load_existing_sessions(self):
        """Load existing session data from disk."""
        for session_file in self.storage_dir.glob("*.json"):
            try:
                with open(session_file, 'r') as f:
                    session_data = json.load(f)
                    session_id = session_file.stem
                    self.sessions[session_id] = session_data
            except Exception as e:
                logger.error(f"Failed to load session {session_file}: {e}")
    
    def create_session(self, user_id: Optional[str] = None) -> str:
        """Create a new session and return session ID."""
        session_id = str(uuid.uuid4())
        session_data = {
            "session_id": session_id,
            "user_id": user_id,
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat(),
            "history": [],
            "results": {},
            "metadata": {}
        }
        self.sessions[session_id] = session_data
        self._save_session(session_id)
        logger.info(f"Created new session: {session_id}")
        return session_id
    
    def get_session(self, session_id: str) -> Optional[Dict[str, Any]]:
        """Get session data by ID."""
        return self.sessions.get(session_id)
    
    def add_history_entry(self, session_id: str, command: str, tool: str, 
                         result: Dict[str, Any], metadata: Optional[Dict[str, Any]] = None):
        """Add a new history entry to a session."""
        if session_id not in self.sessions:
            raise ValueError(f"Session {session_id} not found")
        
        # Serialize the result to handle LangChain message objects
        serialized_result = serialize_langchain_messages(result)
        
        entry = {
            "timestamp": datetime.now().isoformat(),
            "command": command,
            "tool": tool,
            "result": serialized_result,
            "metadata": metadata or {}
        }
        
        self.sessions[session_id]["history"].append(entry)
        self.sessions[session_id]["updated_at"] = datetime.now().isoformat()
        
        # Store result with a unique key for later reference
        result_key = f"{tool}_{len(self.sessions[session_id]['history'])}"
        self.sessions[session_id]["results"][result_key] = serialized_result
        
        self._save_session(session_id)
        logger.info(f"Added history entry to session {session_id}: {tool}")
    
    def get_latest_result(self, session_id: str, tool: str) -> Optional[Dict[str, Any]]:
        """Get the most recent result for a specific tool in a session."""
        if session_id not in self.sessions:
            return None
        
        # Find the most recent entry for this tool
        for entry in reversed(self.sessions[session_id]["history"]):
            if entry["tool"] == tool:
                return entry["result"]
        return None
    
    def get_all_results(self, session_id: str, tool: str) -> List[Dict[str, Any]]:
        """Get all results for a specific tool in a session."""
        if session_id not in self.sessions:
            return []
        
        results = []
        for entry in self.sessions[session_id]["history"]:
            if entry["tool"] == tool:
                results.append(entry["result"])
        return results
    
    def get_session_summary(self, session_id: str) -> Dict[str, Any]:
        """Get a summary of session activity."""
        if session_id not in self.sessions:
            return {}
        
        session = self.sessions[session_id]
        tool_counts = {}
        for entry in session["history"]:
            tool = entry["tool"]
            tool_counts[tool] = tool_counts.get(tool, 0) + 1
        
        return {
            "session_id": session_id,
            "user_id": session.get("user_id"),
            "created_at": session["created_at"],
            "updated_at": session["updated_at"],
            "total_operations": len(session["history"]),
            "tool_usage": tool_counts,
            "available_results": list(session["results"].keys())
        }
    
    def _save_session(self, session_id: str):
        """Save session data to disk."""
        if session_id not in self.sessions:
            raise ValueError(f"Session {session_id} not found")
        session_file = self.storage_dir / f"{session_id}.json"
        try:
            # Only write JSON-serializable data
            import json
            with open(session_file, "w") as f:
                json.dump(self.sessions[session_id], f, indent=2)
        except Exception as e:
            print(f"[ERROR] Failed to save session {session_id}: {e}")
            # Optionally, remove the corrupted file
            # import os
            # if os.path.exists(session_file):
            #     os.remove(session_file)
    
    def cleanup_old_sessions(self, max_age_days: int = 30):
        """Clean up old sessions."""
        cutoff = datetime.now().timestamp() - (max_age_days * 24 * 3600)
        to_delete = []
        
        for session_id, session_data in self.sessions.items():
            try:
                created_at = datetime.fromisoformat(session_data["created_at"]).timestamp()
                if created_at < cutoff:
                    to_delete.append(session_id)
            except:
                pass
        
        for session_id in to_delete:
            self.sessions.pop(session_id, None)
            session_file = self.storage_dir / f"{session_id}.json"
            if session_file.exists():
                session_file.unlink()
        
        logger.info(f"Cleaned up {len(to_delete)} old sessions")

# Global history manager instance
history_manager = HistoryManager() 