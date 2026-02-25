import json
import uuid
from datetime import datetime
from typing import Dict, Any, List, Optional
from pathlib import Path
import logging
import os
import threading

# Try to import fcntl for file locking (not available on Windows)
try:
    import fcntl
    HAS_FCNTL = True
except ImportError:
    HAS_FCNTL = False

logger = logging.getLogger(__name__)

# S3 Configuration
S3_BUCKET_NAME = os.getenv("S3_BUCKET_NAME", "helix-ai-frontend-794270057041-us-west-1")

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

def serialize_for_json(obj: Any, _seen: Optional[frozenset] = None) -> Any:
    """Convert objects to JSON-serializable format.

    Handles numpy types, LangChain message objects, and — critically — circular
    references.  Without cycle detection, LangGraph/LangChain message objects
    (whose `additional_kwargs` / `response_metadata` can contain self-referencing
    Pydantic model dicts) cause infinite recursion and a RecursionError.
    """
    import numpy as np

    if _seen is None:
        _seen = frozenset()

    if isinstance(obj, dict):
        obj_id = id(obj)
        if obj_id in _seen:
            # Break the cycle — return a sentinel rather than recursing forever
            return "[circular reference]"
        _seen = _seen | {obj_id}
        return {key: serialize_for_json(value, _seen) for key, value in obj.items()}
    elif isinstance(obj, list):
        obj_id = id(obj)
        if obj_id in _seen:
            return "[circular reference]"
        _seen = _seen | {obj_id}
        return [serialize_for_json(item, _seen) for item in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (int, float, str, bool)) or obj is None:
        return obj
    else:
        # Unknown non-serializable type — convert to string rather than letting
        # json.dumps fail later with a TypeError.
        try:
            import json
            json.dumps(obj)   # test serialisability cheaply
            return obj
        except (TypeError, ValueError):
            return str(obj)

class HistoryManager:
    """Manages user session history and results for bioinformatics operations."""
    
    def __init__(self, storage_dir: str = "sessions"):
        self.storage_dir = Path(storage_dir)
        self.storage_dir.mkdir(exist_ok=True)
        self.sessions: Dict[str, Dict[str, Any]] = {}
        self.s3_bucket_name = S3_BUCKET_NAME
        self._s3_client = None
        self._sessions_loaded = False  # Lazy loading flag
        # Don't load sessions at import time - load lazily on first access
        # This speeds up server startup when there are many session files
        # Locks for serializing access to each session file (using RLock for reentrancy)
        self._session_locks: Dict[str, threading.RLock] = {}
        self._locks_lock = threading.Lock()  # Lock for accessing the _session_locks dict
    
    def _get_s3_client(self):
        """Get or create S3 client. Returns None if boto3 is not available or AWS credentials are missing."""
        if self._s3_client is not None:
            return self._s3_client
        
        try:
            import boto3
            # Try to create S3 client - will use default AWS credentials from environment/IAM role
            self._s3_client = boto3.client('s3')
            return self._s3_client
        except ImportError:
            logger.warning("boto3 not available - S3 path creation will be skipped")
            return None
        except Exception as e:
            logger.warning(f"Failed to initialize S3 client: {e} - S3 path creation will be skipped")
            return None
    
    def _create_s3_session_path(self, session_id: str) -> Optional[str]:
        """Create an S3 path prefix for a session by uploading a marker file.
        
        Returns the S3 key path if successful, None otherwise.
        The path structure is: {session_id}/.session_initialized
        """
        s3_client = self._get_s3_client()
        if s3_client is None:
            logger.debug(f"S3 client not available - skipping S3 path creation for session {session_id}")
            return None
        
        try:
            # Create a marker file to initialize the session path in S3
            # Path structure: {session_id}/.session_initialized
            s3_key = f"{session_id}/.session_initialized"
            
            # Upload an empty marker file to create the path
            s3_client.put_object(
                Bucket=self.s3_bucket_name,
                Key=s3_key,
                Body=b"",
                Metadata={
                    "session_id": session_id,
                    "created_at": datetime.now().isoformat(),
                    "purpose": "session_path_initialization"
                }
            )
            
            logger.info(f"Created S3 path for session {session_id}: s3://{self.s3_bucket_name}/{session_id}/")
            return f"{session_id}/"
        except Exception as e:
            logger.error(f"Failed to create S3 path for session {session_id}: {e}")
            return None
    
    def _ensure_sessions_loaded(self):
        """Lazily load existing session data from disk on first access."""
        if self._sessions_loaded:
            return
        self._load_existing_sessions()
        self._sessions_loaded = True
    
    def _load_existing_sessions(self):
        """Load existing session data from disk."""
        for session_file in self.storage_dir.glob("*.json"):
            # Skip temporary files
            if session_file.name.endswith('.tmp'):
                continue
            
            try:
                with open(session_file, 'r') as f:
                    session_data = json.load(f)
                    session_id = session_file.stem
                    self.sessions[session_id] = session_data
            except json.JSONDecodeError as e:
                # Handle corrupted JSON files
                logger.error(f"Failed to load session {session_file} due to JSON corruption: {e}")
                # Try to create a backup of the corrupted file
                backup_file = session_file.with_suffix('.json.corrupted')
                try:
                    session_file.rename(backup_file)
                    logger.warning(f"Moved corrupted session file to {backup_file}")
                except Exception as backup_error:
                    logger.error(f"Failed to backup corrupted file {session_file}: {backup_error}")
            except Exception as e:
                logger.error(f"Failed to load session {session_file}: {e}")
    
    def create_session(self, user_id: Optional[str] = None) -> str:
        """Create a new session and return session ID.
        
        Also creates:
        - A local directory structure: sessions/{session_id}/
        - An S3 path for the session in the configured bucket.
        
        Thread-safe: Uses locks to prevent concurrent modifications.
        """
        session_id = str(uuid.uuid4())
        
        # Get lock for the new session (will create the lock)
        session_lock = self._get_session_lock(session_id)
        
        with session_lock:
            # Create local directory for the session
            session_dir = self.storage_dir / session_id
            session_dir.mkdir(exist_ok=True)
            logger.info(f"Created session directory: {session_dir}")
            
            # Create S3 path for the session
            s3_path = self._create_s3_session_path(session_id)
            
            session_data = {
                "session_id": session_id,
                "user_id": user_id,
                "created_at": datetime.now().isoformat(),
                "updated_at": datetime.now().isoformat(),
                "history": [],
                "results": {},
                "metadata": {
                    "s3_path": s3_path,
                    "s3_bucket": self.s3_bucket_name if s3_path else None,
                    "local_path": str(session_dir)
                }
            }
            self.sessions[session_id] = session_data
            self._save_session(session_id)
            logger.info(f"Created new session: {session_id}" + (f" with S3 path: {s3_path}" if s3_path else " (S3 path creation skipped)"))
            return session_id
    
    def get_session(self, session_id: str) -> Optional[Dict[str, Any]]:
        """Get session data by ID.
        
        If session doesn't exist in memory but directory exists, load it.
        Also ensures the session directory exists.
        """
        self._ensure_sessions_loaded()
        if session_id in self.sessions:
            session = self.sessions.get(session_id)
            # Ensure directory exists even if session is in memory
            self._ensure_session_directory(session_id, session)
            return session
        
        # Try to load from disk if file exists
        session_file = self.storage_dir / f"{session_id}.json"
        if session_file.exists():
            try:
                with open(session_file, 'r') as f:
                    session_data = json.load(f)
                    # Ensure local directory exists
                    self._ensure_session_directory(session_id, session_data)
                    self.sessions[session_id] = session_data
                    return session_data
            except json.JSONDecodeError as e:
                # Handle corrupted JSON files
                logger.error(f"Failed to load session {session_file} due to JSON corruption: {e}")
                # Try to create a backup of the corrupted file
                backup_file = session_file.with_suffix('.json.corrupted')
                try:
                    session_file.rename(backup_file)
                    logger.warning(f"Moved corrupted session file to {backup_file}")
                except Exception as backup_error:
                    logger.error(f"Failed to backup corrupted file {session_file}: {backup_error}")
                # Return None to allow session recreation
                return None
            except Exception as e:
                logger.error(f"Failed to load session {session_id} from disk: {e}")
        
        return None
    
    def _ensure_session_directory(self, session_id: str, session_data: Optional[Dict[str, Any]] = None) -> None:
        """Ensure the session directory exists and update metadata if needed.
        
        Thread-safe: Uses locks to prevent concurrent modifications.
        """
        session_lock = self._get_session_lock(session_id)
        
        with session_lock:
            session_dir = self.storage_dir / session_id
            if not session_dir.exists():
                session_dir.mkdir(parents=True, exist_ok=True)
                logger.info(f"Created missing session directory: {session_dir}")
            
            # Update metadata if needed
            if session_data:
                metadata = session_data.setdefault("metadata", {})
                if not metadata.get("local_path"):
                    metadata["local_path"] = str(session_dir)
                    # Save updated session data
                    self.sessions[session_id] = session_data
                    self._save_session(session_id)
    
    def add_history_entry(self, session_id: str, command: str, tool: str, 
                         result: Dict[str, Any], metadata: Optional[Dict[str, Any]] = None):
        """Add a new history entry to a session.
        
        Thread-safe: Uses locks to prevent concurrent modifications.
        """
        self._ensure_sessions_loaded()
        import time
        start_time = time.time()
        
        # Use per-session lock to serialize access
        session_lock = self._get_session_lock(session_id)
        
        with session_lock:
            if session_id not in self.sessions:
                # Auto-create session if missing to avoid hard failures
                # Create local directory for the session
                session_dir = self.storage_dir / session_id
                session_dir.mkdir(exist_ok=True)
                logger.info(f"Auto-created session directory: {session_dir}")
                
                # Also create S3 path for auto-created sessions
                s3_path = self._create_s3_session_path(session_id)
                self.sessions[session_id] = {
                    "session_id": session_id,
                    "user_id": None,
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat(),
                    "history": [],
                    "results": {},
                    "metadata": {
                        "s3_path": s3_path,
                        "s3_bucket": self.s3_bucket_name if s3_path else None,
                        "local_path": str(session_dir)
                    }
                }
                self._save_session(session_id)
            
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
            
            save_start = time.time()
            self._save_session(session_id)
            save_duration = time.time() - save_start
            total_duration = time.time() - start_time
            print(f"✅ [PERF] add_history_entry took {total_duration*1000:.2f}ms (save: {save_duration*1000:.2f}ms)")
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
        self._ensure_sessions_loaded()
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
    
    def _get_session_lock(self, session_id: str) -> threading.RLock:
        """Get or create a reentrant lock for a specific session."""
        with self._locks_lock:
            if session_id not in self._session_locks:
                self._session_locks[session_id] = threading.RLock()
            return self._session_locks[session_id]
    
    def ensure_session_exists(self, session_id: str) -> None:
        """Ensure a session exists, creating it if necessary.
        
        Thread-safe: Uses locks to prevent concurrent modifications.
        """
        self._ensure_sessions_loaded()
        session_lock = self._get_session_lock(session_id)
        
        with session_lock:
            if session_id not in self.sessions:
                # Auto-create session if missing
                session_dir = self.storage_dir / session_id
                session_dir.mkdir(exist_ok=True)
                logger.info(f"Auto-created session directory: {session_dir}")
                
                # Also create S3 path for auto-created sessions
                s3_path = self._create_s3_session_path(session_id)
                self.sessions[session_id] = {
                    "session_id": session_id,
                    "user_id": None,
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat(),
                    "history": [],
                    "results": {},
                    "metadata": {
                        "s3_path": s3_path,
                        "s3_bucket": self.s3_bucket_name if s3_path else None,
                        "local_path": str(session_dir)
                    }
                }
                self._save_session(session_id)
    
    def _save_session(self, session_id: str):
        """Save session data to disk using atomic writes and file locking to prevent corruption.
        
        Thread-safe: Uses per-session locks to serialize access.
        """
        import time
        save_start = time.time()
        if session_id not in self.sessions:
            raise ValueError(f"Session {session_id} not found")
        
        # Use per-session lock to serialize access (handles in-process concurrency)
        session_lock = self._get_session_lock(session_id)
        
        with session_lock:
            session_file = self.storage_dir / f"{session_id}.json"
            temp_file = self.storage_dir / f"{session_id}.json.tmp"
            
            try:
                # Serialize the session data to handle numpy types and other non-JSON-serializable objects
                serialize_start = time.time()
                serialized_session = serialize_for_json(self.sessions[session_id])
                serialize_duration = time.time() - serialize_start
                
                # Write to a temporary file first (atomic write pattern)
                write_start = time.time()
                
                # Use file locking for the temp file write to prevent concurrent writes
                with open(temp_file, "w") as f:
                    # Acquire exclusive lock on the temp file if available
                    if HAS_FCNTL:
                        try:
                            fcntl.flock(f.fileno(), fcntl.LOCK_EX)
                        except (IOError, OSError) as lock_error:
                            # Lock acquisition failed, but continue with write
                            # Thread lock should handle most concurrency issues
                            logger.debug(f"Failed to acquire file lock for {session_id}: {lock_error}")
                    
                    json.dump(serialized_session, f, indent=2)
                    # Ensure data is written to disk
                    f.flush()
                    os.fsync(f.fileno())
                    
                    # Lock is automatically released when file is closed
                
                # Atomically replace the old file with the new one
                # This is atomic on most filesystems (including macOS and Linux)
                temp_file.replace(session_file)
                
                write_duration = time.time() - write_start
                total_duration = time.time() - save_start
                print(f"✅ [PERF] _save_session took {total_duration*1000:.2f}ms (serialize: {serialize_duration*1000:.2f}ms, write: {write_duration*1000:.2f}ms)")
            except Exception as e:
                logger.error(f"Failed to save session {session_id}: {e}", exc_info=True)
                # Clean up temporary file if it exists
                if temp_file.exists():
                    try:
                        temp_file.unlink()
                    except Exception as cleanup_error:
                        logger.warning(f"Failed to clean up temporary file {temp_file}: {cleanup_error}")
                raise
    
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
    
    def add_dataset_reference(self, session_id: str, dataset_ref: Dict[str, Any]) -> bool:
        """Add a dataset reference to session metadata.
        
        Args:
            session_id: Session ID
            dataset_ref: Dictionary with s3_bucket, s3_key, filename, size, dataset_id
        
        Returns:
            True if successful, False otherwise
            
        Thread-safe: Uses locks to prevent concurrent modifications.
        """
        session_lock = self._get_session_lock(session_id)
        
        with session_lock:
            if session_id not in self.sessions:
                logger.warning(f"Session {session_id} not found, cannot add dataset reference")
                return False
            
            session = self.sessions[session_id]
            if "metadata" not in session:
                session["metadata"] = {}
            if "dataset_references" not in session["metadata"]:
                session["metadata"]["dataset_references"] = []
            
            # Check if reference already exists
            existing_refs = session["metadata"]["dataset_references"]
            if any(ref.get("s3_key") == dataset_ref.get("s3_key") for ref in existing_refs):
                logger.info(f"Dataset reference already exists: {dataset_ref.get('s3_key')}")
                return True
            
            # Add timestamp
            dataset_ref["linked_at"] = datetime.now().isoformat()
            session["metadata"]["dataset_references"].append(dataset_ref)
            session["updated_at"] = datetime.now().isoformat()
            
            self._save_session(session_id)
            logger.info(f"Added dataset reference to session {session_id}: {dataset_ref.get('s3_key')}")
            return True
    
    def get_dataset_references(self, session_id: str) -> List[Dict[str, Any]]:
        """Get all dataset references for a session."""
        if session_id not in self.sessions:
            return []
        
        return self.sessions[session_id].get("metadata", {}).get("dataset_references", [])
    
    def get_session_files(self, session_id: str) -> List[Dict[str, Any]]:
        """Get all files available to a session (both uploaded and referenced).
        
        Returns a list of file dictionaries with:
        - type: "uploaded" or "dataset_reference"
        - s3_bucket: S3 bucket name
        - s3_key: S3 object key
        - filename: Original filename
        - size: File size in bytes
        - dataset_id: (optional) Dataset ID for references
        """
        self._ensure_sessions_loaded()
        if session_id not in self.sessions:
            return []
        
        files = []
        metadata = self.sessions[session_id].get("metadata", {})
        
        # Add small uploaded files
        for file_info in metadata.get("uploaded_files", []):
            files.append({
                "type": "uploaded",
                "s3_bucket": file_info.get("s3_bucket", self.s3_bucket_name),
                "s3_key": file_info.get("s3_key", ""),
                "filename": file_info.get("filename", ""),
                "size": file_info.get("size", 0)
            })
        
        # Add large dataset references
        for ref in metadata.get("dataset_references", []):
            files.append({
                "type": "dataset_reference",
                "s3_bucket": ref.get("s3_bucket", ""),
                "s3_key": ref.get("s3_key", ""),
                "filename": ref.get("filename", ""),
                "size": ref.get("size", 0),
                "dataset_id": ref.get("dataset_id", "")
            })
        
        return files

# Global history manager instance
history_manager = HistoryManager() 
