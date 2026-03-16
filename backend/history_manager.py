import json
import uuid
from datetime import datetime
from typing import Dict, Any, List, Optional
from pathlib import Path
import logging
import os
import threading
import re

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

    # -----------------------------
    # Run Ledger (iteration model)
    # -----------------------------
    #
    # We keep the legacy `history[]` + `results{}` structure for backwards compatibility,
    # but also maintain a first-class `runs[]` list and `artifacts{}` registry.
    #
    # A "run" is a single tool invocation (or agent invocation) in a session. Runs are
    # append-only and versioned via parent_run_id links.

    @staticmethod
    def _now_iso() -> str:
        return datetime.now().isoformat()

    def _ensure_run_ledger_structures(self, session_data: Dict[str, Any]) -> None:
        if not isinstance(session_data, dict):
            return
        session_data.setdefault("runs", [])
        session_data.setdefault("artifacts", {})

    @staticmethod
    def _safe_str(v: Any) -> str:
        try:
            return str(v)
        except Exception:
            return ""

    def _new_run_id(self) -> str:
        return f"run_{uuid.uuid4().hex[:12]}"

    def _new_artifact_id(self) -> str:
        return f"art_{uuid.uuid4().hex[:12]}"

    def register_artifact(
        self,
        session_id: str,
        *,
        run_id: str,
        artifact_type: str,
        title: str,
        uri: str,
        format: Optional[str] = None,
        derived_from: Optional[str] = None,
        params: Optional[Dict[str, Any]] = None,
        extra: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Register an artifact in the session-level artifact registry.

        This is metadata-only. Tools may additionally write artifact files under
        sessions/{session_id}/runs/{run_id}/... (local-only) or to S3 (cloud mode).
        """
        self._ensure_sessions_loaded()
        session_lock = self._get_session_lock(session_id)
        with session_lock:
            session = self.sessions.get(session_id)
            if not session:
                # Ensure session exists, then fetch it
                self.ensure_session_exists(session_id)
                session = self.sessions.get(session_id)
            if not session:
                raise ValueError(f"Session {session_id} not found")

            self._ensure_run_ledger_structures(session)
            artifact_id = self._new_artifact_id()
            record = {
                "artifact_id": artifact_id,
                "run_id": run_id,
                "type": artifact_type,
                "artifact_kind": artifact_type,
                "title": title,
                "uri": uri,
                "format": format,
                "derived_from": derived_from,
                "parent_artifact_ids": [derived_from] if derived_from else [],
                "source_run_id": run_id,
                "version": None,
                "analysis_branch": "main",
                "semantic_aliases": [],
                "state_tags": [],
                "params": params or {},
                "created_at": self._now_iso(),
            }
            if extra:
                record["extra"] = extra
            session["artifacts"][artifact_id] = record
            session["updated_at"] = self._now_iso()
            self._save_session(session_id)
            return record

    def list_runs(self, session_id: str) -> List[Dict[str, Any]]:
        self._ensure_sessions_loaded()
        session = self.get_session(session_id)
        if not session:
            return []
        self._ensure_run_ledger_structures(session)
        return list(session.get("runs", []) or [])

    def get_run(self, session_id: str, run_id: str) -> Optional[Dict[str, Any]]:
        self._ensure_sessions_loaded()
        session = self.get_session(session_id)
        if not session:
            return None
        self._ensure_run_ledger_structures(session)
        for r in session.get("runs", []) or []:
            if isinstance(r, dict) and r.get("run_id") == run_id:
                return r
        return None

    def list_artifacts(self, session_id: str) -> List[Dict[str, Any]]:
        self._ensure_sessions_loaded()
        session = self.get_session(session_id)
        if not session:
            return []
        self._ensure_run_ledger_structures(session)
        arts = session.get("artifacts", {}) or {}
        if not isinstance(arts, dict):
            return []
        return list(arts.values())

    def get_artifact(self, session_id: str, artifact_id: str) -> Optional[Dict[str, Any]]:
        self._ensure_sessions_loaded()
        session = self.get_session(session_id)
        if not session:
            return None
        self._ensure_run_ledger_structures(session)
        arts = session.get("artifacts", {}) or {}
        if not isinstance(arts, dict):
            return None
        art = arts.get(artifact_id)
        return art if isinstance(art, dict) else None

    def get_artifact_aliases(self, session_id: str) -> Dict[str, Dict[str, Any]]:
        from backend.artifact_resolver import build_alias_index

        session = self.get_session(session_id)
        if not session:
            return {}
        idx = build_alias_index(session)
        aliases = idx.get("aliases", {})
        return aliases if isinstance(aliases, dict) else {}

    def get_historical_states(self, session_id: str) -> List[Dict[str, Any]]:
        from backend.artifact_resolver import build_alias_index

        session = self.get_session(session_id)
        if not session:
            return []
        idx = build_alias_index(session)
        states = idx.get("historical_states", [])
        return states if isinstance(states, list) else []

    def get_lineage_edges(self, session_id: str) -> List[Dict[str, Any]]:
        artifacts = self.list_artifacts(session_id)
        edges: List[Dict[str, Any]] = []
        for art in artifacts:
            if not isinstance(art, dict):
                continue
            child = art.get("artifact_id")
            parents = art.get("parent_artifact_ids")
            if not isinstance(parents, list):
                parents = [art.get("derived_from")] if art.get("derived_from") else []
            for p in parents:
                if p:
                    edges.append({"from_artifact_id": p, "to_artifact_id": child})
        return edges

    def resolve_run_reference(self, session_id: str, ref: str) -> Optional[Dict[str, Any]]:
        """
        Resolve user-facing run references:
        - "first" / "latest" / "last"
        - "run_..." id
        - "iteration:1" / "#1" / "1"
        """
        runs = self.list_runs(session_id)
        if not runs:
            return None
        if not ref:
            return runs[-1]

        r = ref.strip().lower()
        if r in {"latest", "last", "current"}:
            return runs[-1]
        if r in {"first", "initial"}:
            return runs[0]
        if r.startswith("run_"):
            return self.get_run(session_id, ref)
        m = re.match(r"^(?:iteration\s*[:#]\s*)?(\d+)$", r)
        if m:
            idx = int(m.group(1))
            for run in runs:
                if isinstance(run, dict) and run.get("iteration_index") == idx:
                    return run
        m2 = re.match(r"^#(\d+)$", r)
        if m2:
            idx = int(m2.group(1))
            for run in runs:
                if isinstance(run, dict) and run.get("iteration_index") == idx:
                    return run
        return None

    def add_run(
        self,
        session_id: str,
        *,
        command: str,
        tool: str,
        result: Dict[str, Any],
        run_id: Optional[str] = None,
        tool_args: Optional[Dict[str, Any]] = None,
        inputs: Optional[List[Dict[str, Any]]] = None,
        outputs: Optional[List[Dict[str, Any]]] = None,
        produced_artifacts: Optional[List[Dict[str, Any]]] = None,
        parent_run_id: Optional[str] = None,
        metadata: Optional[Dict[str, Any]] = None,
    ) -> Dict[str, Any]:
        """
        Append a run record to session['runs'] and return the run record.
        """
        self._ensure_sessions_loaded()
        session_lock = self._get_session_lock(session_id)
        with session_lock:
            if session_id not in self.sessions:
                self.ensure_session_exists(session_id)
            session = self.sessions.get(session_id)
            if not session:
                raise ValueError(f"Session {session_id} not found")

            self._ensure_run_ledger_structures(session)

            run_id = run_id or self._new_run_id()
            iteration_index = len(session.get("runs", []) or []) + 1

            # Register produced artifacts (if any) into the session-level registry.
            # We accept a loose dict shape and normalize to:
            # {artifact_id, run_id, type, title, uri, format?, derived_from?, params?}
            artifacts_in: List[Dict[str, Any]] = []
            if isinstance(produced_artifacts, list):
                artifacts_in = [a for a in produced_artifacts if isinstance(a, dict)]

            enriched_artifacts: List[Dict[str, Any]] = []
            for a in artifacts_in:
                a_type = a.get("type") or a.get("artifact_type") or "unknown"
                a_title = a.get("title") or a.get("label") or a.get("name") or "artifact"
                a_uri = a.get("uri") or a.get("url") or ""
                if not a_uri:
                    # Skip malformed artifact records
                    continue
                a_format = a.get("format") or a.get("kind")
                a_derived = a.get("derived_from") or a.get("derived_from_artifact_id")
                a_params = a.get("params") if isinstance(a.get("params"), dict) else {}
                a_extra = a.get("extra") if isinstance(a.get("extra"), dict) else None

                artifact_id = a.get("artifact_id") or self._new_artifact_id()
                record = {
                    "artifact_id": artifact_id,
                    "run_id": run_id,
                    "type": a_type,
                    "artifact_kind": a.get("artifact_kind") or a_type,
                    "title": a_title,
                    "uri": a_uri,
                    "format": a_format,
                    "derived_from": a_derived,
                    "parent_artifact_ids": (
                        a.get("parent_artifact_ids")
                        if isinstance(a.get("parent_artifact_ids"), list)
                        else ([a_derived] if a_derived else [])
                    ),
                    "source_run_id": a.get("source_run_id") or run_id,
                    "version": a.get("version"),
                    "analysis_branch": a.get("analysis_branch") or "main",
                    "semantic_aliases": (
                        a.get("semantic_aliases")
                        if isinstance(a.get("semantic_aliases"), list)
                        else []
                    ),
                    "state_tags": (
                        a.get("state_tags")
                        if isinstance(a.get("state_tags"), list)
                        else []
                    ),
                    "params": a_params,
                    "created_at": self._now_iso(),
                }
                if a_extra:
                    record["extra"] = a_extra

                session["artifacts"][artifact_id] = record
                enriched_artifacts.append(record)

            run = {
                "run_id": run_id,
                "iteration_index": iteration_index,
                "timestamp": self._now_iso(),
                "command": command,
                "tool": tool,
                "tool_args": tool_args or {},
                "inputs": inputs or [],
                "outputs": outputs or [],
                "produced_artifacts": enriched_artifacts,
                "parent_run_id": parent_run_id,
                "result_key": None,  # filled by add_history_entry for legacy results map
                "metadata": metadata or {},
            }

            session["runs"].append(run)
            session["updated_at"] = self._now_iso()
            self._save_session(session_id)
            return run
    
    def _get_s3_client(self):
        """Get or create S3 client. Returns None if boto3 is not available or AWS credentials are missing."""
        if os.getenv("HELIX_MOCK_MODE") == "1":
            return None
        if self._s3_client is not None:
            return self._s3_client
        
        try:
            import boto3
            from botocore.config import Config as _BotoConfig
            # Short timeouts so blocking S3 calls don't stall the asyncio event loop.
            _cfg = _BotoConfig(connect_timeout=3, read_timeout=5, retries={"max_attempts": 1})
            self._s3_client = boto3.client('s3', config=_cfg)
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
        """Mark sessions as ready without scanning all files.

        Loading every session JSON at startup is O(n_sessions) and can block
        the event loop for minutes when thousands of session files accumulate.
        Instead, sessions are loaded on-demand in ensure_session_exists via
        _load_single_session_from_disk.
        """
        self._sessions_loaded = True

    def _load_single_session_from_disk(self, session_id: str) -> bool:
        """Attempt to load a single session JSON file into self.sessions.

        Returns True if the file existed and was loaded successfully.
        """
        session_file = self.storage_dir / f"{session_id}.json"
        if not session_file.exists():
            return False
        try:
            with open(session_file, "r") as f:
                self.sessions[session_id] = json.load(f)
            return True
        except json.JSONDecodeError as e:
            logger.error(f"Corrupted session file {session_file}: {e}")
            backup = session_file.with_suffix(".json.corrupted")
            try:
                session_file.rename(backup)
            except Exception:
                pass
            return False
        except Exception as e:
            logger.error(f"Failed to load session {session_file}: {e}")
            return False

    def _load_existing_sessions(self):
        """Load existing session data from disk (used by tests / maintenance tasks only)."""
        for session_file in self.storage_dir.glob("*.json"):
            if session_file.name.endswith('.tmp'):
                continue
            try:
                with open(session_file, 'r') as f:
                    session_data = json.load(f)
                    session_id = session_file.stem
                    self.sessions[session_id] = session_data
            except json.JSONDecodeError as e:
                logger.error(f"Failed to load session {session_file} due to JSON corruption: {e}")
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
                "runs": [],
                "artifacts": {},
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

    def get_session_storage_paths(self, session_id: str) -> Dict[str, str]:
        """Return absolute paths where this session's inputs, prompts, and results are stored.

        - storage_root: base directory (e.g. sessions/)
        - session_file: session state JSON (history, results, prompts)
        - session_dir: per-session directory (artifacts, runs)
        """
        root = self.storage_dir.resolve()
        return {
            "storage_root": str(root),
            "session_file": str(root / f"{session_id}.json"),
            "session_dir": str(root / session_id),
        }

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
                    "runs": [],
                    "artifacts": {},
                    "metadata": {
                        "s3_path": s3_path,
                        "s3_bucket": self.s3_bucket_name if s3_path else None,
                        "local_path": str(session_dir)
                    }
                }
                self._save_session(session_id)

            # Ensure run-ledger structures exist (back-compat for older session files)
            self._ensure_run_ledger_structures(self.sessions[session_id])
            
            # Serialize the result to handle LangChain message objects
            serialized_result = serialize_langchain_messages(result)

            # Build run record (iteration-aware).
            tool_args = None
            inputs = None
            outputs = None
            produced_artifacts = None
            parent_run_id = None
            pre_run_id = None
            run_metadata = {}
            if isinstance(metadata, dict):
                tool_args = metadata.get("tool_args")
                inputs = metadata.get("inputs")
                outputs = metadata.get("outputs")
                produced_artifacts = metadata.get("produced_artifacts") or metadata.get("artifacts")
                parent_run_id = metadata.get("parent_run_id")
                pre_run_id = metadata.get("run_id")
                # store the rest
                run_metadata = {k: v for k, v in metadata.items() if k not in {
                    "tool_args", "inputs", "outputs", "produced_artifacts", "artifacts", "parent_run_id", "run_id"
                }}

            run_record = self.add_run(
                session_id,
                command=command,
                tool=tool,
                result=serialized_result if isinstance(serialized_result, dict) else {"value": serialized_result},
                run_id=self._safe_str(pre_run_id) if pre_run_id else None,
                tool_args=tool_args if isinstance(tool_args, dict) else None,
                inputs=inputs if isinstance(inputs, list) else None,
                outputs=outputs if isinstance(outputs, list) else None,
                produced_artifacts=produced_artifacts if isinstance(produced_artifacts, list) else None,
                parent_run_id=self._safe_str(parent_run_id) if parent_run_id else None,
                metadata=run_metadata if isinstance(run_metadata, dict) else None,
            )

            # If the run produced artifacts, inject "results_viewer" payload into the *live* result dict
            # so /execute can render plots/links immediately without reading the session file.
            try:
                if isinstance(result, dict):
                    arts = run_record.get("produced_artifacts") or []
                    if isinstance(arts, list) and len(arts) > 0:
                        links_in = result.get("links") if isinstance(result.get("links"), list) else []
                        visuals_in = result.get("visuals") if isinstance(result.get("visuals"), list) else []

                        def _seen_urls(items):
                            seen = set()
                            for it in items:
                                if isinstance(it, dict) and isinstance(it.get("url"), str):
                                    seen.add(it["url"])
                            return seen

                        seen_link_urls = _seen_urls(links_in)
                        seen_visual_urls = _seen_urls(visuals_in)

                        new_links = []
                        new_visuals = []
                        for a in arts:
                            if not isinstance(a, dict):
                                continue
                            artifact_id = a.get("artifact_id")
                            if not artifact_id:
                                continue
                            url = f"/session/{session_id}/artifacts/{artifact_id}/download"
                            title = a.get("title") or a.get("name") or "artifact"
                            a_type = a.get("type") or "artifact"
                            a_format = a.get("format")
                            if url not in seen_link_urls:
                                new_links.append(
                                    {
                                        "label": title,
                                        "url": url,
                                        "artifact_id": artifact_id,
                                        "type": a_type,
                                        "format": a_format,
                                    }
                                )
                                seen_link_urls.add(url)

                            # Heuristic: treat certain artifact types/formats as images for UI preview.
                            fmt = (a_format or "").lower() if isinstance(a_format, str) else ""
                            uri = (a.get("uri") or "")
                            is_image = (
                                a_type in {"plot", "image"}
                                and (fmt in {"png", "jpg", "jpeg", "gif", "webp", "svg"} or str(uri).lower().endswith((".png", ".jpg", ".jpeg", ".gif", ".webp", ".svg")))
                            )
                            if is_image and url not in seen_visual_urls:
                                new_visuals.append({"type": "image", "url": url, "title": title})
                                seen_visual_urls.add(url)

                        if new_links:
                            result["links"] = links_in + new_links
                        if new_visuals:
                            result["visuals"] = visuals_in + new_visuals

                        # Hint frontend to render a results viewer if we have visuals or links
                        if (new_links or new_visuals) and result.get("visualization_type") != "results_viewer":
                            result["visualization_type"] = "results_viewer"
            except Exception:
                pass
            
            entry = {
                "timestamp": datetime.now().isoformat(),
                "command": command,
                "tool": tool,
                "result": serialized_result,
                "metadata": metadata or {},
                "run_id": run_record.get("run_id"),
                "iteration_index": run_record.get("iteration_index"),
            }
            
            self.sessions[session_id]["history"].append(entry)
            self.sessions[session_id]["updated_at"] = datetime.now().isoformat()
            
            # Store result with a unique key for later reference
            result_key = f"{tool}_{len(self.sessions[session_id]['history'])}"
            self.sessions[session_id]["results"][result_key] = serialized_result

            # Link legacy result key back into the run record (best-effort)
            try:
                runs = self.sessions[session_id].get("runs", [])
                if isinstance(runs, list) and runs:
                    runs[-1]["result_key"] = result_key
            except Exception:
                pass
            
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
                # Try to restore from disk before creating a fresh session
                if self._load_single_session_from_disk(session_id):
                    return  # Session was already on disk – nothing more to do

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

    # ── Workflow Checkpoint API ───────────────────────────────────────────────

    def save_checkpoint(self, session_id: str, checkpoint: "WorkflowCheckpoint") -> None:  # noqa: F821
        """Persist a WorkflowCheckpoint into the session store and flush to disk."""
        from backend.workflow_checkpoint import WorkflowCheckpoint  # local import to avoid circular
        self._ensure_sessions_loaded()
        session_lock = self._get_session_lock(session_id)
        with session_lock:
            if session_id not in self.sessions:
                self.ensure_session_exists(session_id)
            self.sessions[session_id]["__checkpoint__"] = checkpoint.to_dict()
            self.sessions[session_id]["updated_at"] = datetime.now().isoformat()
        self._save_session(session_id)

    def load_checkpoint(self, session_id: str) -> "WorkflowCheckpoint":  # noqa: F821
        """Load the current WorkflowCheckpoint for a session.

        Returns an IDLE checkpoint if none exists.
        """
        from backend.workflow_checkpoint import WorkflowCheckpoint
        self._ensure_sessions_loaded()
        session = self.sessions.get(session_id) or {}
        raw = session.get("__checkpoint__")
        if raw and isinstance(raw, dict):
            try:
                return WorkflowCheckpoint.from_dict(raw)
            except Exception:
                pass
        return WorkflowCheckpoint.idle()

    def clear_checkpoint(self, session_id: str) -> None:
        """Remove the checkpoint and return session to IDLE."""
        from backend.workflow_checkpoint import WorkflowCheckpoint
        self.save_checkpoint(session_id, WorkflowCheckpoint.idle())


# Global history manager instance
history_manager = HistoryManager() 
