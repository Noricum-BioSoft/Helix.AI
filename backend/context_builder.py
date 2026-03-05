"""
Utility functions to assemble session-aware context that can be injected into
prompts before calling the LLM. Keeps the formatting consistent and easy to
extend.

Token-efficient design: Only inject thin Session Brief (≤800 tokens) with
pointer IDs. Full data is stored in session DB and retrieved on-demand.
"""
from typing import Dict, Any, List, Optional
import json
import hashlib


def _truncate_sequence(seq: str, max_length: int = 100) -> str:
    """Truncate a sequence if it's too long, showing preview and length.
    
    More aggressive truncation (100 chars instead of 200) to reduce LLM payload size.
    """
    if len(seq) <= max_length:
        return seq
    return seq[:max_length] + f"... (truncated, full length: {len(seq)} bp)"


def _format_fasta_from_list(values: List[str], prefix: str) -> str:
    lines = []
    for idx, seq in enumerate(values):
        lines.append(f">{prefix}_{idx+1}")
        # Truncate sequences to avoid sending very long sequences to LLM
        lines.append(_truncate_sequence(seq))
    return "\n".join(lines)


def build_context_snippet(session_context: Dict[str, Any]) -> str:
    """
    Build a short, human-readable context block from session state to guide
    downstream tool calls.
    """
    if not session_context:
        return ""

    parts = []

    mutated = session_context.get("mutated_sequences") or []
    if mutated:
        parts.append("Mutated sequences:\n" + _format_fasta_from_list(mutated, "mutant"))

    aligned = session_context.get("aligned_sequences")
    if aligned:
        # Truncate aligned sequences if they're very long
        # Parse FASTA and truncate each sequence
        if isinstance(aligned, str):
            lines = aligned.split('\n')
            truncated_lines = []
            current_seq = ""
            for line in lines:
                if line.startswith('>'):
                    # Save previous sequence if exists
                    if current_seq:
                        truncated_lines.append(_truncate_sequence(current_seq))
                        current_seq = ""
                    truncated_lines.append(line)
                else:
                    current_seq += line
            # Add last sequence
            if current_seq:
                truncated_lines.append(_truncate_sequence(current_seq))
            truncated_aligned = "\n".join(truncated_lines)
        else:
            truncated_aligned = str(aligned)
        parts.append("Aligned sequences (FASTA, truncated):\n" + truncated_aligned)

    selected = session_context.get("selected_sequences") or []
    if selected:
        parts.append("Selected sequences:\n" + _format_fasta_from_list(selected, "selected"))

    uploaded = session_context.get("uploaded_files") or []
    if uploaded:
        file_names = ", ".join([f.get("name", f.get("filename", "unknown")) for f in uploaded])
        parts.append(f"Uploaded files in session: {file_names}")
    
    # Include dataset references (large files in S3)
    metadata = session_context.get("metadata", {})
    dataset_refs = metadata.get("dataset_references", [])
    if dataset_refs:
        dataset_files = ", ".join([ref.get("filename", "unknown") for ref in dataset_refs])
        dataset_ids = set([ref.get("dataset_id", "") for ref in dataset_refs if ref.get("dataset_id")])
        dataset_info = f"Dataset files: {dataset_files}"
        if dataset_ids:
            dataset_info += f" (from datasets: {', '.join(dataset_ids)})"
        parts.append(dataset_info)

    if not parts:
        return ""

    return "\n\nContext:\n" + "\n\n".join(parts)


def _compute_hash(data: Any) -> str:
    """Compute SHA256 hash of data for reference."""
    if isinstance(data, str):
        return hashlib.sha256(data.encode()).hexdigest()[:16]
    try:
        return hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()[:16]
    except:
        return "unknown"


def _truncate_text(text: str, max_chars: int = 100) -> str:
    """Truncate text to max_chars, showing length."""
    if len(text) <= max_chars:
        return text
    return text[:max_chars] + f"... ({len(text)} chars)"


def build_session_brief(session_context: Dict[str, Any], max_tokens: int = 800) -> str:
    """
    Build a thin Session Brief (≤800 tokens) with pointer IDs only.
    
    Returns compact JSON with:
    - session_id
    - active_goal (1-2 lines)
    - decisions (organism/build/paired/strandedness)
    - inputs (IDs + filenames + sha256 + size)
    - latest_artifacts (IDs + type + title + derived_from + params_hash)
    - latest_runs (IDs + tool + version + params_hash + outputs)
    - open_questions (if any)
    
    Full data is NOT included - only references.
    """
    if not session_context:
        brief = {
            "session_id": None,
            "active_goal": None,
            "decisions": {},
            "inputs": [],
            "latest_artifacts": [],
            "latest_runs": [],
            "open_questions": []
        }
        return json.dumps(brief, indent=2)
    
    brief = {}
    
    # Session ID
    brief["session_id"] = session_context.get("session_id", "unknown")
    
    # Active goal (from most recent command or metadata)
    history = session_context.get("history", [])
    if history:
        latest = history[-1]
        command = latest.get("command", "")
        brief["active_goal"] = _truncate_text(command, 200)
    else:
        brief["active_goal"] = None
    
    # Decisions (extract from metadata or infer from history)
    decisions = {}
    metadata = session_context.get("metadata", {})
    
    # Try to extract common decisions from history
    for entry in reversed(history[-10:]):  # Check last 10 entries
        tool = entry.get("tool", "")
        result = entry.get("result", {})
        
        # Extract organism/build if mentioned
        if "organism" in str(result).lower() or "build" in str(result).lower():
            # Try to extract from result metadata
            if isinstance(result, dict):
                if "organism" in result:
                    decisions["organism"] = result["organism"]
                if "reference_build" in result:
                    decisions["reference_build"] = result["reference_build"]
                if "genome_build" in result:
                    decisions["genome_build"] = result["genome_build"]
        
        # Extract paired/single-end info
        if "paired" in tool.lower() or "single" in tool.lower():
            if "paired" in str(result).lower():
                decisions["library_layout"] = "paired-end"
            elif "single" in str(result).lower():
                decisions["library_layout"] = "single-end"
        
        # Extract strandedness
        if "stranded" in str(result).lower():
            decisions["strandedness"] = "detected"
    
    brief["decisions"] = decisions if decisions else {}
    
    # Inputs (IDs + filenames + hashes + sizes)
    inputs = []
    
    # Uploaded files
    uploaded = session_context.get("uploaded_files", [])
    for i, file_info in enumerate(uploaded[:10]):  # Max 10 files
        input_id = f"uploaded_{i+1}"
        filename = file_info.get("name") or file_info.get("filename", "unknown")
        content = file_info.get("content", "")
        size = len(content) if isinstance(content, str) else file_info.get("size", 0)
        content_hash = _compute_hash(content) if content else "no_content"
        
        inputs.append({
            "input_id": input_id,
            "filename": filename,
            "sha256": content_hash,
            "size_bytes": size,
            "type": "uploaded"
        })
    
    # Dataset references
    dataset_refs = metadata.get("dataset_references", [])
    for i, ref in enumerate(dataset_refs[:10]):  # Max 10 refs
        input_id = f"dataset_{i+1}"
        filename = ref.get("filename", "unknown")
        s3_key = ref.get("s3_key", "")
        size = ref.get("size", 0)
        s3_hash = _compute_hash(s3_key)
        
        inputs.append({
            "input_id": input_id,
            "filename": filename,
            "sha256": s3_hash,
            "size_bytes": size,
            "type": "dataset_reference",
            "s3_key": s3_key[:100]  # Truncate long keys
        })
    
    brief["inputs"] = inputs
    
    # Latest artifacts (IDs + type + title + derived_from + params_hash)
    artifacts = []
    results = session_context.get("results", {})
    
    # Get latest 5 results as artifacts
    result_items = list(results.items())[-5:]
    for result_key, result_data in result_items:
        artifact_type = "unknown"
        title = result_key
        
        if isinstance(result_data, dict):
            artifact_type = result_data.get("type", "unknown")
            title = result_data.get("title", result_key)
        
        # Check if derived from previous step
        derived_from = None
        if history:
            # Try to find parent in previous entries
            for entry in reversed(history[-5:]):
                if entry.get("tool") in result_key:
                    derived_from = entry.get("tool")
                    break
        
        params_hash = _compute_hash(result_data)
        
        artifacts.append({
            "artifact_id": result_key,
            "type": artifact_type,
            "title": title[:100],  # Truncate
            "derived_from": derived_from,
            "params_hash": params_hash
        })
    
    brief["latest_artifacts"] = artifacts
    
    # Latest runs (IDs + tool + version + params_hash + output_refs)
    runs = []
    for i, entry in enumerate(history[-5:]):  # Last 5 runs
        run_id = f"{entry.get('tool', 'unknown')}_{i+1}"
        tool = entry.get("tool", "unknown")
        timestamp = entry.get("timestamp", "")
        result = entry.get("result", {})
        
        # Try to extract version from result
        version = "unknown"
        if isinstance(result, dict):
            version = result.get("version", result.get("tool_version", "unknown"))
        
        params_hash = _compute_hash(entry)
        output_refs = []
        
        # Reference outputs by ID
        if isinstance(result, dict) and "artifact_id" in result:
            output_refs.append(result["artifact_id"])
        
        runs.append({
            "run_id": run_id,
            "tool": tool,
            "version": str(version)[:50],
            "params_hash": params_hash,
            "outputs": output_refs,
            "timestamp": timestamp[:50]
        })
    
    brief["latest_runs"] = runs
    
    # Open questions (from metadata or infer from context)
    brief["open_questions"] = metadata.get("open_questions", [])
    
    # Convert to JSON and check token count (rough estimate: 1 token ≈ 4 chars)
    brief_json = json.dumps(brief, indent=2)
    estimated_tokens = len(brief_json) // 4
    
    # If over limit, truncate aggressively
    if estimated_tokens > max_tokens:
        # Keep only most recent items
        brief["latest_artifacts"] = brief["latest_artifacts"][-3:]
        brief["latest_runs"] = brief["latest_runs"][-3:]
        brief["inputs"] = brief["inputs"][-5:]
        brief_json = json.dumps(brief, indent=2)
    
    return brief_json


def build_session_context_summary(session_context: Dict[str, Any]) -> str:
    """
    DEPRECATED: Use build_session_brief() instead for token-efficient injection.
    Kept for backward compatibility.
    """
    return build_session_brief(session_context)
