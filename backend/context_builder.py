"""
Utility functions to assemble session-aware context that can be injected into
prompts before calling the LLM. Keeps the formatting consistent and easy to
extend.
"""
from typing import Dict, Any, List


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
