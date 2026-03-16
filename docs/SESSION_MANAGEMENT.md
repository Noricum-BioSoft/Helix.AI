# Session Management and Context Awareness

## Overview

Helix.AI operates within a **session-aware environment** where each user interaction occurs within a persistent session. This document defines how session context is structured, maintained, and used across the system to enable multi-turn conversations and workflow continuity.

---

## Session Context Availability

Each user session maintains:

- **Previous results**: Outputs from prior tool executions in the session
- **History**: Chronological record of all commands and their results
- **Uploaded files**: Files provided by the user in this session
- **Dataset references**: Links to large datasets stored in S3
- **Intermediate artifacts**: Sequences, alignments, tables, and other data generated during workflows

---

## Session Context Structure

The session context dictionary contains:

```python
session_context = {
    # Core identifiers
    "session_id": "sess_abc123",  # Unique session identifier
    
    # Workflow artifacts
    "mutated_sequences": ["ATGC...", "ATGA..."],  # From mutation operations
    "aligned_sequences": "ATGC-\nATGA-\n...",     # FASTA-formatted alignment
    "selected_sequences": ["seq1", "seq2"],        # Selected/variant sequences
    "phylogenetic_tree": "(A:0.1,B:0.2);",        # Newick format tree
    
    # File management
    "uploaded_files": [
        {
            "name": "sample.fastq",
            "size": 500000000,
            "s3_bucket": "helix-uploads",
            "s3_key": "uploads/sess_abc123/sample.fastq",
            "upload_time": "2026-01-18T10:00:00Z"
        }
    ],
    
    # Execution history
    "history": [
        {
            "timestamp": "2026-01-18T10:00:00Z",
            "command": "Align these sequences",
            "tool": "mafft",
            "status": "success",
            "result_id": "alignment_1"
        }
    ],
    
    # Results registry
    "results": {
        "alignment_1": {
            "tool": "mafft",
            "inputs": ["seq1.fasta"],
            "outputs": ["aligned.fasta"],
            "artifacts": [...]
        },
        "tree_2": {
            "tool": "iqtree",
            "inputs": ["aligned.fasta"],
            "outputs": ["tree.nwk"],
            "artifacts": [...]
        }
    },
    
    # Metadata
    "metadata": {
        "s3_bucket": "helix-data",
        "s3_path": "sessions/sess_abc123/",
        "dataset_references": [
            {
                "name": "Reference genome hg38",
                "uri": "s3://public-datasets/genomes/hg38.fa",
                "size": 3000000000
            }
        ]
    }
}
```

---

## Using Session Context

### Rule 1: Always Check Session Context First

Before requesting data from the user, check if it already exists in session context:

```python
# ✅ Good: Check session context first
if "aligned_sequences" in session_context:
    alignment = session_context["aligned_sequences"]
else:
    # Only ask user if not in session
    request_input("Please provide aligned sequences")

# ❌ Bad: Always ask user
request_input("Please provide aligned sequences")
```

### Rule 2: Reference Previous Results

When the user refers to previous operations, look up results in session context:

```python
# User says: "Build a tree from the alignment"

# ✅ Good: Look up alignment from session
if "aligned_sequences" in session_context:
    alignment = session_context["aligned_sequences"]
    build_tree(alignment)
else:
    # Check results registry
    if "alignment_1" in session_context.get("results", {}):
        alignment = session_context["results"]["alignment_1"]["outputs"][0]
        build_tree(alignment)
    else:
        clarify("I don't see an alignment in the session. Please provide sequences.")
```

### Rule 3: Continue Workflows

Chain operations using session data:

```python
# User's first command: "Mutate these sequences"
result1 = mutate_sequences(sequences)
session_context["mutated_sequences"] = result1

# User's second command: "Now align the mutants"
# ✅ Good: Use mutated sequences from session
sequences_to_align = session_context["mutated_sequences"]
result2 = align_sequences(sequences_to_align)
session_context["aligned_sequences"] = result2

# User's third command: "Build a phylogenetic tree"
# ✅ Good: Use aligned sequences from session
alignment = session_context["aligned_sequences"]
result3 = build_tree(alignment)
```

### Rule 4: Leverage Uploaded Files

Check session metadata for uploaded files before requesting new uploads:

```python
# ✅ Good: Check for uploaded files
uploaded_fastq = [
    f for f in session_context.get("metadata", {}).get("uploaded_files", [])
    if f["name"].endswith((".fastq", ".fq"))
]

if uploaded_fastq:
    # Use uploaded files
    analyze_fastq(uploaded_fastq[0])
else:
    # Request upload
    request_upload("Please upload FASTQ files")
```

### Rule 5: Maintain Continuity

Reference previous steps in explanations and provenance:

```python
# ✅ Good: Reference previous operations
message = "Using the alignment from step 2 (MAFFT, completed at 10:05 AM), I built a phylogenetic tree using IQ-TREE."

# ❌ Bad: No reference to previous context
message = "I built a phylogenetic tree."
```

---

## Session Brief and Retrieval Rules

### Session Brief Injection

A compact **Session Brief** (≤800 tokens) is prepended to each user message:

```python
session_brief = {
    "session_id": "sess_abc123",
    "active_goal": "RNA-seq differential expression analysis",
    "decisions": ["organism=human", "reference=hg38"],
    "input_ids": ["fastq_001", "fastq_002"],
    "latest_artifact_ids": ["qc_report_1", "alignment_2"],
    "latest_run_ids": ["mafft_0102", "star_0103"],
    "open_questions": ["Which sample groups for DE?"]
}
```

**Session Brief contains only pointers (IDs, filenames, hashes, sizes) — not full data.**

### Retrieval Rules

1. **If a needed detail is not in the Session Brief, retrieve it via session tools rather than asking the user or guessing.**

   ```python
   # ✅ Good: Retrieve from session
   if artifact_id in session_brief["latest_artifact_ids"]:
       full_data = session.get_artifact(artifact_id)
   
   # ❌ Bad: Ask user for data that exists in session
   request_input("Please provide the QC report")
   ```

2. **If retrieval is expensive, ask a minimal question instead of fetching everything.**

   ```python
   # ✅ Good: Ask specific question
   ask_user("Which samples should I use for differential expression?")
   
   # ❌ Bad: Fetch all data unnecessarily
   all_samples = session.get_all_samples()  # Expensive if 1000s of samples
   ```

3. **Never include large blobs in the final JSON — attach as artifacts by URI or store references.**

   ```python
   # ✅ Good: Reference by URI
   output = {
       "artifact": {
           "type": "alignment",
           "uri": "s3://helix-data/sessions/sess_abc123/alignment.fasta",
           "size_bytes": 500000000
       }
   }
   
   # ❌ Bad: Embed large data
   output = {
       "artifact": {
           "type": "alignment",
           "content": "ATGC..." * 1000000  # 500MB inline
       }
   }
   ```

4. **Use pointer IDs everywhere**: Instead of including data, use references.

   ```python
   # ✅ Good: Use IDs
   {
       "artifact_id": "aln_0239",
       "run_id": "mafft_0102",
       "input_id": "fastq_001"
   }
   
   # ❌ Bad: Embed data
   {
       "artifact": "ATGC...",
       "run_log": "Full log...",
       "input_content": "Full FASTQ..."
   }
   ```

5. **When you need full data**, use session retrieval tools.

   ```python
   # Retrieve specific artifacts
   alignment = session.get_artifact("aln_0239")
   run_log = session.get_run("mafft_0102")
   fastq_metadata = session.get_input("fastq_001")
   ```

---

## Token Limits (Hard Caps)

To keep responses performant and cost-effective:

- **Session Brief**: max 800 tokens
- **Sequences inline**: max 10–50 KB; else use URI
- **Tables inline**: max 500–2,000 rows; else use URI
- **Charts**: downsample points (max 20k points)

### What NOT to Include in Session Brief

❌ Never include:
- Raw FASTA/FASTQ/VCF/BAM content
- Full DE tables / marker gene lists
- Long tool logs
- Full prior chat transcripts

✅ Instead, store in S3 and reference by URI:
- `s3://helix-data/sessions/sess_abc123/alignment.fasta`
- `s3://helix-data/sessions/sess_abc123/de_results.csv`
- `s3://helix-data/sessions/sess_abc123/logs/star_alignment.log`

---

## Context-Aware Workflow Examples

### Example 1: Sequential Operations

```
User: "Mutate these sequences"
→ System: Mutates sequences → Stores in session_context["mutated_sequences"]

User: "Align the mutants"
→ System: Retrieves session_context["mutated_sequences"] → Aligns → Stores as session_context["aligned_sequences"]

User: "Build a phylogenetic tree"
→ System: Retrieves session_context["aligned_sequences"] → Builds tree → Stores as session_context["phylogenetic_tree"]
```

### Example 2: Iterative Refinement

```
User: "Select the top 10 variants"
→ System: Selects top 10 → Stores in session_context["selected_sequences"]

User: "Now show me the top 5"
→ System: Retrieves session_context["selected_sequences"] → Filters to top 5 → Returns subset
```

### Example 3: File Reuse

```
User: "Analyze this FASTQ file" [uploads sample.fastq]
→ System: Uploads to S3 → Stores metadata in session_context["uploaded_files"]

User: "Now trim the same file"
→ System: Retrieves file from session_context["uploaded_files"] → Trims → Returns results
```

### Example 4: Cross-Tool Workflows

```
User: "Align these sequences"
→ System: Aligns → Stores alignment in session_context["aligned_sequences"]

User: "Select representative sequences"
→ System: Retrieves alignment → Selects representatives → Stores in session_context["selected_sequences"]

User: "Visualize as plasmid map"
→ System: Retrieves selected_sequences → Creates plasmid map → Returns visualization
```

---

## When to Request New Inputs

Only request new inputs when:

1. **The session context doesn't contain the required data** (check Session Brief first)
2. **The user explicitly provides new data**
3. **The user wants to start a new analysis branch**

### Examples

✅ **Should request input:**
```
User: "Start a new RNA-seq analysis"
→ System: Starts fresh, requests new FASTQ files
```

❌ **Should NOT request input:**
```
User: "Visualize the previous results"
→ System should check session context for results, not ask user to provide them again
```

---

## Session Management Best Practices

### For Agent Developers

1. **Always check session context before tool execution**
   - Check for existing artifacts that can be reused
   - Avoid redundant operations

2. **Store intermediate results**
   - After each significant operation, store results in session context
   - Use descriptive keys (e.g., `aligned_sequences`, not `temp_data`)

3. **Use pointer-based references**
   - Store large data in S3, keep only URIs/IDs in session context
   - Keep session context lean (<800 tokens for brief)

4. **Maintain provenance**
   - Record which operations produced which artifacts
   - Link downstream operations to upstream results

5. **Clean up expired data**
   - Remove old artifacts when no longer needed
   - Implement session expiration policies

### For System Architects

1. **Session storage backend**
   - Use Redis/DynamoDB for session metadata (fast lookup)
   - Use S3 for large artifacts (cost-effective storage)

2. **Session expiration**
   - Set reasonable TTL (e.g., 24 hours for active sessions)
   - Archive old sessions to cold storage

3. **Concurrency handling**
   - Use locks for concurrent session updates
   - Implement optimistic concurrency control

4. **Garbage collection**
   - Clean up orphaned artifacts
   - Remove unused S3 objects after session expiration

---

## Related Documents

- `docs/TASK_ARCHITECTURE.md` - Micro-/Macroflow pattern
- `docs/OUTPUT_SCHEMA.md` - Output format and artifact schemas
- `backend/history_manager.py` - Session history implementation
- `agents/workflow-planner-agent.md` - Workflow planning with session awareness

---

## Glossary

- **Session Context**: Persistent state maintained across user interactions
- **Session Brief**: Compact summary (≤800 tokens) of session state
- **Artifact**: Output from a tool execution (alignment, table, plot, etc.)
- **Pointer**: Reference to data (ID, URI, hash) instead of embedding data
- **Provenance**: Complete history of operations and their relationships
