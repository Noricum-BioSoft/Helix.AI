"""
Central registry for tool metadata (inputs/outputs) to keep MCP listings and
agent wiring in sync.
"""
from typing import Dict, List, Any


TOOL_SCHEMAS: List[Dict[str, Any]] = [
    {
        "name": "toolbox_inventory",
        "description": "List all tools Helix.AI has access to (registered MCP tools, discovered @tool functions, and local/EC2 CLI tools).",
        "inputs": {
            "type": "object",
            "properties": {},
        },
    },
    {
        "name": "sequence_alignment",
        "description": "Perform multiple sequence alignment on DNA/RNA sequences.",
        "inputs": {
            "type": "object",
            "properties": {
                "sequences": {"type": "string", "description": "FASTA formatted sequences"},
                "algorithm": {"type": "string", "enum": ["clustal", "muscle", "mafft"], "default": "clustal"},
            },
            "required": ["sequences"],
        },
    },
    {
        "name": "mutate_sequence",
        "description": "Generate sequence variants.",
        "inputs": {
            "type": "object",
            "properties": {
                "sequence": {"type": "string"},
                "num_variants": {"type": "integer", "default": 96},
                "mutation_rate": {"type": "number", "default": 0.1},
            },
            "required": ["sequence"],
        },
    },
    {
        "name": "variant_selection",
        "description": "Select variants from session mutations.",
        "inputs": {
            "type": "object",
            "properties": {
                "session_id": {"type": "string"},
                "selection_criteria": {"type": "string", "default": "diversity"},
                "num_variants": {"type": "integer", "default": 10},
                "custom_filters": {"type": "object"},
            },
            "required": ["session_id"],
        },
    },
    {
        "name": "phylogenetic_tree",
        "description": "Create a phylogenetic tree from aligned sequences.",
        "inputs": {
            "type": "object",
            "properties": {
                "aligned_sequences": {"type": "string", "description": "FASTA formatted sequences"},
            },
            "required": ["aligned_sequences"],
        },
    },
    {
        "name": "fetch_ncbi_sequence",
        "description": "Fetch sequence by accession from NCBI.",
        "inputs": {
            "type": "object",
            "properties": {
                "accession": {"type": "string"},
                "database": {"type": "string", "enum": ["nucleotide", "protein"], "default": "nucleotide"},
            },
            "required": ["accession"],
        },
    },
    {
        "name": "query_uniprot",
        "description": "Query UniProt for protein sequences and metadata.",
        "inputs": {
            "type": "object",
            "properties": {
                "query": {"type": "string"},
                "format": {"type": "string", "default": "fasta"},
                "limit": {"type": "integer", "default": 10},
            },
            "required": ["query"],
        },
    },
    {
        "name": "lookup_go_term",
        "description": "Lookup Gene Ontology term details.",
        "inputs": {
            "type": "object",
            "properties": {
                "go_id": {"type": "string"},
            },
            "required": ["go_id"],
        },
    },
    {
        "name": "bulk_rnaseq_analysis",
        "description": "Run DESeq2 bulk RNA-seq differential expression.",
        "inputs": {
            "type": "object",
            "properties": {
                "count_matrix": {"type": "string"},
                "sample_metadata": {"type": "string"},
                "design_formula": {"type": "string", "default": "~condition"},
                "alpha": {"type": "number", "default": 0.05},
            },
            "required": ["count_matrix", "sample_metadata"],
        },
    },
    {
        "name": "plasmid_visualization",
        "description": "Generate plasmid visualization data.",
        "inputs": {
            "type": "object",
            "properties": {
                "vector_name": {"type": "string"},
                "cloning_sites": {"type": "string"},
                "insert_sequence": {"type": "string"},
                "full_plasmid_sequence": {"type": "string"},
                "insert_position": {"type": "integer"},
            },
        },
    },
    {
        "name": "single_cell_analysis",
        "description": "Perform single-cell RNA-seq analysis or answer SC questions.",
        "inputs": {
            "type": "object",
            "properties": {
                "data_file": {"type": "string"},
                "data_format": {"type": "string", "enum": ["10x", "h5", "csv", "seurat"], "default": "10x"},
                "steps": {"type": "array", "items": {"type": "string"}},
                "resolution": {"type": "number", "default": 0.5},
                "nfeatures": {"type": "integer", "default": 2000},
            },
        },
    },
    {
        "name": "dna_vendor_research",
        "description": "Research DNA synthesis vendors and testing options.",
        "inputs": {
            "type": "object",
            "properties": {
                "command": {"type": "string"},
                "sequence_length": {"type": "integer"},
                "quantity": {"type": "string", "default": "standard"},
            },
            "required": ["command"],
        },
    },
    {
        "name": "read_trimming",
        "description": "Trim adapters/low-quality bases from sequencing reads.",
        "inputs": {
            "type": "object",
            "properties": {
                "reads": {"type": "string"},
                "forward_reads": {"type": "string"},
                "reverse_reads": {"type": "string"},
                "adapter": {"type": "string"},
                "quality_threshold": {"type": "integer", "default": 20},
            },
        },
    },
    {
        "name": "read_merging",
        "description": "Merge paired-end reads using overlap consensus.",
        "inputs": {
            "type": "object",
            "properties": {
                "forward_reads": {"type": "string"},
                "reverse_reads": {"type": "string"},
                "min_overlap": {"type": "integer", "default": 12},
            },
            "required": ["forward_reads", "reverse_reads"],
        },
    },
    {
        "name": "sequence_selection",
        "description": "Select sequences from aligned sequences based on various criteria.",
        "inputs": {
            "type": "object",
            "properties": {
                "aligned_sequences": {"type": "string", "description": "FASTA formatted aligned sequences"},
                "selection_type": {"type": "string", "enum": ["random", "best_conservation", "lowest_gaps", "highest_gc", "longest", "shortest"], "default": "random"},
                "num_sequences": {"type": "integer", "default": 1},
            },
            "required": ["aligned_sequences"],
        },
    },
    {
        "name": "synthesis_submission",
        "description": "Submit sequences for DNA synthesis and get pricing quote.",
        "inputs": {
            "type": "object",
            "properties": {
                "sequences": {"type": "string", "description": "FASTA formatted sequences"},
                "vendor_preference": {"type": "string"},
                "quantity": {"type": "string", "enum": ["standard", "large", "custom"], "default": "standard"},
                "delivery_time": {"type": "string", "enum": ["rush", "standard", "economy"], "default": "standard"},
            },
            "required": ["sequences"],
        },
    },
    {
        "name": "plasmid_for_representatives",
        "description": "Create plasmid visualizations for representative sequences from clustering analysis.",
        "inputs": {
            "type": "object",
            "properties": {
                "representatives": {"type": "array", "items": {"type": "string"}},
                "aligned_sequences": {"type": "string", "description": "FASTA formatted sequences"},
                "vector_name": {"type": "string", "default": "pUC19"},
                "cloning_sites": {"type": "string", "default": "EcoRI, BamHI, HindIII"},
            },
            "required": ["representatives", "aligned_sequences"],
        },
    },
    {
        "name": "fastqc_quality_analysis",
        "description": "Run FastQC quality control analysis on paired-end FASTQ files stored in S3 using AWS EMR. This is an asynchronous long-running job that typically takes 10-30 minutes. Returns immediately with a job_id for tracking progress.",
        "inputs": {
            "type": "object",
            "properties": {
                "input_r1": {
                    "type": "string",
                    "description": "S3 URI for the forward/read 1 FASTQ file (e.g., s3://my-bucket/data/sample_R1.fastq)"
                },
                "input_r2": {
                    "type": "string",
                    "description": "S3 URI for the reverse/read 2 FASTQ file (e.g., s3://my-bucket/data/sample_R2.fastq)"
                },
                "output": {
                    "type": "string",
                    "description": "Optional S3 URI for output directory (defaults to same directory as input files with /fastqc-results suffix)"
                },
            },
            "required": ["input_r1", "input_r2"],
        },
    },
]


def list_tool_schemas() -> List[Dict[str, Any]]:
    return TOOL_SCHEMAS


def get_tool_schema(name: str) -> Dict[str, Any]:
    for schema in TOOL_SCHEMAS:
        if schema["name"] == name:
            return schema
    return {}
