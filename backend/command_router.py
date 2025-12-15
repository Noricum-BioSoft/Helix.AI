import re
from typing import Dict, Any, Optional, Tuple
import os

# Import directed evolution handler
from backend.directed_evolution_handler import DirectedEvolutionHandler

class CommandRouter:
    """Simple keyword-based command router to replace ChatGPT."""
    
    def __init__(self):
        self.tool_mappings = {
            'sequence_alignment': {
                'keywords': ['align', 'alignment', 'compare sequences', 'multiple sequence alignment'],
                'description': 'Align DNA/RNA sequences'
            },
            'mutate_sequence': {
                'keywords': ['mutate', 'mutation', 'variant', 'create variants', 'generate variants'],
                'description': 'Create sequence variants'
            },
            'sequence_selection': {
                'keywords': ['pick', 'select', 'choose', 'from the', 'randomly pick', 'select from', 'choose from'],
                'description': 'Select sequences from alignment'
            },
            'phylogenetic_tree': {
                'keywords': ['phylogenetic', 'tree', 'evolutionary tree', 'phylogeny'],
                'description': 'Create phylogenetic tree'
            },
            'dna_vendor_research': {
                'keywords': ['order', 'vendor', 'synthesis', 'test', 'assay', 'function', 'binding'],
                'description': 'Research DNA vendors and testing'
            },
            'synthesis_submission': {
                'keywords': ['submit', 'synthesis', 'order sequences'],
                'description': 'Submit sequences for synthesis'
            },
            'plasmid_visualization': {
                'keywords': ['plasmid', 'vector', 'cloning'],
                'description': 'Visualize plasmid constructs'
            },
            'plasmid_for_representatives': {
                'keywords': ['insert representatives', 'express representatives', 'clone representatives', 'vector representatives', 'plasmid representatives'],
                'description': 'Create plasmid visualizations for representative sequences from clustering'
            },
            'directed_evolution': {
                'keywords': ['directed evolution', 'evolution', 'protein engineering', 'dbtl', 'design build test learn'],
                'description': 'Directed evolution for protein engineering'
            },
            'single_cell_analysis': {
                'keywords': ['single cell', 'scRNA-seq', 'scRNAseq', 'single-cell', 'cell type', 'marker genes', 'differential expression', 'pathway analysis', 'batch correction', 'seurat', 'scpipeline'],
                'description': 'Single-cell RNA-seq analysis using scPipeline'
            },
            'fetch_ncbi_sequence': {
                'keywords': ['fetch sequence', 'get sequence', 'ncbi', 'accession', 'genbank', 'refseq', 'download sequence'],
                'description': 'Fetch sequences from NCBI by accession'
            },
            'query_uniprot': {
                'keywords': ['uniprot', 'protein database', 'query protein', 'get protein'],
                'description': 'Query UniProt for protein sequences'
            },
            'lookup_go_term': {
                'keywords': ['go term', 'gene ontology', 'lookup go', 'go:'],
                'description': 'Lookup Gene Ontology terms'
            },
            'bulk_rnaseq_analysis': {
                'keywords': ['bulk rna-seq', 'deseq2', 'differential expression', 'rna-seq analysis', 'transcriptomics'],
                'description': 'Run bulk RNA-seq differential expression analysis'
            }
        }
        
        # Initialize directed evolution handler
        self.de_handler = DirectedEvolutionHandler()
    
    def route_command(self, command: str, session_context: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
        """
        Route a command to the appropriate tool based on keywords.
        Returns (tool_name, parameters)
        """
        command_lower = command.lower()
        
        # Tool inventory / "what tools do you have?" (HIGHEST PRIORITY)
        if any(
            phrase in command_lower
            for phrase in [
                "what tools do you have",
                "what tools are available",
                "list tools",
                "show tools",
                "your toolbox",
                "toolbox",
                "capabilities",
                "what can you do",
            ]
        ):
            print("🔧 Command router: Matched 'toolbox inventory' -> toolbox_inventory")
            return "toolbox_inventory", {}

        # Priority-based matching to avoid conflicts
        # Check for specific phrases first, then general keywords
        
        # Check for clustering analysis (HIGHEST PRIORITY - before variant selection)
        if any(phrase in command_lower for phrase in ['cluster', 'clustering', 'group sequences', 'cluster sequences', 'representative sequences']):
            print(f"🔧 Command router: Matched 'clustering' -> clustering_analysis")
            return 'clustering_analysis', self._extract_parameters(command, 'clustering_analysis', session_context)
        
        # Check for variant selection (HIGH PRIORITY - before phylogenetic tree)
        # Check for "select X sequences" patterns (but not representative sequences from clustering)
        if re.search(r'select\s+\d+\s+sequences?', command_lower) and 'representative sequences' not in command_lower:
            print(f"🔧 Command router: Matched 'select sequences' -> variant_selection")
            return 'variant_selection', self._extract_parameters(command, 'variant_selection', session_context)
        if any(phrase in command_lower for phrase in ['select variants', 'top variants', 'representative variants', 'diverse variants', 'select top']):
            print(f"🔧 Command router: Matched 'variant selection' -> variant_selection")
            return 'variant_selection', self._extract_parameters(command, 'variant_selection', session_context)
        
        # Check for phylogenetic tree first (highest priority)
        if any(phrase in command_lower for phrase in ['phylogenetic tree', 'evolutionary tree', 'phylogeny']):
            print(f"🔧 Command router: Matched 'phylogenetic tree' -> phylogenetic_tree")
            return 'phylogenetic_tree', self._extract_parameters(command, 'phylogenetic_tree', session_context)
        
        # Check for mutation/variant generation
        if any(phrase in command_lower for phrase in ['mutate', 'mutation', 'variant', 'create variants', 'generate variants']):
            print(f"🔧 Command router: Matched 'mutation' -> mutate_sequence")
            return 'mutate_sequence', self._extract_parameters(command, 'mutate_sequence', session_context)
        
        # Check for sequence alignment (more specific to avoid false matches)
        if any(phrase in command_lower for phrase in ['align sequences', 'sequence alignment', 'compare sequences', 'multiple sequence alignment', 'perform alignment', 'perform multiple sequence alignment']):
            print(f"🔧 Command router: Matched 'alignment' -> sequence_alignment")
            return 'sequence_alignment', self._extract_parameters(command, 'sequence_alignment', session_context)
        
        # Check for bulk RNA-seq analysis (HIGH PRIORITY to avoid misrouting to vendor due to 'expression')
        if any(phrase in command_lower for phrase in ['bulk rna-seq', 'deseq2', 'differential expression', 'rna-seq analysis', 'transcriptomics']):
            print(f"🔧 Command router: Matched 'bulk RNA-seq' -> bulk_rnaseq_analysis")
            return 'bulk_rnaseq_analysis', self._extract_parameters(command, 'bulk_rnaseq_analysis', session_context)
        
        # Check for FastQC analysis (HIGH PRIORITY - before vendor check to avoid false matches with "test" in paths)
        if any(phrase in command_lower for phrase in ['fastqc', 'fastqc analysis', 'quality control analysis', 'perform fastqc', 'run fastqc']):
            print(f"🔧 Command router: Matched 'FastQC' -> fastqc_quality_analysis")
            return 'fastqc_quality_analysis', self._extract_parameters(command, 'fastqc_quality_analysis', session_context)
        
        # Check for vendor research (exclude "test" if it's in a file path)
        # Check if "test" appears in an S3 path or file path context
        is_test_in_path = bool(re.search(r's3://[^/]+/[^/]*test[^/]*/', command_lower) or 
                               re.search(r'[/\\][^/\\]*test[^/\\]*[/\\]', command_lower))
        
        vendor_keywords = ['order', 'vendor', 'synthesis', 'assay', 'function', 'binding']
        # Only include "test" if it's not in a file path
        if not is_test_in_path:
            vendor_keywords.append('test')
        
        if any(phrase in command_lower for phrase in vendor_keywords):
            print(f"🔧 Command router: Matched 'vendor' -> dna_vendor_research")
            return 'dna_vendor_research', self._extract_parameters(command, 'dna_vendor_research', session_context)

        # Check for plasmid for representatives (HIGH PRIORITY - before general plasmid)
        if any(phrase in command_lower for phrase in ['insert representatives', 'express representatives', 'clone representatives', 'vector representatives', 'plasmid representatives']):
            print(f"🔧 Command router: Matched 'plasmid representatives' -> plasmid_for_representatives")
            return 'plasmid_for_representatives', self._extract_parameters(command, 'plasmid_for_representatives', session_context)
        
        # Check for plasmid visualization (HIGH PRIORITY - before alignment fallback)
        if any(phrase in command_lower for phrase in ['plasmid', 'vector', 'cloning', 'insert', 'express', 'into vector', 'into plasmid']):
            print(f"🔧 Command router: Matched 'plasmid' -> plasmid_visualization")
            return 'plasmid_visualization', self._extract_parameters(command, 'plasmid_visualization', session_context)
        
        # Check for synthesis submission
        if any(phrase in command_lower for phrase in ['submit', 'synthesis', 'order sequences']):
            print(f"🔧 Command router: Matched 'synthesis' -> synthesis_submission")
            return 'synthesis_submission', self._extract_parameters(command, 'synthesis_submission', session_context)
        
        # Check for session creation
        if any(phrase in command_lower for phrase in ['create session', 'new session', 'start session', 'initialize session']):
            print(f"🔧 Command router: Matched 'session creation' -> session_creation")
            return 'session_creation', {"command": command}
        
        # Check for single-cell analysis (HIGH PRIORITY - before other checks)
        if any(phrase in command_lower for phrase in ['single cell', 'scrna-seq', 'scrnaseq', 'single-cell', 'cell type', 'marker genes', 'differential expression', 'pathway analysis', 'batch correction', 'seurat', 'scpipeline']):
            print(f"🔧 Command router: Matched 'single-cell analysis' -> single_cell_analysis")
            return 'single_cell_analysis', self._extract_parameters(command, 'single_cell_analysis', session_context)
        
        # Check for NCBI sequence fetch (HIGH PRIORITY)
        if any(phrase in command_lower for phrase in ['fetch sequence', 'get sequence', 'ncbi', 'accession', 'genbank', 'refseq', 'download sequence']):
            print(f"🔧 Command router: Matched 'NCBI sequence fetch' -> fetch_ncbi_sequence")
            return 'fetch_ncbi_sequence', self._extract_parameters(command, 'fetch_ncbi_sequence', session_context)
        
        # Check for UniProt queries
        if any(phrase in command_lower for phrase in ['uniprot', 'protein database', 'query protein', 'get protein']):
            print(f"🔧 Command router: Matched 'UniProt query' -> query_uniprot")
            return 'query_uniprot', self._extract_parameters(command, 'query_uniprot', session_context)
        
        # Check for GO term lookups
        if any(phrase in command_lower for phrase in ['go term', 'gene ontology', 'go:', 'lookup go']):
            print(f"🔧 Command router: Matched 'GO term lookup' -> lookup_go_term")
            return 'lookup_go_term', self._extract_parameters(command, 'lookup_go_term', session_context)
        
        # Check for directed evolution
        if any(phrase in command_lower for phrase in ['directed evolution', 'evolution', 'protein engineering', 'dbtl', 'design build test learn']):
            print(f"🔧 Command router: Matched 'directed evolution' -> directed_evolution")
            return 'directed_evolution', self._extract_parameters(command, 'directed_evolution', session_context)
        
        # Check for quality assessment (HIGH PRIORITY - before other checks)
        if any(phrase in command_lower for phrase in ['quality report', 'quality assessment', 'generate quality', 'quality metrics', 'qc report', 'quality check', 'assess quality']):
            print(f"🔧 Command router: Matched 'quality assessment' -> quality_assessment")
            return 'quality_assessment', self._extract_parameters(command, 'quality_assessment', session_context)
        
        # Check for read merging (HIGH PRIORITY - before trimming to avoid conflicts)
        if any(phrase in command_lower for phrase in ['merge', 'merging', 'merge reads', 'merge paired', 'merge paired-end', 'merge my paired']):
            print(f"🔧 Command router: Matched 'read merging' -> read_merging")
            return 'read_merging', self._extract_parameters(command, 'read_merging', session_context)
        
        # Check for read trimming (HIGH PRIORITY - before alignment fallback)
        if any(phrase in command_lower for phrase in ['trim', 'trimming', 'trim reads', 'trim low-quality', 'quality trim', 'quality threshold', 'remove adapter', 'adapter removal']):
            print(f"🔧 Command router: Matched 'read trimming' -> read_trimming")
            return 'read_trimming', self._extract_parameters(command, 'read_trimming', session_context)
        
        # Default fallback - try to route to sequence_alignment for alignment-like commands
        if any(phrase in command_lower for phrase in ['align', 'alignment', 'sequences']):
            print(f"🔧 Command router: Defaulting to sequence_alignment for alignment command")
            return "sequence_alignment", self._extract_parameters(command, 'sequence_alignment', session_context)
        
        # For other commands, try to use the natural command handler
        print(f"🔧 Command router: No specific match, using natural command handler")
        return "handle_natural_command", {"command": command, "session_id": session_context.get("session_id", "")}
        
    def route_plan(self, command: str, session_context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Phase 3: Emit a minimal Plan IR for multi-step workflows.

        This is intentionally heuristic-based (non-LLM) so it works in mock/CI.
        We split on common workflow delimiters (then, and then, ->, ;, newlines) and map each chunk
        through the existing router.
        """
        from backend.plan_ir import Plan, PlanStep

        parts = self._split_workflow_command(command)
        steps = []
        for idx, part in enumerate(parts, start=1):
            tool_name, params = self.route_command(part, session_context)
            steps.append(
                PlanStep(
                    id=f"step{idx}",
                    tool_name=tool_name,
                    arguments=params or {},
                    description=part.strip(),
                )
            )
        return Plan(steps=steps).dict()

    def _split_workflow_command(self, command: str) -> list[str]:
        if not command:
            return []
        # Normalize separators
        text = command.replace("->", "\n").replace("→", "\n")
        # Split on "then"/"and then" plus newlines/semicolons
        chunks = re.split(r"(?:\bthen\b|\band then\b|;|\n)+", text, flags=re.IGNORECASE)
        parts = [c.strip() for c in chunks if c and c.strip()]
        # If no real split happened, return the original as a single step
        return parts or [command.strip()]
    
    def _extract_parameters(self, command: str, tool_name: str, session_context: Dict[str, Any]) -> Dict[str, Any]:
        """Extract parameters for the specific tool."""
        command_lower = command.lower()
        
        if tool_name == "mutate_sequence":
            # Extract sequence and number of variants
            # Try multiple patterns to extract the full sequence
            sequence = None
            
            # Pattern 1: "Sequence: ATGCGATCG..." (with colon) - allow spaces in sequence
            sequence_match = re.search(r'sequence\s*:\s*([ATCG\s]+?)(?:\s|$|\.|,|;|variants|variant)', command, re.IGNORECASE)
            if sequence_match:
                sequence = sequence_match.group(1).strip()
            
            # Pattern 2: "sequence ATGCGATCG..." (without colon, but with space) - allow spaces
            if not sequence:
                sequence_match = re.search(r'sequence\s+([ATCG\s]+?)(?:\s|$|\.|,|;|variants|variant)', command, re.IGNORECASE)
                if sequence_match:
                    sequence = sequence_match.group(1).strip()
            
            # Pattern 3: Look for long DNA sequences anywhere in the command (no spaces, for fallback)
            if not sequence:
                # Find the longest DNA sequence in the command
                sequences = re.findall(r'[ATCG]{10,}', command.upper())
                if sequences:
                    # Use the longest sequence found
                    sequence = max(sequences, key=len)
            
            # Fallback to default if nothing found
            if not sequence:
                sequence = "ATGCGATCG"
            
            # Clean the sequence (remove spaces) before returning
            sequence = sequence.replace(" ", "").upper()
            
            variants_match = re.search(r'(\d+)\s+variants?', command)
            num_variants = int(variants_match.group(1)) if variants_match else 96
            
            return {
                "sequence": sequence,
                "num_variants": num_variants
            }
        
        elif tool_name == "sequence_alignment":
            # Use sequences from session context or extract from command
            if session_context.get("mutated_sequences"):
                # Format as proper FASTA with each sequence on its own line after header
                fasta_sequences = []
                for i, seq in enumerate(session_context["mutated_sequences"]):
                    fasta_sequences.append(f">seq_{i+1}")
                    fasta_sequences.append(seq)
                sequences = "\n".join(fasta_sequences)
                print(f"[DEBUG] FASTA string sent to alignment tool: {len(sequences)} chars, {len(session_context['mutated_sequences'])} sequences")
                return {"sequences": sequences}
            else:
                # Extract sequences from command - improved pattern to handle FASTA format
                # Try multiple patterns to extract sequences
                sequences = None
                
                # Pattern 1: Extract FASTA format sequences (handles blank lines)
                # Match from first > to end, including blank lines and DNA sequences
                # Try to find FASTA content anywhere in the command
                fasta_match = re.search(r'(>[\w\s\n>ATCGUatcgu]+)', command, re.DOTALL)
                if fasta_match:
                    sequences = fasta_match.group(1).strip()
                    # Clean up blank lines - remove blank lines between headers and sequences
                    lines = sequences.split('\n')
                    cleaned_lines = []
                    for line in lines:
                        stripped = line.strip()
                        if stripped:  # Only keep non-empty lines
                            cleaned_lines.append(stripped)
                    sequences = '\n'.join(cleaned_lines)
                    print(f"🔧 Extracted FASTA sequences: {len(sequences)} chars, {sequences.count('>')} sequences")
                
                # Pattern 2: Extract simple sequences without headers
                if not sequences:
                    simple_match = re.search(r'sequences?[:\s]+([ATCG\s\n]+)', command, re.IGNORECASE)
                    if simple_match:
                        sequences = simple_match.group(1).strip()
                        # Convert to FASTA format
                        seq_list = sequences.split()
                        if len(seq_list) >= 2:
                            sequences = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                        print(f"🔧 Converted simple sequences to FASTA: {len(sequences)} chars, {len(seq_list)} sequences")
                
                # Pattern 3: Extract sequences after "align" or "perform" keywords
                if not sequences:
                    # Try pattern for "perform multiple sequence alignment" followed by sequences
                    perform_match = re.search(r'perform[^>]*(>[\w\s\n>ATCGUatcgu]+)', command, re.IGNORECASE | re.DOTALL)
                    if perform_match:
                        sequences = perform_match.group(1).strip()
                        # Clean up blank lines
                        lines = sequences.split('\n')
                        cleaned_lines = [line.strip() for line in lines if line.strip()]
                        sequences = '\n'.join(cleaned_lines)
                        print(f"🔧 Extracted sequences after 'perform': {len(sequences)} chars, {sequences.count('>')} sequences")
                    
                    # Also try pattern for "align" keyword
                    if not sequences:
                        align_match = re.search(r'align[^:]*[:\s]+([ATCG\s\n>]+)', command, re.IGNORECASE | re.DOTALL)
                        if align_match:
                            sequences = align_match.group(1).strip()
                            # If not in FASTA format, convert it
                            if not sequences.startswith('>'):
                                seq_list = sequences.split()
                                if len(seq_list) >= 2:
                                    sequences = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                            else:
                                # Clean up blank lines in FASTA format
                                lines = sequences.split('\n')
                                cleaned_lines = [line.strip() for line in lines if line.strip()]
                                sequences = '\n'.join(cleaned_lines)
                            print(f"🔧 Extracted sequences after align: {len(sequences)} chars, {sequences.count('>')} sequences")
                
                # Pattern 4: Extract sequences with target/reference format
                if not sequences:
                    target_ref_match = re.search(r'>target\s+([ATCG]+)\s+>reference\s+([ATCG]+)', command, re.IGNORECASE)
                    if target_ref_match:
                        seq1 = target_ref_match.group(1)
                        seq2 = target_ref_match.group(2)
                        sequences = f">target\n{seq1}\n>reference\n{seq2}"
                        print(f"🔧 Extracted target/reference sequences: {len(sequences)} chars, 2 sequences")
                
                # Pattern 5: Extract inline FASTA format and convert to proper format
                if not sequences:
                    inline_fasta_match = re.search(r'sequences?[:\s]+(>[\w]+\s+[ATCG]+(?:\s+>[\w]+\s+[ATCG]+)*)', command, re.IGNORECASE)
                    if inline_fasta_match:
                        inline_fasta = inline_fasta_match.group(1)
                        # Convert inline format to proper FASTA format
                        # Split by > and process each sequence
                        parts = inline_fasta.split('>')
                        sequences = []
                        for part in parts:
                            if part.strip():
                                # Split on whitespace to separate name and sequence
                                name_seq = part.strip().split(None, 1)
                                if len(name_seq) == 2:
                                    name, seq = name_seq
                                    sequences.append(f">{name}\n{seq}")
                        sequences = "\n".join(sequences)
                        print(f"🔧 Converted inline FASTA to proper format: {len(sequences)} chars, {len(parts)} sequences")
                
                if sequences:
                    return {"sequences": sequences}
                else:
                    # Use longer default sequences for better demo
                    return {"sequences": ">seq1\nATGCGATCGATCGATCG\n>seq2\nATGCGATCGATCGATCG"}
        
        elif tool_name == "select_variants":
            # Extract selection parameters
            selection_criteria = "diversity"
            if "best" in command.lower() or "conservation" in command.lower():
                selection_criteria = "conservation"
            elif "random" in command.lower():
                selection_criteria = "random"
            elif "length" in command.lower():
                selection_criteria = "length"
            
            num_match = re.search(r'(\d+)\s+(variants?|sequences?)', command)
            num_variants = int(num_match.group(1)) if num_match else 5
            
            # Get session ID from session context
            session_id = session_context.get("session_id", "default")
            
            return {
                "session_id": session_id,
                "selection_criteria": selection_criteria,
                "num_variants": num_variants,
                "custom_filters": None
            }
        
        elif tool_name == "phylogenetic_tree":
            # First try to extract sequences from the command
            sequences_text = None
            
            # Pattern 1: Extract sequences after "sequences:" or "File content:"
            sequences_match = re.search(r'(?:sequences?|File content?)[:\s]+([^"]+)', command, re.IGNORECASE | re.DOTALL)
            if sequences_match:
                sequences_text = sequences_match.group(1).strip()
                print(f"🔧 Extracted sequences for phylogenetic tree: {len(sequences_text)} chars, {sequences_text.count('>')} sequences")
            
            # Pattern 2: Extract FASTA content directly (if no sequences: prefix)
            if not sequences_text:
                fasta_match = re.search(r'(>[\w\s\n>]+)', command, re.DOTALL)
                if fasta_match:
                    sequences_text = fasta_match.group(1).strip()
                    print(f"🔧 Extracted FASTA content directly: {len(sequences_text)} chars, {sequences_text.count('>')} sequences")
            
            if sequences_text:
                # Convert simple space-separated sequences to FASTA format
                if not sequences_text.startswith('>'):
                    seq_list = sequences_text.split()
                    if len(seq_list) >= 2:
                        sequences_text = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                        print(f"🔧 Converted to FASTA format: {len(sequences_text)} chars, {len(seq_list)} sequences")
                
                # Convert inline FASTA format to proper format
                elif '\n' not in sequences_text:
                    seq_count = sequences_text.count('>')
                    print(f"🔧 Converting inline FASTA format: {len(sequences_text)} chars, {seq_count} sequences")
                    # Convert inline format to proper FASTA format
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                    print(f"🔧 After conversion: {len(sequences_text)} chars")
                    # Clean up any remaining inline sequences
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                    print(f"🔧 Final result: {len(sequences_text)} chars, {sequences_text.count('>')} sequences")
                
                return {"aligned_sequences": sequences_text}
            
            # Use aligned sequences from session context
            aligned_sequences = session_context.get("aligned_sequences", "")
            
            # If not found in session context, try to get from previous alignment result in history
            if not aligned_sequences:
                history = session_context.get("history", [])
                for entry in reversed(history):  # Start from most recent
                    if entry.get("tool") == "sequence_alignment":
                        result = entry.get("result", {})
                        # Check for nested result structure
                        if isinstance(result, dict):
                            actual_result = result.get("result", result)
                            if isinstance(actual_result, dict):
                                # Extract alignment list and convert to FASTA
                                alignment = actual_result.get("alignment", [])
                                if isinstance(alignment, list) and len(alignment) > 0:
                                    fasta_lines = []
                                    for seq in alignment:
                                        if isinstance(seq, dict):
                                            name = seq.get("name", "sequence")
                                            sequence = seq.get("sequence", "")
                                            fasta_lines.append(f">{name}")
                                            fasta_lines.append(sequence)
                                    aligned_sequences = "\n".join(fasta_lines)
                                    print(f"🔧 Extracted aligned sequences from history: {len(alignment)} sequences")
                                    break  # Use the most recent alignment result
            
            if not aligned_sequences and session_context.get("mutated_sequences"):
                # If we have mutated sequences but no aligned sequences, create FASTA format
                aligned_sequences = "\n".join([f">mutant_{i+1}\n{seq}" for i, seq in enumerate(session_context["mutated_sequences"])])
            
            return {"aligned_sequences": aligned_sequences}
        
        elif tool_name == "clustering_analysis":
            # Extract number of clusters or representatives
            clusters_match = re.search(r'(\d+)\s+clusters?', command, re.IGNORECASE)
            representatives_match = re.search(r'(\d+)\s+representative\s+sequences?', command, re.IGNORECASE)
            
            if representatives_match:
                num_clusters = int(representatives_match.group(1))
            elif clusters_match:
                num_clusters = int(clusters_match.group(1))
            else:
                num_clusters = 5
            
            # Extract sequences (same logic as phylogenetic_tree)
            sequences_text = None
            
            # Pattern 1: Extract sequences after "sequences:" or "File content:"
            sequences_match = re.search(r'(?:sequences?|File content?)[:\s]+([^"]+)', command, re.IGNORECASE | re.DOTALL)
            if sequences_match:
                sequences_text = sequences_match.group(1).strip()
                print(f"🔧 Extracted sequences for clustering: {len(sequences_text)} chars, {sequences_text.count('>')} sequences")
            
            # Pattern 2: Extract FASTA content directly
            if not sequences_text:
                fasta_match = re.search(r'(>[\w\s\n>]+)', command, re.DOTALL)
                if fasta_match:
                    sequences_text = fasta_match.group(1).strip()
                    print(f"🔧 Extracted FASTA content for clustering: {len(sequences_text)} chars, {sequences_text.count('>')} sequences")
            
            if sequences_text:
                # Convert simple space-separated sequences to FASTA format
                if not sequences_text.startswith('>'):
                    seq_list = sequences_text.split()
                    if len(seq_list) >= 2:
                        sequences_text = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                        print(f"🔧 Converted to FASTA format for clustering: {len(sequences_text)} chars, {len(seq_list)} sequences")
                
                # Convert inline FASTA format to proper format
                elif '\n' not in sequences_text:
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                
                return {"aligned_sequences": sequences_text, "num_clusters": num_clusters}
            
            # Use aligned sequences from session context
            aligned_sequences = session_context.get("aligned_sequences", "")
            
            # If not found in session context, try to get from previous alignment result in history
            if not aligned_sequences:
                history = session_context.get("history", [])
                for entry in reversed(history):  # Start from most recent
                    if entry.get("tool") == "sequence_alignment":
                        result = entry.get("result", {})
                        # Check for nested result structure
                        if isinstance(result, dict):
                            actual_result = result.get("result", result)
                            if isinstance(actual_result, dict):
                                # Extract alignment list and convert to FASTA
                                alignment = actual_result.get("alignment", [])
                                if isinstance(alignment, list) and len(alignment) > 0:
                                    fasta_lines = []
                                    for seq in alignment:
                                        if isinstance(seq, dict):
                                            name = seq.get("name", "sequence")
                                            sequence = seq.get("sequence", "")
                                            fasta_lines.append(f">{name}")
                                            fasta_lines.append(sequence)
                                    aligned_sequences = "\n".join(fasta_lines)
                                    print(f"🔧 Extracted aligned sequences from history for clustering: {len(alignment)} sequences")
                                    break  # Use the most recent alignment result
            
            if not aligned_sequences and session_context.get("mutated_sequences"):
                aligned_sequences = "\n".join([f">mutant_{i+1}\n{seq}" for i, seq in enumerate(session_context["mutated_sequences"])])
            
            # If no sequences found in session context, use default dataset
            if not aligned_sequences:
                print("🔧 No sequences found in session context, using default phylogenetic dataset")
                # Use a comprehensive default dataset for clustering
                aligned_sequences = ">Random_Sequence_01\nACTCGATCACAAAGCTTAGGTCCGATCAATTTTGATAGTTACCCCCCACGGTCCAATCCGTTGGGTGAACACCGAGAAATTCGACAGATTTGCACTGCAAGTGCAGTCAGTAGGAGTTGCTGACTTACGGGCCGGGATGTCGTACGTCCACGG\n>Random_Sequence_02\nGTGCCGAACTAAGGAGACGTTACAGTACGCACCAGCAGACTCTCACAAAGACTCTGGCTAGTCCGTCGAAACGGCCTGCTAGAACAATGAAAGAGCCACGTCAAAAGAAAACTTCGTTGTACCTAGCGTCAGGTTTCTGCTAGAAACAGCAAGATCGCAGTCGTATGATTGATGGGGTACTCAGCC\n>Random_Sequence_03\nTTCTCACACTGTGTAAAAATTACACAAAAGATACGCCCAGTATTGGGGTGGGTATCCTCCGGGATGGGTAACTGGGGGTTCCCTTATGGTCAATGGAAAACCAGCCAAGATACATCTCATTGTTATAGGATGTTGAGCGCCATTAGCCTGCGATCACTGGGCGCCGTTTTTTCAACGTTTCTCCTCAC\n>Random_Sequence_04\nTTTGTGTGTTCACCTGGTGTCCAACAATTCGATGGATCATTGGGCCGATCCGTTAGCGCCGAACGCGAGTGTTGGGAGTTTTCTGCCGTCGACGTCGTCGAGTGAAATATCAAGCCCTGCAGGTCGACTGCGGCGTGTTGACCGTTAGTGGTTTACAATGGCTGTTAACGTTAATTGCAGGTACCCTGCAG\n>Random_Sequence_05\nTAAATGACGTCAGACTCTCTTAGTTATGCTCCGACTGGCTTTTACAGTTTCTTATAATAGGCTAGCAGCAAGAGGGTCCCGGTGTCCGTTTGGTTATCCTCGTTCTCAGTTGGTGAATAGGGACGCGGCATATGTACGGCAACGTATAATGT\n>Random_Sequence_06\nCGTATTCCAATTCGTCGGATCAAAGTGAGCTACCGAAGCCTGGAGTAGTTGGTTCACACAACGCCCTATCCGGTACCAGCAGGGAGGTTCCTGGGTCCACCTGAACTGGAAGCTATTTTCCAGCCATTTCCAGCCCCTGCCACCTCAAGATAAACACAAACAGGTTACTGTTAAGCACG\n>Random_Sequence_07\nCGGACGTGTAATTCCGAGCCCTCGCGTCACTCACAGGTGGATTCTGAACGGACGTCACAGTACGGAGGCCACAGTGGCGCAGCGTTGCCTGGACGAGGTCAACGGGGGATGTTTGCCGCAGAAAATAATGCAAGAGACAGTCACTTTAACATCGAAAATCTATTATCTAACCGCGATGACCG\n>Random_Sequence_08\nTGATTGATCAACAACCCTCACCCTATCTAATCAAAACCCCTGGGAAAAAGCCCTTCGGCTGAACTTGTCAATAACGGCAATCCAGGCTCGTACCACAAAGTAACGGGCTCCTTGTGGTCCCCTTCCTTGGACTCATTTAGCACGATTCGTTTCTTCCGTCAACCCCT\n>Random_Sequence_09\nAATTAGCTCATTGCGTAAGAGCAGGTTTCGCTGCCCATTCGTTCGGACAGGTACACTTGGAGGAGTTGCCTGTCGGAAAGGGTGGAACCAGCCGCTTGACGAGCTTTTGTTCTAGTAACGAGCTGTGCGGTTTTTGGGCCAAGGTCTCAGCATATACCTGTGCAAGTCAA\n>Random_Sequence_10\nTTAACTTAAATGCGCACTCTTCAGAGAGAGGCTCTGCCTTGAAACTCTGCGTGCATAGATTCCAGGCGAACATGGTATATGTGTCACGTCACAAACTACCATCCGCAGGCCGGCTTCTGTGAAAGCTATAATAACTGCGAGCTATCCGTATA"
            
            return {"aligned_sequences": aligned_sequences, "num_clusters": num_clusters}
        
        elif tool_name == "variant_selection":
            # Extract number of variants
            variants_match = re.search(r'(\d+)\s+variants?', command, re.IGNORECASE)
            num_variants = int(variants_match.group(1)) if variants_match else 10
            
            # Extract sequences (same logic as phylogenetic_tree)
            sequences_text = None
            
            # Pattern 1: Extract sequences after "sequences:" or "File content:"
            sequences_match = re.search(r'(?:sequences?|File content?)[:\s]+([^"]+)', command, re.IGNORECASE | re.DOTALL)
            if sequences_match:
                sequences_text = sequences_match.group(1).strip()
                print(f"🔧 Extracted sequences for variant selection: {len(sequences_text)} chars, {sequences_text.count('>')} sequences")
            
            # Pattern 2: Extract FASTA content directly
            if not sequences_text:
                fasta_match = re.search(r'(>[\w\s\n>]+)', command, re.DOTALL)
                if fasta_match:
                    sequences_text = fasta_match.group(1).strip()
                    print(f"🔧 Extracted FASTA content for variant selection: {len(sequences_text)} chars, {sequences_text.count('>')} sequences")
            
            if sequences_text:
                # Convert simple space-separated sequences to FASTA format
                if not sequences_text.startswith('>'):
                    seq_list = sequences_text.split()
                    if len(seq_list) >= 2:
                        sequences_text = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                        print(f"🔧 Converted to FASTA format for variant selection: {len(sequences_text)} chars, {len(seq_list)} sequences")
                
                # Convert inline FASTA format to proper format
                elif '\n' not in sequences_text:
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                
                return {"aligned_sequences": sequences_text, "num_variants": num_variants}
            
            # Use aligned sequences from session context
            aligned_sequences = session_context.get("aligned_sequences", "")
            if not aligned_sequences and session_context.get("mutated_sequences"):
                aligned_sequences = "\n".join([f">mutant_{i+1}\n{seq}" for i, seq in enumerate(session_context["mutated_sequences"])])
            
            return {"aligned_sequences": aligned_sequences, "num_variants": num_variants}
        
        elif tool_name == "fastqc_quality_analysis":
            # Extract R1 and R2 paths from command
            # Look for S3 paths or file paths
            r1_path = None
            r2_path = None
            output_path = None
            
            # Pattern 1: "forward reads are available here: s3://..."
            r1_match = re.search(r'(?:forward|r1|read\s*1)[^:]*:\s*(s3://[^\s]+|/[^\s]+)', command, re.IGNORECASE)
            if r1_match:
                r1_path = r1_match.group(1).strip()
            
            # Pattern 2: "reverse reads are available here: s3://..."
            r2_match = re.search(r'(?:reverse|r2|read\s*2)[^:]*:\s*(s3://[^\s]+|/[^\s]+)', command, re.IGNORECASE)
            if r2_match:
                r2_path = r2_match.group(1).strip()
            
            # Pattern 3: Look for any S3 paths containing R1 or R2
            if not r1_path:
                r1_match = re.search(r's3://[^\s]*(?:R1|r1|_1\.fq|mate_R1)[^\s]*', command, re.IGNORECASE)
                if r1_match:
                    r1_path = r1_match.group(0).strip()
            
            if not r2_path:
                r2_match = re.search(r's3://[^\s]*(?:R2|r2|_2\.fq|mate_R2)[^\s]*', command, re.IGNORECASE)
                if r2_match:
                    r2_path = r2_match.group(0).strip()
            
            # Extract output path if mentioned
            output_match = re.search(r'(?:output|results|save)[^:]*:\s*(s3://[^\s]+)', command, re.IGNORECASE)
            if output_match:
                output_path = output_match.group(1).strip()
            
            return {
                "input_r1": r1_path or "",
                "input_r2": r2_path or "",
                "output": output_path
            }
        
        elif tool_name == "dna_vendor_research":
            return {
                "command": command,
                "sequence_length": 1000,
                "quantity": "large"
            }
        
        elif tool_name == "synthesis_submission":
            return {
                "sequences": session_context.get("selected_sequences", ""),
                "vendor_preference": None,
                "quantity": "standard",
                "delivery_time": "standard"
            }
        
        elif tool_name == "plasmid_visualization":
            # Extract sequence from command
            sequence_match = re.search(r'sequence\s+([ATCG\s]+)', command, re.IGNORECASE)
            sequence = sequence_match.group(1).replace(" ", "").upper() if sequence_match else None
            
            # Check if this is a full plasmid (no "into" or "insert" keywords) or an insert
            is_full_plasmid = not re.search(r'\b(?:into|insert|in)\s+', command, re.IGNORECASE)
            
            # Extract vector name if specified
            vector_match = re.search(r'\b(?:into|in)\s+(\w+)', command, re.IGNORECASE)
            vector_name = vector_match.group(1) if vector_match else None
            
            # Extract position if specified
            position_match = re.search(r'\bat\s+position\s+(\d+)', command, re.IGNORECASE)
            insert_position = int(position_match.group(1)) if position_match else None
            
            if is_full_plasmid and sequence:
                # Full plasmid sequence provided
                return {
                    "full_plasmid_sequence": sequence,
                    "vector_name": None,
                    "cloning_sites": "",
                    "insert_sequence": ""
                }
            else:
                # Insert sequence into vector
                if not sequence:
                    sequence = session_context.get("selected_sequences", ["ATGCGATCG"])[0] if session_context.get("selected_sequences") else "ATGCGATCG"
                
                return {
                    "vector_name": vector_name or "pUC19",
                    "cloning_sites": "EcoRI, BamHI, HindIII",
                    "insert_sequence": sequence,
                    "insert_position": insert_position
                }
        
        elif tool_name == "plasmid_for_representatives":
            # Extract representatives from session context (from clustering results)
            representatives = session_context.get("clustering_result", {}).get("representatives", [])
            aligned_sequences = session_context.get("aligned_sequences", "")
            
            # Extract vector name from command
            vector_name = "pUC19"  # Default
            vector_match = re.search(r'(\w+)\s+vector', command, re.IGNORECASE)
            if vector_match:
                vector_name = vector_match.group(1)
            
            return {
                "representatives": representatives,
                "aligned_sequences": aligned_sequences,
                "vector_name": vector_name,
                "cloning_sites": "EcoRI, BamHI, HindIII"
            }
        
        elif tool_name == "directed_evolution":
            # Extract directed evolution parameters
            command_lower = command.lower()
            target_property = "thermal_stability"
            if "activity" in command_lower:
                target_property = "activity"
            elif "expression" in command_lower:
                target_property = "expression"
            
            library_size = 50
            size_match = re.search(r'(\d+)\s+(mutants?|variants?)', command)
            if size_match:
                library_size = int(size_match.group(1))
            
            strategy = "rational"
            if "random" in command_lower:
                strategy = "random"
            
            num_cycles = 1
            cycles_match = re.search(r'(\d+)\s+cycles?', command)
            if cycles_match:
                num_cycles = int(cycles_match.group(1))
            
            return {
                "target_property": target_property,
                "library_size": library_size,
                "strategy": strategy,
                "num_cycles": num_cycles,
                "command": command
            }
        
        elif tool_name == "read_trimming":
            # Extract FASTQ reads from command
            # Split command into sections and extract FASTQ content
            forward_reads = ""
            reverse_reads = ""
            reads = ""
            
            # Pattern 1: Extract from "Forward reads (filename):\n{content}"
            # Use a simpler approach: find the section header and extract everything until next section or end
            forward_section = re.search(r'Forward reads\s*\([^)]+\):\s*\n', command, re.IGNORECASE)
            if forward_section:
                start_pos = forward_section.end()
                # Find the next section or end of string
                next_section = re.search(r'\n\n(?:Reverse reads|Forward reads|File content)', command[start_pos:], re.IGNORECASE)
                if next_section:
                    forward_reads = command[start_pos:start_pos + next_section.start()].strip()
                else:
                    forward_reads = command[start_pos:].strip()
                print(f"🔧 [DEBUG] Extracted forward reads: {len(forward_reads)} characters")
            
            reverse_section = re.search(r'Reverse reads\s*\([^)]+\):\s*\n', command, re.IGNORECASE)
            if reverse_section:
                start_pos = reverse_section.end()
                # Find the next section or end of string
                next_section = re.search(r'\n\n(?:Forward reads|Reverse reads|File content)', command[start_pos:], re.IGNORECASE)
                if next_section:
                    reverse_reads = command[start_pos:start_pos + next_section.start()].strip()
                else:
                    reverse_reads = command[start_pos:].strip()
                print(f"🔧 [DEBUG] Extracted reverse reads: {len(reverse_reads)} characters")
            
            if not forward_reads and not reverse_reads:
                # Pattern 2: Extract FASTQ content directly (look for @ headers)
                # Match FASTQ records: @header\nsequence\n+\nquality (can repeat)
                fastq_match = re.search(r'(@[^\n]+\n[^\n]+\n\+[^\n]+\n[^\n]+(?:\n@[^\n]+\n[^\n]+\n\+[^\n]+\n[^\n]+)*)', command, re.MULTILINE)
                if fastq_match:
                    reads = fastq_match.group(0).strip()
                    print(f"🔧 [DEBUG] Extracted combined reads: {len(reads)} characters")
                else:
                    # Try to get trimmed reads from session history (for adapter removal on already-trimmed reads)
                    history = session_context.get("history", [])
                    for entry in reversed(history):  # Start from most recent
                        if entry.get("tool") == "read_trimming":
                            result = entry.get("result", {})
                            # Check for nested result structure
                            if isinstance(result, dict):
                                actual_result = result.get("result", result)
                                if isinstance(actual_result, dict):
                                    # Check for paired-end results
                                    if "forward_reads" in actual_result and "reverse_reads" in actual_result:
                                        forward_data = actual_result["forward_reads"]
                                        reverse_data = actual_result["reverse_reads"]
                                        if isinstance(forward_data, dict) and "trimmed_reads" in forward_data:
                                            forward_reads = forward_data["trimmed_reads"]
                                        if isinstance(reverse_data, dict) and "trimmed_reads" in reverse_data:
                                            reverse_reads = reverse_data["trimmed_reads"]
                                        print(f"🔧 [DEBUG] Retrieved paired-end trimmed reads from session history")
                                        break
                                    # Check for single-end results
                                    elif "trimmed_reads" in actual_result:
                                        reads = actual_result["trimmed_reads"]
                                        print(f"🔧 [DEBUG] Retrieved trimmed reads from session history: {len(reads)} characters")
                                        break
                    
                    # Also check results dict
                    if not forward_reads and not reverse_reads and not reads:
                        results = session_context.get("results", {})
                        trimming_keys = [k for k in results.keys() if k.startswith("read_trimming_")]
                        if trimming_keys:
                            latest_key = sorted(trimming_keys, key=lambda x: int(x.split("_")[-1]) if x.split("_")[-1].isdigit() else 0)[-1]
                            trimmed_result = results.get(latest_key, {})
                            if isinstance(trimmed_result, dict):
                                actual_result = trimmed_result.get("result", trimmed_result)
                                if isinstance(actual_result, dict):
                                    if "forward_reads" in actual_result and "reverse_reads" in actual_result:
                                        forward_data = actual_result["forward_reads"]
                                        reverse_data = actual_result["reverse_reads"]
                                        if isinstance(forward_data, dict) and "trimmed_reads" in forward_data:
                                            forward_reads = forward_data["trimmed_reads"]
                                        if isinstance(reverse_data, dict) and "trimmed_reads" in reverse_data:
                                            reverse_reads = reverse_data["trimmed_reads"]
                                    elif "trimmed_reads" in actual_result:
                                        reads = actual_result["trimmed_reads"]
            
            # Extract quality threshold
            quality_threshold = 20  # default
            threshold_match = re.search(r'quality\s+threshold\s+(\d+)', command, re.IGNORECASE)
            if threshold_match:
                quality_threshold = int(threshold_match.group(1))
            else:
                # Look for "quality 20" or "Q20" patterns
                simple_threshold = re.search(r'quality\s+(\d+)', command, re.IGNORECASE)
                if simple_threshold:
                    quality_threshold = int(simple_threshold.group(1))
            
            # Extract adapter sequence
            adapter = None
            # Try multiple patterns to catch different command formats
            # Pattern 1: "remove adapter sequences AGATCGGAAGAGC" or "remove adapter AGATCGGAAGAGC"
            remove_adapter = re.search(r'remove\s+adapter\s+(?:sequences?|sequence)?\s+([ATCGN]+)', command, re.IGNORECASE)
            if remove_adapter:
                adapter = remove_adapter.group(1).strip()
                print(f"🔧 [DEBUG] Extracted adapter from 'remove adapter' pattern: {adapter}")
            else:
                # Pattern 2: "adapter sequences AGATCGGAAGAGC" or "adapter AGATCGGAAGAGC"
                adapter_match = re.search(r'adapter\s+(?:sequences?|sequence)?\s+([ATCGN]+)', command, re.IGNORECASE)
                if adapter_match:
                    adapter = adapter_match.group(1).strip()
                    print(f"🔧 [DEBUG] Extracted adapter from 'adapter' pattern: {adapter}")
            
            if adapter:
                print(f"🔧 [DEBUG] Adapter sequence to remove: {adapter}")
            else:
                print(f"🔧 [DEBUG] No adapter sequence found in command")
            
            # Return parameters - prefer forward/reverse if both present
            if forward_reads and reverse_reads:
                return {
                    "forward_reads": forward_reads,
                    "reverse_reads": reverse_reads,
                    "adapter": adapter,
                    "quality_threshold": quality_threshold
                }
            elif reads:
                return {
                    "reads": reads,
                    "adapter": adapter,
                    "quality_threshold": quality_threshold
                }
            else:
                return {
                    "forward_reads": forward_reads,
                    "reverse_reads": reverse_reads,
                    "adapter": adapter,
                    "quality_threshold": quality_threshold
                }
        
        elif tool_name == "read_merging":
            # Extract forward and reverse reads
            # Pattern 1: Extract S3 paths or file paths from "R1: s3://..." or "R2: s3://..."
            forward_reads = ""
            reverse_reads = ""
            
            # Pattern 1a: Extract S3 paths from "R1: s3://..." format
            r1_match = re.search(r'R1\s*:\s*(s3://[^\s]+|/[^\s]+)', command, re.IGNORECASE)
            if r1_match:
                forward_reads = r1_match.group(1).strip()
            
            r2_match = re.search(r'R2\s*:\s*(s3://[^\s]+|/[^\s]+)', command, re.IGNORECASE)
            if r2_match:
                reverse_reads = r2_match.group(1).strip()
            
            # Pattern 1b: Extract S3 paths containing R1 or R2 in the path
            if not forward_reads:
                r1_path_match = re.search(r's3://[^\s]*(?:R1|r1|_1\.fq|mate_R1)[^\s]*', command, re.IGNORECASE)
                if r1_path_match:
                    forward_reads = r1_path_match.group(0).strip()
            
            if not reverse_reads:
                r2_path_match = re.search(r's3://[^\s]*(?:R2|r2|_2\.fq|mate_R2)[^\s]*', command, re.IGNORECASE)
                if r2_path_match:
                    reverse_reads = r2_path_match.group(0).strip()
            
            # Pattern 2: Extract from "Forward reads (filename):\n{content}" and "Reverse reads (filename):\n{content}"
            if not forward_reads:
                forward_section = re.search(r'Forward reads\s*\([^)]+\):\s*\n', command, re.IGNORECASE)
                if forward_section:
                    start_pos = forward_section.end()
                    next_section = re.search(r'\n\n(?:Reverse reads|Forward reads|File content)', command[start_pos:], re.IGNORECASE)
                    if next_section:
                        forward_reads = command[start_pos:start_pos + next_section.start()].strip()
                    else:
                        forward_reads = command[start_pos:].strip()
            
            if not reverse_reads:
                reverse_section = re.search(r'Reverse reads\s*\([^)]+\):\s*\n', command, re.IGNORECASE)
                if reverse_section:
                    start_pos = reverse_section.end()
                    next_section = re.search(r'\n\n(?:Forward reads|Reverse reads|File content)', command[start_pos:], re.IGNORECASE)
                    if next_section:
                        reverse_reads = command[start_pos:start_pos + next_section.start()].strip()
                    else:
                        reverse_reads = command[start_pos:].strip()
            
            # If not found in command, try to get from session context (trimmed reads from previous step)
            if not forward_reads or not reverse_reads:
                # Try to get trimmed reads from session history
                # Look for the latest read_trimming result in history
                history = session_context.get("history", [])
                for entry in reversed(history):  # Start from most recent
                    if entry.get("tool") == "read_trimming":
                        result = entry.get("result", {})
                        # Check for nested result structure
                        if isinstance(result, dict):
                            actual_result = result.get("result", result)
                            if isinstance(actual_result, dict):
                                if "forward_reads" in actual_result:
                                    forward_data = actual_result["forward_reads"]
                                    if isinstance(forward_data, dict) and "trimmed_reads" in forward_data:
                                        forward_reads = forward_data["trimmed_reads"]
                                if "reverse_reads" in actual_result:
                                    reverse_data = actual_result["reverse_reads"]
                                    if isinstance(reverse_data, dict) and "trimmed_reads" in reverse_data:
                                        reverse_reads = reverse_data["trimmed_reads"]
                        break  # Use the most recent trimming result
                
                # Also check results dict
                if not forward_reads or not reverse_reads:
                    results = session_context.get("results", {})
                    # Find the latest read_trimming result
                    trimming_keys = [k for k in results.keys() if k.startswith("read_trimming_")]
                    if trimming_keys:
                        # Sort by number to get the latest
                        latest_key = sorted(trimming_keys, key=lambda x: int(x.split("_")[-1]) if x.split("_")[-1].isdigit() else 0)[-1]
                        trimmed_result = results.get(latest_key, {})
                        if isinstance(trimmed_result, dict):
                            actual_result = trimmed_result.get("result", trimmed_result)
                            if isinstance(actual_result, dict):
                                if "forward_reads" in actual_result:
                                    forward_data = actual_result["forward_reads"]
                                    if isinstance(forward_data, dict) and "trimmed_reads" in forward_data:
                                        forward_reads = forward_data["trimmed_reads"]
                                if "reverse_reads" in actual_result:
                                    reverse_data = actual_result["reverse_reads"]
                                    if isinstance(reverse_data, dict) and "trimmed_reads" in reverse_data:
                                        reverse_reads = reverse_data["trimmed_reads"]
            
            # Extract minimum overlap
            min_overlap = 12  # default
            overlap_match = re.search(r'overlap\s+of\s+(\d+)', command, re.IGNORECASE)
            if overlap_match:
                min_overlap = int(overlap_match.group(1))
            else:
                # Look for "minimum overlap 12" pattern
                min_overlap_match = re.search(r'minimum\s+overlap\s+(\d+)', command, re.IGNORECASE)
                if min_overlap_match:
                    min_overlap = int(min_overlap_match.group(1))
            
            return {
                "forward_reads": forward_reads,
                "reverse_reads": reverse_reads,
                "min_overlap": min_overlap,
                "command": command,  # Include original command for tool-generator-agent
                "original_command": command  # Also include as original_command
            }
        
        elif tool_name == "quality_assessment":
            # Extract merged sequences from command or session context
            merged_sequences = ""
            
            # Try to extract FASTA sequences from command
            fasta_match = re.search(r'>(?:merged_|sequence_)[^\n]+\n[ATCGN\n]+', command, re.IGNORECASE | re.MULTILINE)
            if fasta_match:
                merged_sequences = fasta_match.group(0).strip()
            
            # If not found in command, try to get from session context (merged sequences from previous step)
            if not merged_sequences:
                # Look for the latest read_merging result in history
                history = session_context.get("history", [])
                for entry in reversed(history):  # Start from most recent
                    if entry.get("tool") == "read_merging":
                        result = entry.get("result", {})
                        # Check for nested result structure
                        if isinstance(result, dict):
                            actual_result = result.get("result", result)
                            if isinstance(actual_result, dict):
                                if "merged_sequences" in actual_result:
                                    merged_sequences = actual_result["merged_sequences"]
                                elif "merged_reads" in actual_result:
                                    merged_sequences = actual_result["merged_reads"]
                        break  # Use the most recent merging result
                
                # Also check results dict
                if not merged_sequences:
                    results = session_context.get("results", {})
                    # Find the latest read_merging result
                    merging_keys = [k for k in results.keys() if k.startswith("read_merging_")]
                    if merging_keys:
                        # Sort by number to get the latest
                        latest_key = sorted(merging_keys, key=lambda x: int(x.split("_")[-1]) if x.split("_")[-1].isdigit() else 0)[-1]
                        merged_result = results.get(latest_key, {})
                        if isinstance(merged_result, dict):
                            actual_result = merged_result.get("result", merged_result)
                            if isinstance(actual_result, dict):
                                if "merged_sequences" in actual_result:
                                    merged_sequences = actual_result["merged_sequences"]
                                elif "merged_reads" in actual_result:
                                    merged_sequences = actual_result["merged_reads"]
            
            print(f"🔧 [DEBUG] Quality assessment: extracted {len(merged_sequences)} characters of merged sequences")
            
            return {
                "sequences": merged_sequences
            }
        
        elif tool_name == "single_cell_analysis":
            # Extract single-cell analysis parameters
            # Determine analysis steps
            steps = []
            if "marker" in command_lower or "markers" in command_lower:
                steps.append("markers")
            if "differential" in command_lower or "deg" in command_lower:
                steps.append("differential")
            if "pathway" in command_lower or "enrichment" in command_lower:
                steps.append("pathways")
            if "annotate" in command_lower or "cell type" in command_lower:
                steps.append("annotation")
            if "batch" in command_lower or "correct" in command_lower:
                steps.append("batch_correction")
            if not steps or "all" in command_lower or "complete" in command_lower or "full" in command_lower:
                steps = ["all"]
            
            # Extract data file from command or session context
            data_file = None
            data_format = "10x"  # default
            
            # Check for file references in command
            file_match = re.search(r'(?:file|data|input)[:\s]+([^\s]+\.(?:h5|rds|csv|h5ad|mtx))', command, re.IGNORECASE)
            if file_match:
                data_file = file_match.group(1)
                if data_file.endswith('.h5') or data_file.endswith('.h5ad'):
                    data_format = "h5"
                elif data_file.endswith('.rds'):
                    data_format = "seurat"
                elif data_file.endswith('.csv'):
                    data_format = "csv"
            
            # Check session context for uploaded files
            if not data_file and session_context:
                uploaded_files = session_context.get("uploaded_files", [])
                for file_info in uploaded_files:
                    file_name = file_info.get("name", "")
                    if any(ext in file_name.lower() for ext in [".h5", ".h5ad", ".rds", "matrix.mtx", ".csv"]):
                        data_file = file_name
                        if file_name.endswith('.h5') or file_name.endswith('.h5ad'):
                            data_format = "h5"
                        elif file_name.endswith('.rds'):
                            data_format = "seurat"
                        elif file_name.endswith('.csv'):
                            data_format = "csv"
                        break
            
            # Extract additional parameters
            resolution = 0.5
            resolution_match = re.search(r'resolution[:\s]+([\d.]+)', command, re.IGNORECASE)
            if resolution_match:
                resolution = float(resolution_match.group(1))
            
            nfeatures = 2000
            nfeatures_match = re.search(r'(?:nfeatures|features|variable features)[:\s]+(\d+)', command, re.IGNORECASE)
            if nfeatures_match:
                nfeatures = int(nfeatures_match.group(1))
            
            return {
                "data_file": data_file,
                "data_format": data_format,
                "steps": steps,
                "resolution": resolution,
                "nfeatures": nfeatures,
                "command": command
            }
        
        elif tool_name == "fetch_ncbi_sequence":
            params = {}
            
            # Extract accession number (common patterns)
            accession_pattern = r'\b([A-Z]{1,2}_?\d+\.?\d*)\b'
            matches = re.findall(accession_pattern, command)
            if matches:
                params["accession"] = matches[0]
            
            # Determine database from command
            if any(word in command_lower for word in ['protein', 'prot', 'aa', 'amino']):
                params["database"] = "protein"
            else:
                params["database"] = "nucleotide"
            
            return params
        
        elif tool_name == "query_uniprot":
            # Use accession if present, otherwise use command text as query
            params = {}
            accession_match = re.search(r'\b[A-NR-Z][0-9]{5}\b|\b[OPQ][0-9][A-Z0-9]{3}[0-9]\b', command, re.IGNORECASE)
            if accession_match:
                params["query"] = accession_match.group(0)
            else:
                # Strip common leading phrases
                cleaned = re.sub(r'\b(query|search|get|lookup|find)\b', '', command_lower, flags=re.IGNORECASE).strip()
                cleaned = re.sub(r'\buniprot\b', '', cleaned, flags=re.IGNORECASE).strip()
                params["query"] = cleaned or command
            params["format"] = "fasta"
            params["limit"] = 10
            return params
        
        elif tool_name == "lookup_go_term":
            params = {}
            go_match = re.search(r'GO:\d{7}', command, re.IGNORECASE)
            if go_match:
                params["go_id"] = go_match.group(0).upper()
            else:
                params["go_id"] = command.strip()
            return params
        
        elif tool_name == "bulk_rnaseq_analysis":
            params = {}
            
            # Extract file paths (first two CSV-like tokens)
            csv_paths = re.findall(r'[\w./-]+\.csv', command)
            if len(csv_paths) >= 1:
                params["count_matrix"] = csv_paths[0]
            if len(csv_paths) >= 2:
                params["sample_metadata"] = csv_paths[1]
            
            # Fallback to session context if available
            if not params.get("count_matrix"):
                params["count_matrix"] = session_context.get("count_matrix", "")
            if not params.get("sample_metadata"):
                params["sample_metadata"] = session_context.get("sample_metadata", "")
            
            # Extract design formula and alpha
            design_match = re.search(r'design\s*[:=]\s*([~\w+\s]+)', command, re.IGNORECASE)
            if design_match:
                params["design_formula"] = design_match.group(1).strip()
            else:
                params["design_formula"] = "~condition"
            
            alpha_match = re.search(r'alpha\s*[:=]\s*([\d.]+)', command, re.IGNORECASE)
            params["alpha"] = float(alpha_match.group(1)) if alpha_match else 0.05
            
            return params
        
        else:
            return {"command": command} 
