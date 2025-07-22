import re
from typing import Dict, Any, Optional, Tuple
import sys
import os

# Add the project root to Python path to access tools
project_root = os.path.join(os.path.dirname(__file__), '..')
sys.path.insert(0, project_root)

# Import directed evolution handler
from directed_evolution_handler import DirectedEvolutionHandler

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
                'keywords': ['order', 'vendor', 'synthesis', 'test', 'assay', 'expression', 'function', 'binding'],
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
        
        # Priority-based matching to avoid conflicts
        # Check for specific phrases first, then general keywords
        
        # Check for clustering analysis (HIGHEST PRIORITY - before variant selection)
        if any(phrase in command_lower for phrase in ['cluster', 'clustering', 'group sequences', 'cluster sequences', 'representative sequences']):
            print(f"🔧 Command router: Matched 'clustering' -> clustering_analysis")
            return 'clustering_analysis', self._extract_parameters(command, 'clustering_analysis', session_context)
        
        # Check for variant selection (HIGH PRIORITY - before phylogenetic tree)
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
        if any(phrase in command_lower for phrase in ['align sequences', 'sequence alignment', 'compare sequences', 'multiple sequence alignment', 'perform alignment']):
            print(f"🔧 Command router: Matched 'alignment' -> sequence_alignment")
            return 'sequence_alignment', self._extract_parameters(command, 'sequence_alignment', session_context)
        
        # Check for vendor research
        if any(phrase in command_lower for phrase in ['order', 'vendor', 'synthesis', 'test', 'assay', 'expression', 'function', 'binding']):
            print(f"🔧 Command router: Matched 'vendor' -> dna_vendor_research")
            return 'dna_vendor_research', self._extract_parameters(command, 'dna_vendor_research', session_context)
        
        # Check for plasmid visualization
        if any(phrase in command_lower for phrase in ['plasmid', 'vector', 'cloning']):
            print(f"🔧 Command router: Matched 'plasmid' -> plasmid_visualization")
            return 'plasmid_visualization', self._extract_parameters(command, 'plasmid_visualization', session_context)
        
        # Check for plasmid for representatives (HIGH PRIORITY - before general plasmid)
        if any(phrase in command_lower for phrase in ['insert representatives', 'express representatives', 'clone representatives', 'vector representatives', 'plasmid representatives']):
            print(f"🔧 Command router: Matched 'plasmid representatives' -> plasmid_for_representatives")
            return 'plasmid_for_representatives', self._extract_parameters(command, 'plasmid_for_representatives', session_context)
        
        # Check for synthesis submission
        if any(phrase in command_lower for phrase in ['submit', 'synthesis', 'order sequences']):
            print(f"🔧 Command router: Matched 'synthesis' -> synthesis_submission")
            return 'synthesis_submission', self._extract_parameters(command, 'synthesis_submission', session_context)
        
        # Check for session creation
        if any(phrase in command_lower for phrase in ['create session', 'new session', 'start session', 'initialize session']):
            print(f"🔧 Command router: Matched 'session creation' -> session_creation")
            return 'session_creation', {"command": command}
        
        # Check for directed evolution
        if any(phrase in command_lower for phrase in ['directed evolution', 'evolution', 'protein engineering', 'dbtl', 'design build test learn']):
            print(f"🔧 Command router: Matched 'directed evolution' -> directed_evolution")
            return 'directed_evolution', self._extract_parameters(command, 'directed_evolution', session_context)
        
        # Default fallback - try to route to sequence_alignment for alignment-like commands
        if any(phrase in command_lower for phrase in ['align', 'alignment', 'sequences']):
            print(f"🔧 Command router: Defaulting to sequence_alignment for alignment command")
            return "sequence_alignment", self._extract_parameters(command, 'sequence_alignment', session_context)
        
        # For other commands, try to use the natural command handler
        print(f"🔧 Command router: No specific match, using natural command handler")
        return "handle_natural_command", {"command": command, "session_id": session_context.get("session_id", "")}
    
    def _extract_parameters(self, command: str, tool_name: str, session_context: Dict[str, Any]) -> Dict[str, Any]:
        """Extract parameters for the specific tool."""
        
        if tool_name == "mutate_sequence":
            # Extract sequence and number of variants
            sequence_match = re.search(r'sequence\s+([ATCG]+)', command, re.IGNORECASE)
            variants_match = re.search(r'(\d+)\s+variants?', command)
            
            sequence = sequence_match.group(1) if sequence_match else "ATGCGATCG"
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
                print(f"[DEBUG] FASTA string sent to alignment tool:\n{sequences}")
                return {"sequences": sequences}
            else:
                # Extract sequences from command - improved pattern to handle FASTA format
                # Try multiple patterns to extract sequences
                sequences = None
                
                # Pattern 1: Extract FASTA format sequences
                fasta_match = re.search(r'sequences?[:\s]+(>[\w\s\n>]+)', command, re.IGNORECASE | re.DOTALL)
                if fasta_match:
                    sequences = fasta_match.group(1).strip()
                    print(f"🔧 Extracted FASTA sequences: {sequences}")
                
                # Pattern 2: Extract simple sequences without headers
                if not sequences:
                    simple_match = re.search(r'sequences?[:\s]+([ATCG\s\n]+)', command, re.IGNORECASE)
                    if simple_match:
                        sequences = simple_match.group(1).strip()
                        # Convert to FASTA format
                        seq_list = sequences.split()
                        if len(seq_list) >= 2:
                            sequences = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                        print(f"🔧 Converted simple sequences to FASTA: {sequences}")
                
                # Pattern 3: Extract sequences after "align" keyword
                if not sequences:
                    align_match = re.search(r'align[^:]*[:\s]+([ATCG\s\n>]+)', command, re.IGNORECASE | re.DOTALL)
                    if align_match:
                        sequences = align_match.group(1).strip()
                        # If not in FASTA format, convert it
                        if not sequences.startswith('>'):
                            seq_list = sequences.split()
                            if len(seq_list) >= 2:
                                sequences = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                        print(f"🔧 Extracted sequences after align: {sequences}")
                
                # Pattern 4: Extract sequences with target/reference format
                if not sequences:
                    target_ref_match = re.search(r'>target\s+([ATCG]+)\s+>reference\s+([ATCG]+)', command, re.IGNORECASE)
                    if target_ref_match:
                        seq1 = target_ref_match.group(1)
                        seq2 = target_ref_match.group(2)
                        sequences = f">target\n{seq1}\n>reference\n{seq2}"
                        print(f"🔧 Extracted target/reference sequences: {sequences}")
                
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
                        print(f"🔧 Converted inline FASTA to proper format: {sequences}")
                
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
                print(f"🔧 Extracted sequences for phylogenetic tree: '{sequences_text}'")
            
            # Pattern 2: Extract FASTA content directly (if no sequences: prefix)
            if not sequences_text:
                fasta_match = re.search(r'(>[\w\s\n>]+)', command, re.DOTALL)
                if fasta_match:
                    sequences_text = fasta_match.group(1).strip()
                    print(f"🔧 Extracted FASTA content directly: '{sequences_text}'")
            
            if sequences_text:
                # Convert simple space-separated sequences to FASTA format
                if not sequences_text.startswith('>'):
                    seq_list = sequences_text.split()
                    if len(seq_list) >= 2:
                        sequences_text = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                        print(f"🔧 Converted to FASTA format: '{sequences_text}'")
                
                # Convert inline FASTA format to proper format
                elif '\n' not in sequences_text:
                    print(f"🔧 Converting inline FASTA format: '{sequences_text}'")
                    # Convert inline format to proper FASTA format
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                    print(f"🔧 After conversion: '{sequences_text}'")
                    # Clean up any remaining inline sequences
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                    print(f"🔧 Final result: '{sequences_text}'")
                
                return {"aligned_sequences": sequences_text}
            
            # Use aligned sequences from session context
            aligned_sequences = session_context.get("aligned_sequences", "")
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
                print(f"🔧 Extracted sequences for clustering: '{sequences_text}'")
            
            # Pattern 2: Extract FASTA content directly
            if not sequences_text:
                fasta_match = re.search(r'(>[\w\s\n>]+)', command, re.DOTALL)
                if fasta_match:
                    sequences_text = fasta_match.group(1).strip()
                    print(f"🔧 Extracted FASTA content for clustering: '{sequences_text}'")
            
            if sequences_text:
                # Convert simple space-separated sequences to FASTA format
                if not sequences_text.startswith('>'):
                    seq_list = sequences_text.split()
                    if len(seq_list) >= 2:
                        sequences_text = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                        print(f"🔧 Converted to FASTA format for clustering: '{sequences_text}'")
                
                # Convert inline FASTA format to proper format
                elif '\n' not in sequences_text:
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                    sequences_text = re.sub(r'>([^\s]+)\s+([^>]+)', r'>\1\n\2', sequences_text)
                
                return {"aligned_sequences": sequences_text, "num_clusters": num_clusters}
            
            # Use aligned sequences from session context
            aligned_sequences = session_context.get("aligned_sequences", "")
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
                print(f"🔧 Extracted sequences for variant selection: '{sequences_text}'")
            
            # Pattern 2: Extract FASTA content directly
            if not sequences_text:
                fasta_match = re.search(r'(>[\w\s\n>]+)', command, re.DOTALL)
                if fasta_match:
                    sequences_text = fasta_match.group(1).strip()
                    print(f"🔧 Extracted FASTA content for variant selection: '{sequences_text}'")
            
            if sequences_text:
                # Convert simple space-separated sequences to FASTA format
                if not sequences_text.startswith('>'):
                    seq_list = sequences_text.split()
                    if len(seq_list) >= 2:
                        sequences_text = "\n".join([f">seq{i+1}\n{seq}" for i, seq in enumerate(seq_list)])
                        print(f"🔧 Converted to FASTA format for variant selection: '{sequences_text}'")
                
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
            return {
                "vector_name": "pUC19",
                "cloning_sites": "EcoRI, BamHI, HindIII",
                "insert_sequence": session_context.get("selected_sequences", ["ATGCGATCG"])[0] if session_context.get("selected_sequences") else "ATGCGATCG"
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
        
        else:
            return {"command": command} 