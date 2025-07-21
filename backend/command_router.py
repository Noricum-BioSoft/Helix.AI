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
        
        # Check for phylogenetic tree first (highest priority)
        if any(phrase in command_lower for phrase in ['phylogenetic tree', 'evolutionary tree', 'phylogeny']):
            print(f"🔧 Command router: Matched 'phylogenetic tree' -> phylogenetic_tree")
            return 'phylogenetic_tree', self._extract_parameters(command, 'phylogenetic_tree', session_context)
        
        # Check for sequence selection (pick, select, choose) - HIGH PRIORITY
        if any(phrase in command_lower for phrase in ['pick', 'select', 'choose']) and not any(phrase in command_lower for phrase in ['phylogenetic', 'evolutionary']):
            print(f"🔧 Command router: Matched 'selection' -> select_variants")
            return 'select_variants', self._extract_parameters(command, 'select_variants', session_context)
        
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
                # Extract sequences from command
                sequences_match = re.search(r'sequences?\s+([ATCG\s\n]+)', command, re.IGNORECASE)
                if sequences_match:
                    return {"sequences": sequences_match.group(1)}
                else:
                    return {"sequences": ">seq1\nATGCGATCG\n>seq2\nATGCGATC"}
        
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
            sequences_match = re.search(r'sequences?[:\s]+([^"]+)', command, re.IGNORECASE | re.DOTALL)
            if sequences_match:
                # Extract sequences from command
                sequences_text = sequences_match.group(1).strip()
                # Convert inline format to proper FASTA format
                # Handle both formats: ">seq1 ATCG >seq2 ATCG" and ">seq1\nATCG\n>seq2\nATCG"
                if '\n' not in sequences_text:
                    print(f"🔧 Converting inline format: '{sequences_text}'")
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