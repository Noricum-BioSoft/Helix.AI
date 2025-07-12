import re
from typing import Dict, Any, Optional, Tuple
import sys
import os

# Add the project root to Python path to access tools
project_root = os.path.join(os.path.dirname(__file__), '..')
sys.path.insert(0, project_root)

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
            }
        }
    
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
        
        # Check for sequence alignment
        if any(phrase in command_lower for phrase in ['align', 'alignment', 'compare sequences', 'multiple sequence alignment']):
            print(f"🔧 Command router: Matched 'alignment' -> sequence_alignment")
            return 'sequence_alignment', self._extract_parameters(command, 'sequence_alignment', session_context)
        
        # Check for mutation/variant generation
        if any(phrase in command_lower for phrase in ['mutate', 'mutation', 'variant', 'create variants', 'generate variants']):
            print(f"🔧 Command router: Matched 'mutation' -> mutate_sequence")
            return 'mutate_sequence', self._extract_parameters(command, 'mutate_sequence', session_context)
        
        # Check for sequence selection (pick, select, choose)
        if any(phrase in command_lower for phrase in ['pick', 'select', 'choose']) and not any(phrase in command_lower for phrase in ['phylogenetic', 'evolutionary']):
            print(f"🔧 Command router: Matched 'selection' -> sequence_selection")
            return 'sequence_selection', self._extract_parameters(command, 'sequence_selection', session_context)
        
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
        
        # Default fallback
        print(f"🔧 Command router: No specific match, using general command execution")
        return "general_command", {"command": command}
    
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
                sequences = "\n".join([f">seq_{i+1}\n{seq}" for i, seq in enumerate(session_context["mutated_sequences"])])
                return {"sequences": sequences}
            else:
                # Extract sequences from command
                sequences_match = re.search(r'sequences?\s+([ATCG\s\n]+)', command, re.IGNORECASE)
                if sequences_match:
                    return {"sequences": sequences_match.group(1)}
                else:
                    return {"sequences": ">seq1\nATGCGATCG\n>seq2\nATGCGATC"}
        
        elif tool_name == "sequence_selection":
            # Extract selection parameters
            selection_type = "random"
            if "best" in command.lower() or "conservation" in command.lower():
                selection_type = "best_conservation"
            elif "gap" in command.lower():
                selection_type = "lowest_gaps"
            elif "gc" in command.lower():
                selection_type = "highest_gc"
            
            num_match = re.search(r'(\d+)\s+sequences?', command)
            num_sequences = int(num_match.group(1)) if num_match else 1
            
            # Use aligned sequences from session context
            aligned_sequences = session_context.get("aligned_sequences", "")
            
            return {
                "aligned_sequences": aligned_sequences,
                "selection_type": selection_type,
                "num_sequences": num_sequences
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
        
        else:
            return {"command": command} 