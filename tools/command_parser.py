import re
import logging
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class ParsedCommand:
    """Represents a parsed command with extracted parameters."""
    action: str
    tool: str
    parameters: Dict[str, Any]
    session_id: Optional[str] = None
    confidence: float = 0.0

class CommandParser:
    """Parses natural language commands into structured operations."""
    
    def __init__(self):
        self.action_patterns = {
            'select_variants': [
                r'select\s+(\d+)\s+variants?',
                r'pick\s+(\d+)\s+variants?',
                r'choose\s+(\d+)\s+variants?',
                r'from\s+.*variants?.*pick\s+(\d+)',
            ],
            'select_variants_diverse': [
                r'select\s+(\d+)\s+diverse\s+variants?',
                r'pick\s+(\d+)\s+diverse\s+variants?',
                r'choose\s+(\d+)\s+diverse\s+variants?',
            ],
            'select_variants_length': [
                r'select\s+(\d+)\s+sequences?\s+by\s+length',
                r'choose\s+(\d+)\s+sequences?\s+by\s+length',
                r'from\s+.*variants?.*pick\s+(\d+)\s+sequences?\s+by\s+length',
            ],
            'mutate_sequence': [
                r'generate\s+(\d+)\s+variants?',
                r'create\s+(\d+)\s+variants?',
                r'mutate\s+.*to\s+generate\s+(\d+)\s+variants?',
                r'generate\s+(\d+)\s+mutations?',
            ],
            'dna_vendor_research': [
                r'order\s+.*variants?.*from\s+.*vendor',
                r'find\s+.*synthesis\s+vendor',
                r'research\s+.*vendor',
                r'what\s+.*testing\s+options',
                r'find\s+.*testing\s+options',
                r'order\s+.*sequences?.*from\s+.*vendor',
                r'find\s+.*dna\s+synthesis',
                r'research\s+.*synthesis',
                r'what\s+.*vendors?.*can\s+.*synthesize',
                r'find\s+.*companies?.*can\s+.*make',
            ],
            'align_sequences': [
                r'align.*sequences?',
                r'perform\s+alignment',
                r'sequence\s+alignment',
            ],
            'analyze_data': [
                r'analyze\s+data',
                r'perform\s+analysis',
                r'data\s+analysis',
            ],
            'plasmid_visualization': [
                r'visualize\s+plasmid',
                r'show\s+plasmid',
                r'create\s+plasmid\s+visualization',
                r'plasmid\s+viewer',
                r'circular\s+plasmid',
                r'insert\s+.*sequence.*into\s+.*vector',
                r'express\s+.*sequence.*in\s+.*vector',
                r'clone\s+.*sequence.*into\s+.*vector',
                r'insert\s+.*into\s+.*vector',
                r'express\s+.*in\s+.*vector',
                r'clone\s+.*into\s+.*vector',
                r'put\s+.*sequence.*in\s+.*vector',
                r'place\s+.*sequence.*in\s+.*vector',
                r'insert\s+.*gene.*into\s+.*vector',
                r'express\s+.*gene.*in\s+.*vector',
            ],
            'phylogenetic_tree': [
                r'create\s+phylogenetic\s+tree',
                r'build\s+phylogenetic\s+tree',
                r'visualize\s+.*phylogenetic\s+tree',
                r'phylogenetic\s+tree',
                r'evolutionary\s+tree',
                r'tree\s+from\s+alignment',
            ],
            'sequence_selection': [
                r'pick\s+(\d+)\s+sequence',
                r'select\s+(\d+)\s+sequence',
                r'choose\s+(\d+)\s+sequence',
                r'pick\s+(\d+)\s+from\s+alignment',
                r'select\s+(\d+)\s+from\s+alignment',
                r'randomly\s+pick\s+(\d+)',
                r'randomly\s+select\s+(\d+)',
            ],
            'synthesis_submission': [
                r'submit\s+.*synthesis',
                r'order\s+.*synthesis',
                r'get\s+quote\s+for\s+synthesis',
                r'synthesis\s+quote',
                r'order\s+sequences',
                r'submit\s+sequences',
            ]
        }
        
        self.parameter_patterns = {
            'selection_criteria': {
                'random': [r'randomly', r'random'],
                'diversity': [r'diverse', r'diversity'],
                'length': [r'by\s+length', r'length\s+based'],
                'custom': [r'custom', r'specific']
            },
            'count_patterns': [
                r'(\d+)\s+sequences?',
                r'(\d+)\s+variants?',
                r'(\d+)\s+items?'
            ]
        }
    
    def parse_command(self, command: str, session_id: Optional[str] = None) -> ParsedCommand:
        """Parse a natural language command into structured parameters."""
        command_lower = command.lower().strip()
        logger.info(f"Parsing command: '{command}'")
        
        # Try to match action patterns
        for action, patterns in self.action_patterns.items():
            for pattern in patterns:
                match = re.search(pattern, command_lower, re.IGNORECASE)
                if match:
                    logger.info(f"Matched action: {action}")
                    return self._build_parsed_command(action, command, match, session_id)
        
        # If no specific action found, try to infer from context
        return self._infer_command(command, session_id)
    
    def _build_parsed_command(self, action: str, original_command: str, 
                             match: re.Match, session_id: Optional[str] = None) -> ParsedCommand:
        """Build a parsed command object from the matched action."""
        command_lower = original_command.lower()
        
        if action == 'select_variants':
            # Extract count
            count = self._extract_count(original_command)
            criteria = self._extract_selection_criteria(command_lower)
            
            return ParsedCommand(
                action="select_variants",
                tool="select_variants",
                parameters={
                    "num_variants": count,
                    "selection_criteria": criteria,
                    "session_id": session_id
                },
                session_id=session_id,
                confidence=0.9
            )
        
        elif action == 'select_variants_diverse':
            count = self._extract_count(original_command)
            
            return ParsedCommand(
                action="select_variants",
                tool="select_variants",
                parameters={
                    "num_variants": count,
                    "selection_criteria": "diversity",
                    "session_id": session_id
                },
                session_id=session_id,
                confidence=0.9
            )
        
        elif action == 'select_variants_length':
            count = self._extract_count(original_command)
            
            return ParsedCommand(
                action="select_variants",
                tool="select_variants",
                parameters={
                    "num_variants": count,
                    "selection_criteria": "length",
                    "session_id": session_id
                },
                session_id=session_id,
                confidence=0.9
            )
        
        elif action == 'mutate_sequence':
            count = self._extract_count(original_command)
            sequence = self._extract_sequence(original_command)
            
            return ParsedCommand(
                action="mutate_sequence",
                tool="mutate_sequence",
                parameters={
                    "sequence": sequence,
                    "num_variants": count,
                    "mutation_rate": 0.1  # default
                },
                session_id=session_id,
                confidence=0.8
            )
        
        elif action == 'dna_vendor_research':
            # Extract sequence length if mentioned
            sequence_length = None
            if '96' in original_command or 'variants' in original_command:
                sequence_length = 1000  # Default for variants
            
            return ParsedCommand(
                action="dna_vendor_research",
                tool="dna_vendor_research",
                parameters={
                    "command": original_command,
                    "sequence_length": sequence_length,
                    "quantity": "large"  # Default for vendor research
                },
                session_id=session_id,
                confidence=0.9
            )
        
        elif action == 'align_sequences':
            sequences = self._extract_sequences(original_command)
            
            return ParsedCommand(
                action="align_sequences",
                tool="sequence_alignment",
                parameters={
                    "sequences": sequences,
                    "algorithm": "clustal"  # default
                },
                session_id=session_id,
                confidence=0.8
            )
        
        elif action == 'phylogenetic_tree':
            return ParsedCommand(
                action="phylogenetic_tree",
                tool="phylogenetic_tree",
                parameters={
                    "aligned_sequences": original_command,
                },
                session_id=session_id,
                confidence=0.8
            )
        
        elif action == 'sequence_selection':
            count = self._extract_count(original_command)
            selection_type = "random"  # default to random selection
            
            return ParsedCommand(
                action="sequence_selection",
                tool="sequence_selection",
                parameters={
                    "aligned_sequences": original_command,
                    "selection_type": selection_type,
                    "num_sequences": count,
                },
                session_id=session_id,
                confidence=0.8
            )
        
        elif action == 'synthesis_submission':
            return ParsedCommand(
                action="synthesis_submission",
                tool="synthesis_submission",
                parameters={
                    "sequences": original_command,
                    "vendor_preference": None,
                    "quantity": "standard",
                    "delivery_time": "standard",
                },
                session_id=session_id,
                confidence=0.8
            )
        
        elif action == 'plasmid_visualization':
            # Extract plasmid parameters from command
            # Check if this is a full plasmid (no "into" or "insert" keywords) or an insert
            is_full_plasmid = not re.search(r'\b(?:into|insert|in)\s+', original_command, re.IGNORECASE)
            
            # Extract sequence
            sequence = self._extract_sequence_from_command(original_command)
            
            # Extract vector name if specified
            vector_name = self._extract_vector_name(original_command)
            
            # Extract position if specified
            position_match = re.search(r'\bat\s+position\s+(\d+)', original_command, re.IGNORECASE)
            insert_position = int(position_match.group(1)) if position_match else None
            
            if is_full_plasmid and sequence:
                # Full plasmid sequence provided
                return ParsedCommand(
                    action="plasmid_visualization",
                    tool="plasmid_visualization",
                    parameters={
                        "full_plasmid_sequence": sequence,
                        "vector_name": None,
                        "cloning_sites": "",
                        "insert_sequence": ""
                    },
                    session_id=session_id,
                    confidence=0.8
                )
            else:
                # Insert mode
                cloning_sites = self._extract_cloning_sites(original_command)
                return ParsedCommand(
                    action="plasmid_visualization",
                    tool="plasmid_visualization",
                    parameters={
                        "vector_name": vector_name,
                        "cloning_sites": cloning_sites,
                        "insert_sequence": sequence,
                        "insert_position": insert_position
                    },
                    session_id=session_id,
                    confidence=0.8
                )
        
        else:
            # Default fallback
            return ParsedCommand(
                action="unknown",
                tool="unknown",
                parameters={},
                session_id=session_id,
                confidence=0.1
            )
    
    def _infer_command(self, command: str, session_id: Optional[str] = None) -> ParsedCommand:
        """Infer command type from context when no specific pattern matches."""
        command_lower = command.lower()
        
        # Check for vendor research keywords
        if any(word in command_lower for word in ['order', 'vendor', 'synthesis', 'test', 'assay', 'expression', 'function', 'binding']):
            sequence_length = None
            if '96' in command or 'variants' in command:
                sequence_length = 1000
            
            return ParsedCommand(
                action="dna_vendor_research",
                tool="dna_vendor_research",
                parameters={
                    "command": command,
                    "sequence_length": sequence_length,
                    "quantity": "large"
                },
                session_id=session_id,
                confidence=0.8
            )
        
        # Check for selection-related keywords
        elif any(word in command_lower for word in ['pick', 'select', 'choose']):
            count = self._extract_count(command)
            criteria = self._extract_selection_criteria(command_lower)
            
            return ParsedCommand(
                action="select_variants",
                tool="select_variants",
                parameters={
                    "num_variants": count,
                    "selection_criteria": criteria,
                    "session_id": session_id
                },
                session_id=session_id,
                confidence=0.7
            )
        
        # Check for mutation-related keywords
        elif any(word in command_lower for word in ['mutate', 'generate', 'create', 'variant']):
            count = self._extract_count(command)
            sequence = self._extract_sequence(command)
            
            return ParsedCommand(
                action="mutate_sequence",
                tool="mutate_sequence",
                parameters={
                    "sequence": sequence,
                    "num_variants": count,
                    "mutation_rate": 0.1
                },
                session_id=session_id,
                confidence=0.7
            )
        
        # Default fallback
        return ParsedCommand(
            action="unknown",
            tool="unknown",
            parameters={},
            session_id=session_id,
            confidence=0.1
        )
    
    def _extract_count(self, command: str) -> int:
        """Extract count/number from command."""
        # Look for numbers followed by words like 'sequences', 'variants', etc.
        patterns = [
            r'(\d+)\s+sequences?',
            r'(\d+)\s+variants?',
            r'(\d+)\s+items?',
            r'pick\s+(\d+)',
            r'select\s+(\d+)',
            r'choose\s+(\d+)',
            r'generate\s+(\d+)',
            r'create\s+(\d+)'
        ]
        
        for pattern in patterns:
            match = re.search(pattern, command, re.IGNORECASE)
            if match:
                try:
                    return int(match.group(1))
                except (ValueError, IndexError):
                    continue
        
        # Default count
        return 10
    
    def _extract_selection_criteria(self, command: str) -> str:
        """Extract selection criteria from command."""
        if any(word in command for word in ['randomly', 'random']):
            return 'random'
        elif any(word in command for word in ['diverse', 'diversity']):
            return 'diversity'
        elif any(word in command for word in ['length', 'by length']):
            return 'length'
        else:
            return 'random'  # default
    
    def _extract_sequence(self, command: str) -> str:
        """Extract DNA/RNA sequence from command."""
        # Look for DNA/RNA sequences (A, T, C, G, U)
        sequence_pattern = r'[ATCGU]{4,}'
        match = re.search(sequence_pattern, command.upper())
        if match:
            return match.group(0)
        
        return ""  # Return empty if no sequence found
    
    def _extract_sequences(self, command: str) -> str:
        """Extract multiple sequences from command."""
        # Try to extract sequences after a colon, separated by 'and', ',' or whitespace
        import re
        command_upper = command.upper()
        if ':' in command_upper:
            after_colon = command_upper.split(':', 1)[1]
            # Split on 'and', ',' or whitespace
            seqs = re.split(r'\s+AND\s+|,|\s{2,}', after_colon)
            seqs = [s.strip() for s in seqs if re.search(r'[ATCGU]{4,}', s)]
            if len(seqs) >= 2:
                # Convert to FASTA
                fasta = ''
                for i, seq in enumerate(seqs, 1):
                    # Extract only the sequence part
                    match = re.search(r'([ATCGU]{4,})', seq)
                    if match:
                        fasta += f">sequence_{i}\n{match.group(1)}\n"
                return fasta
        # Fallback: find all long DNA/RNA substrings
        sequences = re.findall(r'[ATCGU]{4,}', command_upper)
        if sequences:
            fasta = ''
            for i, seq in enumerate(sequences, 1):
                fasta += f">sequence_{i}\n{seq}\n"
            return fasta
        return ""
    
    def _extract_cloning_sites(self, command: str) -> str:
        """Extract cloning sites from command."""
        # Look for patterns like "BsaI:123-456, EcoRI:789-1012"
        sites_pattern = r'([A-Za-z]+:\d+-\d+(?:\s*,\s*[A-Za-z]+:\d+-\d+)*)'
        match = re.search(sites_pattern, command)
        if match:
            return match.group(1)
        return "BsaI:100-200, EcoRI:300-400"  # Default sites
    
    def _extract_sequence_from_command(self, command: str) -> str:
        """Extract sequence from plasmid visualization command."""
        # Look for sequence patterns
        sequence_patterns = [
            r'sequence\s+([ATCG]+)',
            r'gene\s+([ATCG]+)',
            r'([ATCG]{10,})',  # Any DNA sequence of 10+ bases
            r'insert\s+([ATCG]+)',
            r'express\s+([ATCG]+)',
        ]
        
        for pattern in sequence_patterns:
            match = re.search(pattern, command, re.IGNORECASE)
            if match:
                return match.group(1).upper()
        
        # If no sequence found, return a default GFP sequence
        return "ATGGTGCACCTGACTGATGCTGAGAAGTCTGCGGTACTGCCTGCTGGGGGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
    
    def _extract_vector_name(self, command: str) -> str:
        """Extract vector name from plasmid visualization command."""
        # Look for vector patterns
        vector_patterns = [
            r'into\s+(\w+)\s+vector',
            r'in\s+(\w+)\s+vector',
            r'vector\s+(\w+)',
            r'(\w+)\s+vector',
        ]
        
        for pattern in vector_patterns:
            match = re.search(pattern, command, re.IGNORECASE)
            if match:
                return match.group(1)
        
        return "pUC19"  # Default vector

def parse_command_raw(command: str, session_id: Optional[str] = None) -> Dict[str, Any]:
    """Raw function for command parsing (for direct calls)."""
    parser = CommandParser()
    parsed = parser.parse_command(command, session_id)
    
    return {
        "action": parsed.action,
        "tool": parsed.tool,
        "parameters": parsed.parameters,
        "session_id": parsed.session_id,
        "confidence": parsed.confidence,
        "original_command": command
    } 