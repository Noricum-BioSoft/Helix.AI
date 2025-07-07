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
                r'pick\s+(\d+)\s+sequences?\s+randomly',
                r'select\s+(\d+)\s+sequences?\s+randomly',
                r'choose\s+(\d+)\s+sequences?\s+randomly',
                r'from\s+.*variants?.*pick\s+(\d+)\s+sequences?\s+randomly',
                r'from\s+.*variants?.*select\s+(\d+)\s+sequences?\s+randomly',
                r'randomly\s+pick\s+(\d+)\s+sequences?',
                r'randomly\s+select\s+(\d+)\s+sequences?',
            ],
            'select_variants_diverse': [
                r'pick\s+(\d+)\s+diverse\s+sequences?',
                r'select\s+(\d+)\s+diverse\s+sequences?',
                r'choose\s+(\d+)\s+diverse\s+sequences?',
                r'from\s+.*variants?.*pick\s+(\d+)\s+diverse\s+sequences?',
                r'from\s+.*variants?.*select\s+(\d+)\s+diverse\s+sequences?',
            ],
            'select_variants_length': [
                r'pick\s+(\d+)\s+sequences?\s+by\s+length',
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
            'align_sequences': [
                r'align\s+sequences?',
                r'perform\s+alignment',
                r'sequence\s+alignment',
            ],
            'analyze_data': [
                r'analyze\s+data',
                r'perform\s+analysis',
                r'data\s+analysis',
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
        
        # Check for selection-related keywords
        if any(word in command_lower for word in ['pick', 'select', 'choose']):
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
        # Look for multiple DNA/RNA sequences
        sequences = re.findall(r'[ATCGU]{4,}', command.upper())
        if sequences:
            # Convert to FASTA format
            fasta = ""
            for i, seq in enumerate(sequences, 1):
                fasta += f">sequence_{i}\n{seq}\n"
            return fasta
        
        return ""

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