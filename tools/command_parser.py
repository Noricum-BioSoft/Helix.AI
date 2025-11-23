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
            ],
            'read_trimming': [
                r'trim\s+.*reads?',
                r'trim\s+.*fastq',
                r'quality\s+trim',
                r'remove\s+.*quality\s+.*bases?',
                r'trim\s+.*low.*quality',
                r'quality\s+.*trimming',
                r'remove\s+adapter',
                r'remove\s+adapter\s+sequences?',
                r'adapter\s+removal',
            ],
            'read_merging': [
                r'merge\s+.*reads?',
                r'merge\s+.*paired.*end',
                r'combine\s+.*reads?',
                r'merge\s+.*r1.*r2',
                r'merge\s+.*forward.*reverse',
            ],
            'quality_assessment': [
                r'quality\s+report',  # "quality report"
                r'quality\s+\w+\s+report',  # "quality a report" or "quality assessment report"
                r'quality\s+assessment',
                r'generate\s+\w+\s+quality\s+report',  # "generate a quality report"
                r'generate\s+quality',
                r'assess\s+quality',
                r'quality\s+check',
            ],
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
        print(f"🔧 [PARSER] Parsing command: '{command[:100]}...'")  # Log first 100 chars
        
        # Try to match action patterns
        for action, patterns in self.action_patterns.items():
            for pattern in patterns:
                match = re.search(pattern, command_lower, re.IGNORECASE)
                if match:
                    print(f"🔧 [PARSER] Matched action: {action} with pattern: {pattern}")
                    return self._build_parsed_command(action, command, match, session_id)
        
        print("🔧 [PARSER] No specific action pattern matched, trying inference")
        # If no specific action found, try to infer from context
        inferred = self._infer_command(command, session_id)
        print(f"🔧 [PARSER] Inferred action: {inferred.action}, tool: {inferred.tool}, confidence: {inferred.confidence}")
        return inferred
    
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
            vector_name = self._extract_vector_name(original_command)
            cloning_sites = self._extract_cloning_sites(original_command)
            insert_sequence = self._extract_sequence_from_command(original_command)
            
            return ParsedCommand(
                action="plasmid_visualization",
                tool="plasmid_visualization",
                parameters={
                    "vector_name": vector_name,
                    "cloning_sites": cloning_sites,
                    "insert_sequence": insert_sequence
                },
                session_id=session_id,
                confidence=0.8
            )
        
        elif action == 'read_trimming':
            # Extract FASTQ data and parameters - support multiple files
            quality_threshold = self._extract_quality_threshold(original_command)
            adapter = self._extract_adapter(original_command)
            print(f"🔧 [PARSER] Extracted adapter: {adapter}")
            
            # Try to extract paired-end reads first
            forward_reads, reverse_reads = self._extract_paired_reads(original_command)
            
            print(f"🔧 [PARSER] After extraction - forward: {len(forward_reads) if forward_reads else 0} chars, reverse: {len(reverse_reads) if reverse_reads else 0} chars")
            
            # If we have both forward and reverse, use them
            if forward_reads and reverse_reads and len(forward_reads) > 10 and len(reverse_reads) > 10:
                print(f"🔧 [PARSER] Using paired-end mode")
                return ParsedCommand(
                    action="read_trimming",
                    tool="read_trimming",
                    parameters={
                        "forward_reads": forward_reads,
                        "reverse_reads": reverse_reads,
                        "quality_threshold": quality_threshold,
                        "adapter": adapter,
                        "file_type": "paired_end"
                    },
                    session_id=session_id,
                    confidence=0.9
                )
            else:
                # Fall back to single file extraction
                print(f"🔧 [PARSER] Falling back to single file mode")
                reads = self._extract_fastq_data(original_command)
                return ParsedCommand(
                    action="read_trimming",
                    tool="read_trimming",
                    parameters={
                        "reads": reads,
                        "quality_threshold": quality_threshold,
                        "adapter": adapter,
                        "file_type": "single"
                    },
                    session_id=session_id,
                    confidence=0.9
                )
        
        elif action == 'read_merging':
            # Extract forward and reverse reads
            forward_reads, reverse_reads = self._extract_paired_reads(original_command)
            min_overlap = self._extract_min_overlap(original_command)
            
            print(f"🔧 [PARSER] After extraction for merging - forward: {len(forward_reads) if forward_reads else 0} chars, reverse: {len(reverse_reads) if reverse_reads else 0} chars")
            
            if not forward_reads or not reverse_reads:
                print(f"🔧 [PARSER] WARNING: Missing paired reads for merging!")
            
            return ParsedCommand(
                action="read_merging",
                tool="read_merging",
                parameters={
                    "forward_reads": forward_reads,
                    "reverse_reads": reverse_reads,
                    "min_overlap": min_overlap
                },
                session_id=session_id,
                confidence=0.9
            )
        
        elif action == 'quality_assessment':
            # Quality assessment - will retrieve merged sequences from history
            return ParsedCommand(
                action="quality_assessment",
                tool="quality_assessment",
                parameters={},
                session_id=session_id,
                confidence=0.9
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
        
        # Check for read trimming keywords
        elif any(word in command_lower for word in ['trim', 'quality', 'fastq', 'reads']):
            reads = self._extract_fastq_data(command)
            quality_threshold = self._extract_quality_threshold(command)
            adapter = self._extract_adapter(command)
            
            return ParsedCommand(
                action="read_trimming",
                tool="read_trimming",
                parameters={
                    "reads": reads,
                    "quality_threshold": quality_threshold,
                    "adapter": adapter
                },
                session_id=session_id,
                confidence=0.8
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
    
    def _extract_fastq_data(self, command: str) -> str:
        """Extract FASTQ data from command (may include 'File content:', 'Forward reads:', 'Reverse reads:', or 'File N:' prefix)."""
        print(f"🔧 [PARSER] Extracting FASTQ data from command (length: {len(command)})")
        
        # Look for "Forward reads:" or "Reverse reads:" (for paired-end)
        if "Forward reads:" in command or "forward reads:" in command.lower():
            parts = command.split("Forward reads:", 1) if "Forward reads:" in command else command.split("forward reads:", 1)
            if len(parts) > 1:
                # Extract until "Reverse reads:" or end
                if "Reverse reads:" in parts[1] or "reverse reads:" in parts[1].lower():
                    extracted = parts[1].split("Reverse reads:", 1)[0].split("reverse reads:", 1)[0].strip()
                else:
                    extracted = parts[1].strip()
                print(f"🔧 [PARSER] Extracted forward FASTQ data (length: {len(extracted)})")
                return extracted
        
        # Look for "Reverse reads:"
        if "Reverse reads:" in command or "reverse reads:" in command.lower():
            parts = command.split("Reverse reads:", 1) if "Reverse reads:" in command else command.split("reverse reads:", 1)
            if len(parts) > 1:
                extracted = parts[1].strip()
                print(f"🔧 [PARSER] Extracted reverse FASTQ data (length: {len(extracted)})")
                return extracted
        
        # Look for "File content:" followed by FASTQ data (legacy single file)
        if "File content:" in command:
            parts = command.split("File content:", 1)
            if len(parts) > 1:
                extracted = parts[1].strip()
                print(f"🔧 [PARSER] Extracted FASTQ data from 'File content:' (length: {len(extracted)})")
                return extracted
        
        # Look for "File 1:" or "File 2:" pattern
        if re.search(r'File\s+\d+.*?:', command, re.IGNORECASE):
            # Extract first file
            match = re.search(r'File\s+\d+.*?:\s*\n(.*?)(?=File\s+\d+.*?:|$)', command, re.IGNORECASE | re.DOTALL)
            if match:
                extracted = match.group(1).strip()
                print(f"🔧 [PARSER] Extracted FASTQ data from 'File N:' pattern (length: {len(extracted)})")
                return extracted
        
        # Also check for FASTQ format directly (@ header lines)
        if "@" in command and ("+" in command or len(command) > 100):
            # Try to extract FASTQ block
            lines = command.split('\n')
            fastq_lines = []
            in_fastq = False
            for line in lines:
                if line.startswith('@'):
                    in_fastq = True
                if in_fastq:
                    fastq_lines.append(line)
                    # FASTQ is 4 lines per read, stop after reasonable number
                    if len(fastq_lines) > 40:
                        break
            if fastq_lines:
                extracted = '\n'.join(fastq_lines)
                print(f"🔧 [PARSER] Extracted FASTQ data from direct format (length: {len(extracted)})")
                return extracted
        
        print(f"🔧 [PARSER] No FASTQ pattern found, returning empty string")
        return ""  # Return empty string when no FASTQ data found - executor will retrieve from history  # Return full command if no FASTQ pattern found
    
    def _extract_quality_threshold(self, command: str) -> int:
        """Extract quality threshold from command."""
        # Look for patterns like "quality threshold 20", "threshold of 20", "Phred 20"
        patterns = [
            r'quality\s+threshold\s+(\d+)',
            r'threshold\s+of\s+(\d+)',
            r'phred\s+(\d+)',
            r'threshold\s+(\d+)',
            r'quality\s+(\d+)',
        ]
        
        for pattern in patterns:
            match = re.search(pattern, command, re.IGNORECASE)
            if match:
                try:
                    return int(match.group(1))
                except (ValueError, IndexError):
                    continue
        
        return 20  # Default quality threshold
    
    def _extract_adapter(self, command: str) -> Optional[str]:
        """Extract adapter sequence from command."""
        # Look for patterns like "adapter AGATCGGAAGAGC", "remove adapter sequences AGATCGGAAGAGC", etc.
        adapter_patterns = [
            r'adapter\s+sequences?\s+([ATCGN]{6,})',  # "adapter sequences AGATCGGAAGAGC" or "adapter sequence AGATCGGAAGAGC"
            r'remove\s+adapter\s+sequences?\s+([ATCGN]{6,})',  # "remove adapter sequences AGATCGGAAGAGC"
            r'adapter\s+([ATCGN]{6,})',  # "adapter AGATCGGAAGAGC"
            r'remove\s+adapter\s+([ATCGN]{6,})',  # "remove adapter AGATCGGAAGAGC"
        ]
        
        for pattern in adapter_patterns:
            match = re.search(pattern, command, re.IGNORECASE)
            if match:
                adapter = match.group(1).upper()
                print(f"🔧 [PARSER] Found adapter: {adapter}")
                return adapter
        
        print(f"🔧 [PARSER] No adapter found in command")
        return None
    
    def _extract_paired_reads(self, command: str) -> Tuple[str, str]:
        """Extract forward and reverse reads from command."""
        forward_reads = ""
        reverse_reads = ""
        
        print(f"🔧 [PARSER] Extracting paired reads from command (length: {len(command)})")
        
        # Split command by "Forward reads" and "Reverse reads" markers to find both sections
        # Handle both orders: Forward first or Reverse first
        
        # Find all occurrences of "Forward reads" and "Reverse reads"
        forward_markers = list(re.finditer(r'Forward\s+reads\s*(?:\([^)]+\))?:\s*\n', command, re.IGNORECASE))
        reverse_markers = list(re.finditer(r'Reverse\s+reads\s*(?:\([^)]+\))?:\s*\n', command, re.IGNORECASE))
        
        print(f"🔧 [PARSER] Found {len(forward_markers)} forward markers, {len(reverse_markers)} reverse markers")
        
        # Extract forward reads
        if forward_markers:
            forward_start = forward_markers[0].end()
            # Find where forward reads end (either at reverse reads marker or end of command)
            forward_end = len(command)
            if reverse_markers:
                # Find the reverse marker that comes after this forward marker
                for rev_marker in reverse_markers:
                    if rev_marker.start() > forward_start:
                        forward_end = rev_marker.start()
                        break
            forward_reads = command[forward_start:forward_end].strip()
            print(f"🔧 [PARSER] Extracted forward reads: {len(forward_reads)} chars")
            if len(forward_reads) > 0:
                print(f"🔧 [PARSER] Forward reads preview: {forward_reads[:100]}...")
        
        # Extract reverse reads
        if reverse_markers:
            reverse_start = reverse_markers[0].end()
            # Find where reverse reads end (either at forward reads marker or end of command)
            reverse_end = len(command)
            if forward_markers:
                # Find the forward marker that comes after this reverse marker
                for fwd_marker in forward_markers:
                    if fwd_marker.start() > reverse_start:
                        reverse_end = fwd_marker.start()
                        break
            reverse_reads = command[reverse_start:reverse_end].strip()
            print(f"🔧 [PARSER] Extracted reverse reads: {len(reverse_reads)} chars")
            if len(reverse_reads) > 0:
                print(f"🔧 [PARSER] Reverse reads preview: {reverse_reads[:100]}...")
        
        # Also check for "File 1:" and "File 2:" with R1/R2 in filename (fallback)
        if not forward_reads or not reverse_reads:
            file_pattern = re.compile(r'(?:Forward\s+reads|Reverse\s+reads|File\s+\d+)\s*\(([^)]+)\):\s*\n(.*?)(?=(?:Forward\s+reads|Reverse\s+reads|File\s+\d+)|$)', re.IGNORECASE | re.DOTALL)
            matches = list(file_pattern.finditer(command))
            
            print(f"🔧 [PARSER] Found {len(matches)} file patterns (fallback)")
            for match in matches:
                filename = match.group(1)
                content = match.group(2).strip()
                print(f"🔧 [PARSER] Checking file: {filename} ({len(content)} chars)")
                
                if re.search(r'[._-]R?1[._-]', filename, re.IGNORECASE) or filename.upper().endswith('_R1') or filename.upper().endswith('.R1'):
                    if not forward_reads:
                        forward_reads = content
                        print(f"🔧 [PARSER] Identified as forward (R1): {filename}")
                elif re.search(r'[._-]R?2[._-]', filename, re.IGNORECASE) or filename.upper().endswith('_R2') or filename.upper().endswith('.R2'):
                    if not reverse_reads:
                        reverse_reads = content
                        print(f"🔧 [PARSER] Identified as reverse (R2): {filename}")
        
        # If we still don't have both, try to get them by order (first file = forward, second = reverse)
        if not forward_reads or not reverse_reads:
            # Try to find any file patterns
            file_pattern = re.compile(r'(?:Forward\s+reads|Reverse\s+reads|File\s+\d+)\s*\(([^)]+)\):\s*\n(.*?)(?=(?:Forward\s+reads|Reverse\s+reads|File\s+\d+)|$)', re.IGNORECASE | re.DOTALL)
            matches = list(file_pattern.finditer(command))
            
            if len(matches) >= 2:
                if not forward_reads:
                    forward_reads = matches[0].group(2).strip()
                    print(f"🔧 [PARSER] Using first file as forward: {matches[0].group(1)}")
                if not reverse_reads:
                    reverse_reads = matches[1].group(2).strip()
                    print(f"🔧 [PARSER] Using second file as reverse: {matches[1].group(1)}")
        
        print(f"🔧 [PARSER] Final extraction - forward: {len(forward_reads)} chars, reverse: {len(reverse_reads)} chars")
        return forward_reads, reverse_reads
    
    def _extract_min_overlap(self, command: str) -> int:
        """Extract minimum overlap from command."""
        patterns = [
            r'overlap\s+of\s+(\d+)',
            r'minimum\s+overlap\s+(\d+)',
            r'min\s+overlap\s+(\d+)',
            r'overlap\s+(\d+)',
        ]
        
        for pattern in patterns:
            match = re.search(pattern, command, re.IGNORECASE)
            if match:
                try:
                    return int(match.group(1))
                except (ValueError, IndexError):
                    continue
        
        return 12  # Default minimum overlap

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