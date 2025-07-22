import re
from typing import Dict, Any, List, Optional
from dataclasses import dataclass

@dataclass
class PlasmidFeature:
    """Represents a feature on a plasmid."""
    name: str
    start: int
    end: int
    type: str
    color: str = "#ff0000"
    description: str = ""

@dataclass
class PlasmidData:
    """Represents complete plasmid data for visualization."""
    name: str
    sequence: str
    features: List[PlasmidFeature]
    size: int
    description: str = ""

def parse_cloning_sites(sites_str: str) -> List[Dict[str, Any]]:
    """Parse cloning sites string into structured data."""
    sites = []
    # Expected format: "BsaI:123-456, EcoRI:789-1012"
    site_pattern = r'(\w+):(\d+)-(\d+)'
    matches = re.findall(site_pattern, sites_str)
    
    for match in matches:
        enzyme, start, end = match
        sites.append({
            "name": enzyme,
            "start": int(start),
            "end": int(end),
            "type": "restriction_site",
            "color": "#00ff00"
        })
    
    return sites

def create_plasmid_data(vector_name: str, cloning_sites: str, insert_sequence: str) -> PlasmidData:
    """Create plasmid data from user inputs."""
    
    # Clean inputs
    vector_name = vector_name.strip()
    cloning_sites = cloning_sites.strip()
    insert_sequence = insert_sequence.strip().upper()
    
    # Validate insert sequence
    valid_bases = set("ATCG")
    if not all(base in valid_bases for base in insert_sequence):
        raise ValueError(f"Invalid insert sequence '{insert_sequence}'. Only A, T, C, G allowed.")
    
    # Parse cloning sites
    sites = parse_cloning_sites(cloning_sites)
    
    # Create features list
    features = []
    
    # Add vector backbone features
    features.extend([
        PlasmidFeature(
            name="Origin of Replication",
            start=0,
            end=432,
            type="rep_origin",
            color="#feca57",
            description="Origin of replication"
        ),
        PlasmidFeature(
            name="Ampicillin Resistance",
            start=433,
            end=1200,
            type="CDS",
            color="#ff6b6b",
            description="Ampicillin resistance gene"
        ),
        PlasmidFeature(
            name="Multiple Cloning Site",
            start=1201,
            end=1400,
            type="misc_feature",
            color="#95a5a6",
            description="Multiple cloning site"
        )
    ])
    
    # Add cloning sites as features
    for site in sites:
        features.append(PlasmidFeature(
            name=f"{site['name']} Site",
            start=site["start"],
            end=site["end"],
            type="restriction_site",
            color="#00ff00",
            description=f"Restriction site for {site['name']}"
        ))
    
    # Add insert sequence as a feature
    if insert_sequence:
        # Place insert after the last cloning site or at position 1401
        insert_start = max([site["end"] for site in sites]) if sites else 1401
        insert_end = insert_start + len(insert_sequence) - 1
        
        features.append(PlasmidFeature(
            name="Inserted Sequence",
            start=insert_start,
            end=insert_end,
            type="insert",
            color="#0000ff",
            description=f"Inserted sequence: {insert_sequence[:50]}{'...' if len(insert_sequence) > 50 else ''}"
        ))
    
    # Create complete sequence (synthetic vector + insert)
    vector_backbone = "ATGGTGCACCTGACTGATGCTGAGAAGTCTGCGGTACTGCCTGCTGGGGGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAA"
    
    # Calculate total size
    total_size = max([f.end for f in features]) if features else len(vector_backbone) + len(insert_sequence)
    
    # Create complete sequence by repeating vector backbone and adding insert
    complete_sequence = vector_backbone * (total_size // len(vector_backbone) + 1)
    complete_sequence = complete_sequence[:total_size]
    
    return PlasmidData(
        name=vector_name,
        sequence=complete_sequence,
        features=features,
        size=total_size,
        description=f"Plasmid {vector_name} with {len(sites)} cloning sites and {len(insert_sequence)}bp insert"
    )

def run_plasmid_visualization_raw(vector_name: str, cloning_sites: str, insert_sequence: str) -> Dict[str, Any]:
    """Generate plasmid visualization data from user inputs."""
    
    try:
        # Create plasmid data
        plasmid_data = create_plasmid_data(vector_name, cloning_sites, insert_sequence)
        
        # Convert to dictionary format for API response
        features_dict = []
        for feature in plasmid_data.features:
            features_dict.append({
                "name": feature.name,
                "start": feature.start,
                "end": feature.end,
                "type": feature.type,
                "color": feature.color,
                "description": feature.description
            })
        
        result = {
            "text": f"Plasmid visualization data generated for {plasmid_data.name}",
            "plasmid_data": {
                "name": plasmid_data.name,
                "sequence": plasmid_data.sequence,
                "features": features_dict,
                "size": plasmid_data.size,
                "description": plasmid_data.description
            },
            "visualization_type": "circular_plasmid",
            "metadata": {
                "vector_name": vector_name,
                "cloning_sites": cloning_sites,
                "insert_sequence": insert_sequence,
                "feature_count": len(plasmid_data.features)
            }
        }
        
        return result
        
    except ValueError as e:
        return {
            "text": f"Error: {str(e)}",
            "error": str(e)
        }
    except Exception as e:
        return {
            "text": f"Unexpected error: {str(e)}",
            "error": str(e)
        }

from langchain.agents import tool

@tool
def run_plasmid_visualization(vector_name: str, cloning_sites: str, insert_sequence: str):
    """Generate plasmid visualization data from vector name, cloning sites, and insert sequence."""
    return run_plasmid_visualization_raw(vector_name, cloning_sites, insert_sequence)

def create_plasmid_for_representatives(representatives: List[str], aligned_sequences: str, vector_name: str = "pUC19", cloning_sites: str = "EcoRI, BamHI, HindIII") -> Dict[str, Any]:
    """Create plasmid visualizations for representative sequences from clustering."""
    
    try:
        # Parse aligned sequences to get the actual sequences
        sequences = parse_aligned_sequences(aligned_sequences)
        sequence_dict = {seq["name"]: seq["sequence"] for seq in sequences}
        
        # Create plasmid data for each representative
        plasmid_results = []
        
        for i, rep_name in enumerate(representatives):
            if rep_name in sequence_dict:
                insert_sequence = sequence_dict[rep_name]
                
                # Create plasmid data for this representative
                plasmid_data = create_plasmid_data(
                    vector_name=f"{vector_name}_{rep_name}",
                    cloning_sites=cloning_sites,
                    insert_sequence=insert_sequence
                )
                
                plasmid_results.append({
                    "representative_name": rep_name,
                    "plasmid_data": plasmid_data,
                    "insert_sequence": insert_sequence,
                    "sequence_length": len(insert_sequence)
                })
        
        # Create summary text
        summary_text = f"Created plasmid visualizations for {len(plasmid_results)} representative sequences:\n\n"
        for i, result in enumerate(plasmid_results, 1):
            summary_text += f"{i}. {result['representative_name']} ({result['sequence_length']} bp)\n"
        
        return {
            "text": summary_text,
            "plasmid_results": plasmid_results,
            "total_representatives": len(plasmid_results),
            "vector_name": vector_name,
            "cloning_sites": cloning_sites
        }
        
    except Exception as e:
        return {
            "text": f"Error creating plasmid visualizations: {str(e)}",
            "error": str(e)
        }

def parse_aligned_sequences(alignment_data: str) -> List[Dict[str, str]]:
    """Parse aligned sequences from string format."""
    sequences = []
    current_name = ""
    current_sequence = ""
    
    # Check if this is colon-separated format (name: sequence)
    if ':' in alignment_data and not alignment_data.startswith('>'):
        # Parse colon-separated format
        for line in alignment_data.strip().split('\n'):
            line = line.strip()
            if ':' in line:
                parts = line.split(':', 1)
                if len(parts) == 2:
                    name = parts[0].strip()
                    sequence = parts[1].strip()
                    sequences.append({
                        "name": name,
                        "sequence": sequence
                    })
        return sequences
    
    # Parse FASTA format (with > headers)
    for line in alignment_data.strip().split('\n'):
        line = line.strip()
        if line.startswith('>'):
            # Save previous sequence if exists
            if current_name and current_sequence:
                sequences.append({
                    "name": current_name,
                    "sequence": current_sequence
                })
            # Start new sequence - check if sequence is on same line
            parts = line[1:].split(' ', 1)  # Split on space, max 1 split
            if len(parts) == 2:
                # Sequence is on same line as header
                current_name = parts[0].strip()
                current_sequence = parts[1].strip()
            else:
                # Sequence is on separate line
                current_name = line[1:].strip()
                current_sequence = ""
        else:
            # Add to current sequence
            current_sequence += line
    
    # Add the last sequence
    if current_name and current_sequence:
        sequences.append({
            "name": current_name,
            "sequence": current_sequence
        })
    
    return sequences 

@tool
def run_plasmid_for_representatives(representatives: List[str], aligned_sequences: str, vector_name: str = "pUC19", cloning_sites: str = "EcoRI, BamHI, HindIII"):
    """Create plasmid visualizations for representative sequences from clustering analysis."""
    return create_plasmid_for_representatives(representatives, aligned_sequences, vector_name, cloning_sites) 