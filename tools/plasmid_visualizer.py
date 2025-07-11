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
    
    # Add cloning sites as features
    for site in sites:
        features.append(PlasmidFeature(
            name=site["name"],
            start=site["start"],
            end=site["end"],
            type="restriction_site",
            color=site["color"],
            description=f"Restriction site for {site['name']}"
        ))
    
    # Add insert sequence as a feature
    if insert_sequence:
        # Place insert after the last cloning site or at position 1
        insert_start = max([site["end"] for site in sites]) if sites else 1
        insert_end = insert_start + len(insert_sequence) - 1
        
        features.append(PlasmidFeature(
            name="Insert",
            start=insert_start,
            end=insert_end,
            type="insert",
            color="#0000ff",
            description=f"Inserted sequence: {insert_sequence}"
        ))
    
    # Create complete sequence (placeholder - in real implementation, this would be the actual vector sequence)
    # For now, we'll create a synthetic sequence based on the features
    total_size = max([f.end for f in features]) if features else len(insert_sequence)
    sequence = "A" * total_size  # Placeholder sequence
    
    return PlasmidData(
        name=vector_name,
        sequence=sequence,
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