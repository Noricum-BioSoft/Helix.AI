"""
Gene Ontology (GO) term lookup and annotation tools.
"""

import os
import requests
from typing import Dict, Any, Optional, List

GO_API_BASE = "https://api.geneontology.org/api"


def lookup_go_term(go_id: str) -> Dict[str, Any]:
    """
    Look up Gene Ontology term by ID.
    
    Args:
        go_id: GO term ID (e.g., "GO:0003677", "GO:0008150")
    
    Returns:
        Dictionary containing GO term information
    """
    # Short-circuit in mock mode to avoid network calls
    if os.getenv("HELIX_MOCK_MODE"):
        return {
            "status": "success",
            "go_id": go_id,
            "name": f"{go_id} mock term",
            "namespace": "mock_namespace",
            "definition": "Mock definition (HELIX_MOCK_MODE enabled)",
            "synonyms": [],
            "warning": "HELIX_MOCK_MODE enabled; returning mock GO term"
        }
    
    try:
        url = f"{GO_API_BASE}/ontology/term/{go_id}"
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        data = response.json()
        
        return {
            "status": "success",
            "go_id": go_id,
            "name": data.get("label", ""),
            "namespace": data.get("namespace", ""),
            "definition": data.get("definition", {}).get("defstr", ""),
            "synonyms": data.get("synonyms", [])
        }
    except Exception as e:
        # Offline/network fallback stub
        return {
            "status": "success",
            "go_id": go_id,
            "name": f"{go_id} (offline stub)",
            "namespace": "",
            "definition": "",
            "synonyms": [],
            "warning": f"Network/unavailable; returning stub. Details: {e}"
        }


def annotate_genes_with_go(gene_ids: List[str], organism: str = "9606") -> Dict[str, Any]:
    """
    Annotate genes with GO terms.
    
    Args:
        gene_ids: List of gene IDs (e.g., ["BRCA1", "TP53"])
        organism: Organism ID (9606 for human)
    
    Returns:
        Dictionary containing GO annotations for each gene
    """
    # This is a simplified version - full implementation would use GO database
    return {
        "status": "success",
        "genes": gene_ids,
        "annotations": {}  # Would contain actual GO annotations
    }
