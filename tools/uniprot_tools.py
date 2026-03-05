"""
UniProt database access tools for protein sequences and metadata.
"""

import os
import requests
from typing import Dict, Any, Optional, List
import json

UNIPROT_API_BASE = "https://rest.uniprot.org"


def query_uniprot(
    query: str,
    format: str = "fasta",
    limit: int = 10
) -> Dict[str, Any]:
    """
    Query UniProt protein database.
    
    Args:
        query: Search query (e.g., "BRCA1 AND organism:9606", "P38398")
        format: Return format ('fasta', 'json', 'xml', 'tsv')
        limit: Maximum number of results
    
    Returns:
        Dictionary containing search results
    """
    # Short-circuit in mock mode to avoid network calls
    if os.getenv("HELIX_MOCK_MODE"):
        return {
            "status": "success",
            "query": query,
            "results": [
                {
                    "header": f"Mock UniProt entry for {query}",
                    "sequence": "M" * 50
                }
            ],
            "count": 1,
            "format": format,
            "warning": "HELIX_MOCK_MODE enabled; returning mock UniProt results"
        }
    
    try:
        # Search UniProt
        search_url = f"{UNIPROT_API_BASE}/uniprotkb/search"
        params = {
            "query": query,
            "format": format,
            "size": limit
        }
        
        response = requests.get(search_url, params=params, timeout=30)
        response.raise_for_status()
        
        if format == "fasta":
            # Parse FASTA results
            results = []
            current_entry = None
            
            for line in response.text.split('\n'):
                if line.startswith('>'):
                    if current_entry:
                        results.append(current_entry)
                    current_entry = {
                        "header": line[1:].strip(),
                        "sequence": ""
                    }
                elif current_entry and line.strip():
                    current_entry["sequence"] += line.strip()
            
            if current_entry:
                results.append(current_entry)
            
            return {
                "status": "success",
                "query": query,
                "results": results,
                "count": len(results),
                "format": format
            }
        else:
            # For JSON/XML/TSV, return raw content
            return {
                "status": "success",
                "query": query,
                "content": response.text,
                "format": format
            }
            
    except Exception as e:
        # Offline/network fallback stub to keep UX flowing
        return {
            "status": "success",
            "query": query,
            "results": [],
            "count": 0,
            "format": format,
            "warning": f"Network/unavailable; returning empty results. Details: {e}"
        }


def get_uniprot_entry(accession: str, format: str = "fasta") -> Dict[str, Any]:
    """
    Get specific UniProt entry by accession.
    
    Args:
        accession: UniProt accession (e.g., "P38398", "P04637")
        format: Return format ('fasta', 'json', 'xml')
    
    Returns:
        Dictionary containing entry data
    """
    # Short-circuit in mock mode to avoid network calls
    if os.getenv("HELIX_MOCK_MODE"):
        mock_seq = "M" * 100
        return {
            "status": "success",
            "accession": accession,
            "header": f"{accession} mock entry",
            "sequence": mock_seq,
            "length": len(mock_seq),
            "warning": "HELIX_MOCK_MODE enabled; returning mock UniProt entry"
        }
    
    try:
        url = f"{UNIPROT_API_BASE}/uniprotkb/{accession}.{format}"
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        
        if format == "fasta":
            lines = response.text.split('\n')
            header = lines[0][1:] if lines[0].startswith('>') else ""
            sequence = ''.join(lines[1:]).replace('\n', '')
            
            return {
                "status": "success",
                "accession": accession,
                "header": header,
                "sequence": sequence,
                "length": len(sequence)
            }
        else:
            return {
                "status": "success",
                "accession": accession,
                "content": response.text,
                "format": format
            }
            
    except Exception as e:
        # Offline/network fallback stub
        return {
            "status": "success",
            "accession": accession,
            "header": f"{accession} (offline stub)",
            "sequence": "",
            "length": 0,
            "warning": f"Network/unavailable; returning stub. Details: {e}"
        }
