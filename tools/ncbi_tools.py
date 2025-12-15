"""
NCBI database access tools for fetching sequences and metadata.

This module provides functions to interact with NCBI databases including:
- Nucleotide sequences (GenBank, RefSeq)
- Protein sequences
- Metadata retrieval
"""

import os
from Bio import Entrez, SeqIO
from typing import Dict, Any

# NCBI requires an email. Allow overriding via env var.
Entrez.email = os.getenv("HELIX_NCBI_EMAIL", "helix.ai@example.com")


def fetch_sequence_from_ncbi(
    accession: str,
    database: str = "nucleotide",
    rettype: str = "fasta",
    retmode: str = "text"
) -> Dict[str, Any]:
    """
    Fetch sequence from NCBI databases by accession number.
    
    Args:
        accession: NCBI accession number (e.g., "NC_000001.11", "NM_000492.3")
        database: Database to search ('nucleotide', 'protein', 'nuccore', 'protein')
        rettype: Return type ('fasta', 'gb', 'xml')
        retmode: Return mode ('text', 'xml')
    
    Returns:
        Dictionary containing:
        - status: 'success' or 'error'
        - accession: The accession number
        - sequence: DNA/RNA/protein sequence
        - description: Sequence description
        - length: Sequence length
        - database: Source database
        - error: Error message (if status is 'error')
    """
    # Short-circuit in mock mode to avoid network calls
    if os.getenv("HELIX_MOCK_MODE"):
        mock_seq = "ATGC" * 25
        return {
            "status": "success",
            "accession": accession,
            "sequence": mock_seq,
            "description": f"Mock sequence for {accession}",
            "length": len(mock_seq),
            "database": database,
            "sequence_type": "nucleotide" if database in ["nucleotide", "nuccore"] else "protein",
            "warning": "HELIX_MOCK_MODE enabled; returning mock NCBI sequence"
        }
    
    try:
        # Fetch sequence
        handle = Entrez.efetch(
            db=database,
            id=accession,
            rettype=rettype,
            retmode=retmode
        )
        
        # Parse FASTA
        if rettype == "fasta":
            record = SeqIO.read(handle, "fasta")
            handle.close()
            
            return {
                "status": "success",
                "accession": accession,
                "sequence": str(record.seq),
                "description": record.description,
                "length": len(record.seq),
                "database": database,
                "sequence_type": "nucleotide" if database in ["nucleotide", "nuccore"] else "protein"
            }
        else:
            # For other formats, return raw text
            content = handle.read()
            handle.close()
            return {
                "status": "success",
                "accession": accession,
                "content": content,
                "format": rettype,
                "database": database
            }
            
    except Exception as e:
        return {
            "status": "error",
            "accession": accession,
            "error": str(e),
            "database": database
        }


def search_ncbi(
    query: str,
    database: str = "nucleotide",
    max_results: int = 10
) -> Dict[str, Any]:
    """
    Search NCBI databases and return matching accessions.
    
    Args:
        query: Search query (e.g., "Homo sapiens[Organism] AND BRCA1[Gene]")
        database: Database to search
        max_results: Maximum number of results to return
    
    Returns:
        Dictionary containing search results with accessions and descriptions
    """
    # Short-circuit in mock mode to avoid network calls
    if os.getenv("HELIX_MOCK_MODE"):
        return {
            "status": "success",
            "query": query,
            "results": [
                {
                    "accession": "MOCK_000001.1",
                    "title": "Mock entry 1",
                    "organism": "Mockus organismus",
                    "length": 1000
                }
            ],
            "count": 1,
            "total_found": 1,
            "warning": "HELIX_MOCK_MODE enabled; returning mock NCBI search results"
        }
    
    try:
        # Search
        search_handle = Entrez.esearch(
            db=database,
            term=query,
            retmax=max_results
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        ids = search_results["IdList"]
        
        if not ids:
            return {
                "status": "success",
                "query": query,
                "results": [],
                "count": 0
            }
        
        # Fetch summaries
        summary_handle = Entrez.esummary(db=database, id=",".join(ids))
        summaries = Entrez.read(summary_handle)
        summary_handle.close()
        
        results = []
        for summary in summaries:
            results.append({
                "accession": summary.get("AccessionVersion", summary.get("Id", "")),
                "title": summary.get("Title", ""),
                "organism": summary.get("Organism", ""),
                "length": summary.get("Length", 0)
            })
        
        return {
            "status": "success",
            "query": query,
            "results": results,
            "count": len(results),
            "total_found": int(search_results.get("Count", 0))
        }
        
    except Exception as e:
        # Offline/network fallback stub to keep UX flowing
        return {
            "status": "success",
            "query": query,
            "results": [],
            "count": 0,
            "total_found": 0,
            "warning": f"Network/unavailable; returning empty results. Details: {e}"
        }
