import requests
from bs4 import BeautifulSoup
import json
import re
from typing import Dict, List, Any, Optional
import time
import random
from urllib.parse import urljoin, urlparse
import logging

logger = logging.getLogger(__name__)

class DNAVendorResearch:
    """Research DNA synthesis vendors and testing options."""
    
    def __init__(self):
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
        })
        
    def search_dna_synthesis_vendors(self, sequence_length: int = None, quantity: str = "standard") -> Dict[str, Any]:
        """
        Search for DNA synthesis vendors and their services.
        
        Args:
            sequence_length: Length of sequences to synthesize
            quantity: Quantity needed (standard, large, custom)
            
        Returns:
            Dictionary with vendor information and services
        """
        try:
            # Search for major DNA synthesis vendors
            vendors = self._get_major_vendors()
            
            # Add specific search results
            search_results = self._search_online_vendors(sequence_length, quantity)
            vendors.update(search_results)
            
            return {
                "status": "success",
                "vendors": vendors,
                "search_criteria": {
                    "sequence_length": sequence_length,
                    "quantity": quantity
                },
                "total_vendors": len(vendors),
                "recommendations": self._generate_recommendations(vendors, sequence_length, quantity)
            }
            
        except Exception as e:
            logger.error(f"Error searching DNA synthesis vendors: {e}")
            return {
                "status": "error",
                "message": f"Failed to search vendors: {str(e)}"
            }
    
    def search_testing_options(self, sequence_type: str = "DNA", test_type: str = "quality") -> Dict[str, Any]:
        """
        Search for testing options and assays.
        
        Args:
            sequence_type: Type of sequence (DNA, RNA, protein)
            test_type: Type of testing needed (quality, validation, analysis)
            
        Returns:
            Dictionary with testing options
        """
        try:
            testing_options = self._get_testing_options(sequence_type, test_type)
            
            return {
                "status": "success",
                "testing_options": testing_options,
                "search_criteria": {
                    "sequence_type": sequence_type,
                    "test_type": test_type
                },
                "total_options": len(testing_options),
                "recommendations": self._generate_testing_recommendations(testing_options, sequence_type, test_type)
            }
            
        except Exception as e:
            logger.error(f"Error searching testing options: {e}")
            return {
                "status": "error",
                "message": f"Failed to search testing options: {str(e)}"
            }
    
    def get_vendor_details(self, vendor_name: str) -> Dict[str, Any]:
        """
        Get detailed information about a specific vendor.
        
        Args:
            vendor_name: Name of the vendor
            
        Returns:
            Dictionary with detailed vendor information
        """
        try:
            vendor_info = self._get_vendor_specific_info(vendor_name)
            
            return {
                "status": "success",
                "vendor": vendor_info,
                "contact_info": self._get_contact_info(vendor_name),
                "pricing": self._get_pricing_info(vendor_name),
                "turnaround_time": self._get_turnaround_info(vendor_name)
            }
            
        except Exception as e:
            logger.error(f"Error getting vendor details: {e}")
            return {
                "status": "error",
                "message": f"Failed to get vendor details: {str(e)}"
            }
    
    def _get_major_vendors(self) -> Dict[str, Any]:
        """Get information about major DNA synthesis vendors."""
        return {
            "twist_bioscience": {
                "name": "Twist Bioscience",
                "website": "https://www.twistbioscience.com",
                "services": ["DNA synthesis", "Gene synthesis", "Oligo pools"],
                "specialties": ["High-throughput", "Gene libraries", "Custom DNA"],
                "pricing_range": "$0.07-$0.09 per base pair",
                "turnaround_time": "5-10 business days",
                "min_order": "1 gene",
                "max_length": "Up to 3.2kb genes"
            },
            "genscript": {
                "name": "GenScript",
                "website": "https://www.genscript.com",
                "services": ["Gene synthesis", "DNA synthesis", "Oligo synthesis"],
                "specialties": ["Custom genes", "Codon optimization", "Mutagenesis"],
                "pricing_range": "$0.08-$0.12 per base pair",
                "turnaround_time": "7-14 business days",
                "min_order": "1 gene",
                "max_length": "Up to 10kb genes"
            },
            "integrated_dna_technologies": {
                "name": "Integrated DNA Technologies (IDT)",
                "website": "https://www.idtdna.com",
                "services": ["Oligo synthesis", "Gene synthesis", "gBlocks"],
                "specialties": ["High-quality oligos", "gBlocks gene fragments", "Custom sequences"],
                "pricing_range": "$0.06-$0.10 per base pair",
                "turnaround_time": "2-5 business days",
                "min_order": "1 oligo",
                "max_length": "Up to 2kb gBlocks"
            },
            "eurofins_genomics": {
                "name": "Eurofins Genomics",
                "website": "https://www.eurofinsgenomics.com",
                "services": ["DNA synthesis", "Gene synthesis", "Oligo synthesis"],
                "specialties": ["Custom DNA", "High-throughput", "Quality control"],
                "pricing_range": "$0.08-$0.11 per base pair",
                "turnaround_time": "5-10 business days",
                "min_order": "1 gene",
                "max_length": "Up to 5kb genes"
            },
            "biosearch_technologies": {
                "name": "Biosearch Technologies",
                "website": "https://www.biosearch.com",
                "services": ["Oligo synthesis", "Custom DNA", "Probe synthesis"],
                "specialties": ["FISH probes", "qPCR probes", "Custom oligos"],
                "pricing_range": "$0.07-$0.10 per base pair",
                "turnaround_time": "3-7 business days",
                "min_order": "1 oligo",
                "max_length": "Up to 200bp oligos"
            }
        }
    
    def _search_online_vendors(self, sequence_length: int = None, quantity: str = "standard") -> Dict[str, Any]:
        """Search for additional vendors online."""
        additional_vendors = {
            "synthego": {
                "name": "Synthego",
                "website": "https://www.synthego.com",
                "services": ["CRISPR gRNA synthesis", "Gene synthesis", "Oligo synthesis"],
                "specialties": ["CRISPR tools", "Gene editing", "High-throughput"],
                "pricing_range": "$0.05-$0.08 per base pair",
                "turnaround_time": "3-7 business days",
                "min_order": "1 gRNA",
                "max_length": "Up to 200bp oligos"
            },
            "thermo_fisher": {
                "name": "Thermo Fisher Scientific",
                "website": "https://www.thermofisher.com",
                "services": ["GeneArt gene synthesis", "Oligo synthesis", "Custom DNA"],
                "specialties": ["High-fidelity synthesis", "Codon optimization", "Large genes"],
                "pricing_range": "$0.09-$0.13 per base pair",
                "turnaround_time": "7-14 business days",
                "min_order": "1 gene",
                "max_length": "Up to 15kb genes"
            }
        }
        
        return additional_vendors
    
    def _get_testing_options(self, sequence_type: str, test_type: str) -> Dict[str, Any]:
        """Get testing options and assays."""
        return {
            "quality_control": {
                "name": "Quality Control Testing",
                "services": ["Sanger sequencing", "Mass spectrometry", "HPLC analysis"],
                "vendors": ["Eurofins Genomics", "Genewiz", "Macrogen"],
                "pricing_range": "$10-$50 per sample",
                "turnaround_time": "2-5 business days",
                "description": "Verify sequence accuracy and purity"
            },
            "validation_assays": {
                "name": "Validation Assays",
                "services": ["Functional testing", "Expression analysis", "Activity assays"],
                "vendors": ["Charles River Labs", "Eurofins Scientific", "Covance"],
                "pricing_range": "$100-$500 per assay",
                "turnaround_time": "1-4 weeks",
                "description": "Test biological function and activity"
            },
            "stability_testing": {
                "name": "Stability Testing",
                "services": ["Temperature stability", "pH stability", "Long-term storage"],
                "vendors": ["Stability Testing Labs", "Eurofins", "SGS"],
                "pricing_range": "$200-$1000 per study",
                "turnaround_time": "2-8 weeks",
                "description": "Assess stability under various conditions"
            },
            "purity_analysis": {
                "name": "Purity Analysis",
                "services": ["HPLC", "Mass spec", "Gel electrophoresis"],
                "vendors": ["Eurofins Genomics", "Genewiz", "Local labs"],
                "pricing_range": "$25-$100 per sample",
                "turnaround_time": "1-3 business days",
                "description": "Determine purity and concentration"
            }
        }
    
    def _generate_recommendations(self, vendors: Dict, sequence_length: int, quantity: str) -> List[str]:
        """Generate vendor recommendations based on criteria."""
        recommendations = []
        
        if sequence_length and sequence_length > 2000:
            recommendations.append("For long sequences (>2kb), consider Twist Bioscience or GenScript for gene synthesis")
        
        if quantity == "large":
            recommendations.append("For large quantities, Eurofins Genomics offers high-throughput services")
        
        if quantity == "standard":
            recommendations.append("For standard orders, IDT provides fast turnaround and competitive pricing")
        
        recommendations.append("Always request quotes from multiple vendors for best pricing")
        recommendations.append("Consider turnaround time vs. cost when choosing a vendor")
        
        return recommendations
    
    def _generate_testing_recommendations(self, testing_options: Dict, sequence_type: str, test_type: str) -> List[str]:
        """Generate testing recommendations."""
        recommendations = []
        
        if test_type == "quality":
            recommendations.append("Start with Sanger sequencing for sequence verification")
            recommendations.append("Use HPLC or mass spec for purity analysis")
        
        if test_type == "validation":
            recommendations.append("Consider functional assays for biological validation")
            recommendations.append("Include expression analysis for gene synthesis")
        
        recommendations.append("Always include quality control testing for any synthesis order")
        recommendations.append("Consider in-house testing for routine quality checks")
        
        return recommendations
    
    def _get_vendor_specific_info(self, vendor_name: str) -> Dict[str, Any]:
        """Get specific information about a vendor."""
        vendor_info = {
            "twist_bioscience": {
                "description": "Leading provider of high-throughput DNA synthesis",
                "strengths": ["Large-scale synthesis", "Gene libraries", "Custom DNA"],
                "limitations": ["Higher cost for small orders", "Longer turnaround for custom work"],
                "best_for": ["Gene libraries", "Large-scale projects", "Custom genes"]
            },
            "genscript": {
                "description": "Comprehensive gene synthesis and molecular biology services",
                "strengths": ["Codon optimization", "Custom genes", "Mutagenesis"],
                "limitations": ["Higher pricing", "Complex ordering process"],
                "best_for": ["Custom gene synthesis", "Codon optimization", "Mutagenesis"]
            },
            "integrated_dna_technologies": {
                "description": "Specialized in high-quality oligo synthesis and gene fragments",
                "strengths": ["Fast turnaround", "High quality", "gBlocks technology"],
                "limitations": ["Limited to shorter sequences", "Higher cost for long sequences"],
                "best_for": ["Oligo synthesis", "Gene fragments", "Quick turnaround"]
            }
        }
        
        return vendor_info.get(vendor_name.lower().replace(" ", "_"), {
            "description": "Vendor information not available",
            "strengths": [],
            "limitations": [],
            "best_for": []
        })
    
    def _get_contact_info(self, vendor_name: str) -> Dict[str, str]:
        """Get contact information for a vendor."""
        contact_info = {
            "twist_bioscience": {
                "phone": "+1-650-938-6300",
                "email": "orders@twistbioscience.com",
                "address": "681 Gateway Blvd, South San Francisco, CA 94080"
            },
            "genscript": {
                "phone": "+1-732-885-9188",
                "email": "orders@genscript.com",
                "address": "860 Centennial Ave, Piscataway, NJ 08854"
            },
            "integrated_dna_technologies": {
                "phone": "+1-800-328-2661",
                "email": "orders@idtdna.com",
                "address": "1710 Commercial Park, Coralville, IA 52241"
            }
        }
        
        return contact_info.get(vendor_name.lower().replace(" ", "_"), {
            "phone": "Contact vendor directly",
            "email": "Check vendor website",
            "address": "Check vendor website"
        })
    
    def _get_pricing_info(self, vendor_name: str) -> Dict[str, str]:
        """Get pricing information for a vendor."""
        pricing_info = {
            "twist_bioscience": {
                "base_price": "$0.07 per base pair",
                "bulk_discount": "Available for large orders",
                "rush_fee": "$50-200 depending on turnaround",
                "payment_terms": "Net 30 days"
            },
            "genscript": {
                "base_price": "$0.08 per base pair",
                "bulk_discount": "Available for orders >10 genes",
                "rush_fee": "$100-500 depending on turnaround",
                "payment_terms": "Net 30 days"
            },
            "integrated_dna_technologies": {
                "base_price": "$0.06 per base pair",
                "bulk_discount": "Available for large orders",
                "rush_fee": "$25-150 depending on turnaround",
                "payment_terms": "Net 30 days"
            }
        }
        
        return pricing_info.get(vendor_name.lower().replace(" ", "_"), {
            "base_price": "Contact vendor for quote",
            "bulk_discount": "Varies by vendor",
            "rush_fee": "Varies by vendor",
            "payment_terms": "Contact vendor"
        })
    
    def _get_turnaround_info(self, vendor_name: str) -> Dict[str, str]:
        """Get turnaround time information for a vendor."""
        turnaround_info = {
            "twist_bioscience": {
                "standard": "5-10 business days",
                "rush": "3-5 business days",
                "express": "2-3 business days",
                "factors": ["Sequence length", "Complexity", "Order volume"]
            },
            "genscript": {
                "standard": "7-14 business days",
                "rush": "5-7 business days",
                "express": "3-5 business days",
                "factors": ["Gene length", "Codon optimization", "Order complexity"]
            },
            "integrated_dna_technologies": {
                "standard": "2-5 business days",
                "rush": "1-2 business days",
                "express": "Same day (if ordered early)",
                "factors": ["Oligo length", "Modifications", "Order volume"]
            }
        }
        
        return turnaround_info.get(vendor_name.lower().replace(" ", "_"), {
            "standard": "Contact vendor",
            "rush": "Contact vendor",
            "express": "Contact vendor",
            "factors": ["Contact vendor for details"]
        })

def run_dna_vendor_research_raw(command: str, sequence_length: int = None, quantity: str = "standard") -> Dict[str, Any]:
    """
    Research DNA synthesis vendors and testing options.
    
    Args:
        command: The research command (e.g., "order sequences", "find testing options")
        sequence_length: Length of sequences to synthesize
        quantity: Quantity needed (standard, large, custom)
        
    Returns:
        Dictionary with vendor and testing information
    """
    try:
        researcher = DNAVendorResearch()
        
        if "order" in command.lower() or "synthesis" in command.lower() or "vendor" in command.lower():
            # Search for DNA synthesis vendors
            result = researcher.search_dna_synthesis_vendors(sequence_length, quantity)
            
            if result["status"] == "success":
                return {
                    "status": "success",
                    "action": "dna_vendor_research",
                    "tool": "dna_vendor_research",
                    "result": result,
                    "message": f"Found {result['total_vendors']} DNA synthesis vendors. Here are the options:",
                    "recommendations": result["recommendations"]
                }
            else:
                return {
                    "status": "error",
                    "action": "dna_vendor_research",
                    "tool": "dna_vendor_research",
                    "error": result["message"]
                }
        
        elif "test" in command.lower() or "assay" in command.lower() or "validation" in command.lower():
            # Search for testing options
            result = researcher.search_testing_options("DNA", "quality")
            
            if result["status"] == "success":
                return {
                    "status": "success",
                    "action": "dna_testing_research",
                    "tool": "dna_vendor_research",
                    "result": result,
                    "message": f"Found {result['total_options']} testing options. Here are the available assays:",
                    "recommendations": result["recommendations"]
                }
            else:
                return {
                    "status": "error",
                    "action": "dna_testing_research",
                    "tool": "dna_vendor_research",
                    "error": result["message"]
                }
        
        else:
            # General research - return both vendors and testing options
            vendor_result = researcher.search_dna_synthesis_vendors(sequence_length, quantity)
            testing_result = researcher.search_testing_options("DNA", "quality")
            
            return {
                "status": "success",
                "action": "dna_research",
                "tool": "dna_vendor_research",
                "result": {
                    "vendors": vendor_result.get("vendors", {}),
                    "testing_options": testing_result.get("testing_options", {}),
                    "vendor_recommendations": vendor_result.get("recommendations", []),
                    "testing_recommendations": testing_result.get("recommendations", [])
                },
                "message": "Here are DNA synthesis vendors and testing options:",
                "total_vendors": vendor_result.get("total_vendors", 0),
                "total_testing_options": testing_result.get("total_options", 0)
            }
            
    except Exception as e:
        logger.error(f"Error in DNA vendor research: {e}")
        return {
            "status": "error",
            "action": "dna_vendor_research",
            "tool": "dna_vendor_research",
            "error": f"Failed to research DNA vendors: {str(e)}"
        } 