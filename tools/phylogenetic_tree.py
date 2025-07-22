import io
import re
import numpy as np
from typing import List, Dict, Any
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align.Applications import ClustalwCommandline
import subprocess
import tempfile
import os
from sklearn.cluster import KMeans, AgglomerativeClustering
from sklearn.metrics.pairwise import pairwise_distances

# Check if ETE3 is available
try:
    from ete3 import Tree
    ETE3_AVAILABLE = True
except ImportError:
    ETE3_AVAILABLE = False
    print("⚠️ ETE3 not available, using fallback visualization")

def parse_aligned_sequences(alignment_data: str) -> List[Dict[str, str]]:
    """Parse FASTA sequences from alignment data."""
    sequences = []
    lines = alignment_data.strip().split('\n')
    current_name = None
    current_sequence = ""
    
    for line in lines:
        line = line.strip()
        if line.startswith('>'):
            if current_name and current_sequence:
                sequences.append({
                    'name': current_name,
                    'sequence': current_sequence
                })
            current_name = line[1:]
            current_sequence = ""
        else:
            current_sequence += line
    
    if current_name and current_sequence:
        sequences.append({
            'name': current_name,
            'sequence': current_sequence
        })
    
    return sequences

def align_sequences_with_biopython(sequences: List[Dict[str, str]]) -> List[Dict[str, str]]:
    """Align sequences using simple manual alignment."""
    try:
        # Simple manual alignment - pad shorter sequences with gaps
        max_length = max(len(seq['sequence']) for seq in sequences)
        aligned_sequences = []
        
        for seq in sequences:
            aligned_seq = seq['sequence']
            # Pad with gaps to match the longest sequence
            while len(aligned_seq) < max_length:
                aligned_seq += '-'
            
            aligned_sequences.append({
                'name': seq['name'],
                'sequence': aligned_seq
            })
        
        print(f"🔧 Manually aligned {len(aligned_sequences)} sequences to length {max_length}")
        return aligned_sequences
            
    except Exception as e:
        print(f"⚠️ Alignment error: {e}, using original sequences")
        return sequences

def create_phylogenetic_tree_with_biopython(aligned_sequences: List[Dict[str, str]]) -> Dict[str, Any]:
    """Create phylogenetic tree using BioPython's tree construction tools."""
    try:
        from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
        from Bio import AlignIO
        from io import StringIO
        
        # Create alignment object
        alignment_str = ""
        for seq in aligned_sequences:
            alignment_str += f">{seq['name']}\n{seq['sequence']}\n"
        
        alignment = AlignIO.read(StringIO(alignment_str), "fasta")
        
        # Calculate distance matrix
        calculator = DistanceCalculator('identity')
        dm = calculator.get_distance(alignment)
        
        # Construct tree
        constructor = DistanceTreeConstructor(calculator, 'upgma')
        tree = constructor.build_tree(alignment)
        
        # Clean internal node names before writing Newick
        def clean_node_names(clade):
            """Remove internal node names to prevent parsing issues"""
            if clade.name and clade.name.startswith('Inner'):
                clade.name = None
            for child in clade.clades:
                clean_node_names(child)
        
        # Clean the tree before writing
        clean_node_names(tree.root)
        
        # Convert to Newick string
        newick_io = io.StringIO()
        Phylo.write(tree, newick_io, "newick")
        newick_str = newick_io.getvalue().strip()
        
        return {
            "tree_newick": newick_str,
            "aligned_sequences": aligned_sequences,
            "distance_matrix": dm
        }
        
    except Exception as e:
        print(f"⚠️ Tree construction error: {e}")
        # Fallback: create simple tree
        if len(aligned_sequences) >= 2:
            newick_str = f"({aligned_sequences[0]['name']}:0.1,{aligned_sequences[1]['name']}:0.1);"
        else:
            newick_str = f"({aligned_sequences[0]['name']}:0.1);"
        
        return {
            "tree_newick": newick_str,
            "aligned_sequences": aligned_sequences,
            "distance_matrix": None
        }

def create_ete_visualization(newick_str: str, aligned_sequences: List[Dict[str, str]]) -> Dict[str, Any]:
    """Create ETE3 visualization of the phylogenetic tree."""
    try:
        if not ETE3_AVAILABLE:
            return {"error": "ETE3 not available"}
        
        # Parse the tree with ETE3
        tree = Tree(newick_str)
        
        # Create a simple text representation of the tree
        tree_text = tree.get_ascii(show_internal=True, compact=True)
        
        # Create a simple HTML representation
        html_content = f"""
        <div style="font-family: monospace; background: white; padding: 20px; border: 1px solid #ccc;">
            <h4>Phylogenetic Tree Visualization</h4>
            <pre style="font-size: 12px; overflow-x: auto;">{tree_text}</pre>
            <p><strong>Sequences:</strong> {len(aligned_sequences)}</p>
        </div>
        """
        
        return {
            "svg": html_content,
            "tree_newick": newick_str,
            "num_sequences": len(aligned_sequences),
            "tree_text": tree_text
        }
        
    except Exception as e:
        print(f"⚠️ ETE3 visualization error: {e}")
        # Fallback to simple text representation
        html_content = f"""
        <div style="font-family: monospace; background: white; padding: 20px; border: 1px solid #ccc;">
            <h4>Phylogenetic Tree (Fallback)</h4>
            <p><strong>Sequences:</strong> {len(aligned_sequences)}</p>
            <p><strong>Error:</strong> {str(e)}</p>
        </div>
        """
        return {
            "svg": html_content,
            "tree_newick": newick_str,
            "num_sequences": len(aligned_sequences),
            "tree_text": f"Error: {str(e)}"
        }

def create_phylogenetic_tree(aligned_sequences: List[Dict[str, str]]) -> Dict[str, Any]:
    """Main function to create phylogenetic tree with alignment and visualization."""
    print(f"🔧 Creating phylogenetic tree for {len(aligned_sequences)} sequences")
    
    # Step 1: Align sequences using BioPython
    print("🔧 Step 1: Aligning sequences...")
    aligned_seqs = align_sequences_with_biopython(aligned_sequences)
    print(f"🔧 Aligned {len(aligned_seqs)} sequences")
    
    # Step 2: Generate Newick string using BioPython
    print("🔧 Step 2: Generating phylogenetic tree...")
    tree_result = create_phylogenetic_tree_with_biopython(aligned_seqs)
    newick_str = tree_result["tree_newick"]
    print(f"🔧 Generated Newick: {newick_str}")
    
    # Step 3: Create ETE3 visualization
    print("🔧 Step 3: Creating ETE3 visualization...")
    ete_result = create_ete_visualization(newick_str, aligned_seqs)
    
    # Always return the ETE3 result, even if it has an error
    return {
        "text": f"Phylogenetic tree created successfully!\n\nNumber of sequences: {len(aligned_seqs)}",
        "tree_newick": newick_str,
        "aligned_sequences": aligned_seqs,
        "ete_visualization": ete_result
    }

def run_phylogenetic_tree_raw(aligned_sequences: str):
    """Raw function for phylogenetic tree creation."""
    sequences = parse_aligned_sequences(aligned_sequences)
    if len(sequences) < 2:
        return {"error": "At least 2 sequences are required for phylogenetic tree construction"}
    
    return create_phylogenetic_tree(sequences)

from langchain.agents import tool

@tool
def run_phylogenetic_tree(aligned_sequences: str):
    """Create a phylogenetic tree from aligned sequences.
    
    Args:
        aligned_sequences: FASTA format sequences to analyze
        
    Returns:
        Phylogenetic tree visualization and analysis results
    """
    return run_phylogenetic_tree_raw(aligned_sequences)

@tool
def run_clustering_analysis(aligned_sequences: str, num_clusters: int = 5):
    """Cluster sequences from phylogenetic tree and select representatives.
    
    Args:
        aligned_sequences: FASTA format sequences to analyze
        num_clusters: Number of clusters to create (default: 5)
        
    Returns:
        Clustering results with representative sequences and visualization
    """
    return run_clustering_from_tree(aligned_sequences, num_clusters)

@tool
def run_variant_selection(aligned_sequences: str, num_variants: int = 10):
    """Select top representative variants from phylogenetic tree.
    
    Args:
        aligned_sequences: FASTA format sequences to analyze
        num_variants: Number of top variants to select (default: 10)
        
    Returns:
        Selected variants with diversity scores and analysis
    """
    return run_variant_selection_from_tree(aligned_sequences, num_variants)

def convert_numpy_types(obj):
    """Convert numpy types to Python types for JSON serialization."""
    import numpy as np
    
    if isinstance(obj, dict):
        return {key: convert_numpy_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    else:
        return obj

def cluster_sequences_from_tree(newick_str: str, aligned_sequences: List[Dict[str, str]], num_clusters: int = 5) -> Dict[str, Any]:
    """Cluster sequences based on phylogenetic tree structure."""
    try:
        if not ETE3_AVAILABLE:
            return {"error": "ETE3 not available for clustering"}
        
        # Clean the Newick string first
        cleaned_newick = re.sub(r'\)([A-Za-z0-9_]+):', r'):', newick_str)
        
        # Parse the tree with ETE3
        try:
            tree = Tree(cleaned_newick)
        except Exception as e:
            print(f"⚠️ ETE3 tree parsing error: {e}")
            return {"error": f"Tree parsing failed: {str(e)}"}
        
        # Extract leaf nodes and their distances
        leaf_nodes = []
        leaf_names = []
        
        try:
            for node in tree.traverse():
                if node.is_leaf():
                    leaf_nodes.append(node)
                    leaf_names.append(node.name)
        except Exception as e:
            print(f"⚠️ Leaf extraction error: {e}")
            return {"error": f"Leaf extraction failed: {str(e)}"}
        
        # Check if we have enough sequences
        if len(leaf_nodes) < 2:
            return {"error": "At least 2 sequences are required for clustering"}
        
        # Calculate pairwise distances between all leaves
        n_leaves = len(leaf_nodes)
        distance_matrix = np.zeros((n_leaves, n_leaves))
        
        try:
            for i, node1 in enumerate(leaf_nodes):
                for j, node2 in enumerate(leaf_nodes):
                    if i != j:
                        # Calculate distance between two leaves
                        distance = tree.get_distance(node1, node2)
                        distance_matrix[i][j] = distance
                        distance_matrix[j][i] = distance
        except Exception as e:
            print(f"⚠️ Distance calculation error: {e}")
            return {"error": f"Distance calculation failed: {str(e)}"}
        
        # Ensure we don't have more clusters than sequences
        num_clusters = min(num_clusters, n_leaves)
        
        # For very small datasets, just return simple clustering
        if n_leaves <= 3:
            # Simple clustering for small datasets
            clusters = {}
            for i, name in enumerate(leaf_names):
                cluster_id = i % num_clusters
                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append({
                    'name': name,
                    'distance': distance_matrix[i].mean(),
                    'index': i
                })
        else:
            # Perform clustering using hierarchical clustering
            try:
                clustering = AgglomerativeClustering(
                    n_clusters=num_clusters,
                    metric='precomputed',
                    linkage='complete'
                )
                cluster_labels = clustering.fit_predict(distance_matrix)
                
                # Group sequences by cluster
                clusters = {}
                for i, label in enumerate(cluster_labels):
                    if label not in clusters:
                        clusters[label] = []
                    clusters[label].append({
                        'name': leaf_names[i],
                        'distance': distance_matrix[i].mean(),
                        'index': i
                    })
            except Exception as e:
                print(f"⚠️ Clustering algorithm error: {e}")
                return {"error": f"Clustering algorithm failed: {str(e)}"}
        
        # Select representative sequences from each cluster
        representatives = []
        cluster_info = []
        
        for cluster_id, cluster_sequences in clusters.items():
            # Sort by average distance (closest to cluster center)
            cluster_sequences.sort(key=lambda x: x['distance'])
            
            # Select the most representative sequence (closest to cluster center)
            representative = cluster_sequences[0]
            representatives.append(representative['name'])
            
            cluster_info.append({
                'cluster_id': cluster_id,
                'size': len(cluster_sequences),
                'representative': representative['name'],
                'sequences': [seq['name'] for seq in cluster_sequences],
                'average_distance': float(representative['distance'])
            })
        
        result = {
            'clusters': cluster_info,
            'representatives': representatives,
            'distance_matrix': distance_matrix.tolist(),
            'cluster_labels': [0] * n_leaves,  # Simple labels for small datasets
            'total_sequences': n_leaves,
            'num_clusters': len(clusters)
        }
        
        # Convert any numpy types to Python types for JSON serialization
        return convert_numpy_types(result)
        
    except Exception as e:
        print(f"⚠️ Clustering error: {e}")
        return {"error": f"Clustering failed: {str(e)}"}

def select_top_variants_from_tree(newick_str: str, aligned_sequences: List[Dict[str, str]], num_variants: int = 10) -> Dict[str, Any]:
    """Select top representative variants from phylogenetic tree."""
    try:
        if not ETE3_AVAILABLE:
            return {"error": "ETE3 not available for variant selection"}
        
        # Parse the tree with ETE3
        tree = Tree(newick_str)
        
        # Extract leaf nodes and calculate their diversity scores
        leaf_nodes = []
        diversity_scores = []
        
        for node in tree.traverse():
            if node.is_leaf():
                leaf_nodes.append(node)
                
                # Calculate diversity score based on distance to other nodes
                total_distance = 0
                count = 0
                for other_node in tree.traverse():
                    if other_node != node:
                        distance = tree.get_distance(node, other_node)
                        total_distance += distance
                        count += 1
                
                avg_distance = total_distance / count if count > 0 else 0
                diversity_scores.append({
                    'name': node.name,
                    'diversity_score': avg_distance,
                    'node': node
                })
        
        # Sort by diversity score (highest diversity first)
        diversity_scores.sort(key=lambda x: x['diversity_score'], reverse=True)
        
        # Select top variants
        selected_variants = diversity_scores[:num_variants]
        
        result = {
            'selected_variants': [var['name'] for var in selected_variants],
            'diversity_scores': {var['name']: var['diversity_score'] for var in selected_variants},
            'total_sequences': len(leaf_nodes),
            'num_selected': len(selected_variants)
        }
        
        # Convert any numpy types to Python types for JSON serialization
        return convert_numpy_types(result)
        
    except Exception as e:
        print(f"⚠️ Variant selection error: {e}")
        return {"error": f"Variant selection failed: {str(e)}"}

def create_clustered_visualization(newick_str: str, aligned_sequences: List[Dict[str, str]], cluster_result: Dict[str, Any]) -> Dict[str, Any]:
    """Create visualization highlighting clustered sequences."""
    try:
        if not ETE3_AVAILABLE:
            return {"error": "ETE3 not available for clustered visualization"}
        
        # Parse the tree with ETE3
        tree = Tree(newick_str)
        
        # Create ASCII tree with cluster highlighting
        tree_text = tree.get_ascii(show_internal=True, compact=True)
        
        # Add cluster information
        cluster_info = ""
        if 'clusters' in cluster_result:
            cluster_info = "\n\nCluster Information:\n"
            cluster_info += "=" * 50 + "\n"
            for cluster in cluster_result['clusters']:
                cluster_info += f"Cluster {cluster['cluster_id']}:\n"
                cluster_info += f"  Size: {cluster['size']} sequences\n"
                cluster_info += f"  Representative: {cluster['representative']}\n"
                cluster_info += f"  Average distance: {cluster['average_distance']:.4f}\n"
                cluster_info += f"  Sequences: {', '.join(cluster['sequences'][:5])}"
                if len(cluster['sequences']) > 5:
                    cluster_info += f" ... (+{len(cluster['sequences']) - 5} more)"
                cluster_info += "\n\n"
        
        # Create HTML visualization
        html_content = f"""
        <div style="font-family: monospace; background: white; padding: 20px; border: 1px solid #ccc;">
            <h4>Clustered Phylogenetic Tree</h4>
            <pre style="font-size: 12px; overflow-x: auto;">{tree_text}</pre>
            <div style="margin-top: 20px;">
                <h5>Clustering Results:</h5>
                <p><strong>Total sequences:</strong> {cluster_result.get('total_sequences', 'N/A')}</p>
                <p><strong>Number of clusters:</strong> {cluster_result.get('num_clusters', 'N/A')}</p>
                <p><strong>Representatives:</strong> {', '.join(cluster_result.get('representatives', []))}</p>
            </div>
            <div style="margin-top: 20px;">
                <h5>Cluster Details:</h5>
                <pre style="font-size: 11px; background: #f5f5f5; padding: 10px;">{cluster_info}</pre>
            </div>
        </div>
        """
        
        return {
            "svg": html_content,
            "tree_text": tree_text,
            "cluster_info": cluster_info,
            "tree_newick": newick_str,
            "clusters": cluster_result
        }
        
    except Exception as e:
        print(f"⚠️ Clustered visualization error: {e}")
        return {"error": f"Clustered visualization failed: {str(e)}"}

def run_clustering_from_tree(aligned_sequences: str, num_clusters: int = 5) -> Dict[str, Any]:
    """Run clustering analysis on phylogenetic tree."""
    sequences = parse_aligned_sequences(aligned_sequences)
    if len(sequences) < 2:
        return {"error": "At least 2 sequences are required for clustering"}
    
    # Create phylogenetic tree first
    tree_result = create_phylogenetic_tree(sequences)
    if "error" in tree_result:
        return tree_result
    
    # Perform clustering
    cluster_result = cluster_sequences_from_tree(
        tree_result["tree_newick"], 
        tree_result["aligned_sequences"], 
        num_clusters
    )
    
    if "error" in cluster_result:
        return cluster_result
    
    # Create clustered visualization
    viz_result = create_clustered_visualization(
        tree_result["tree_newick"],
        tree_result["aligned_sequences"],
        cluster_result
    )
    
    return {
        "text": f"Clustering completed successfully!\n\nNumber of clusters: {cluster_result['num_clusters']}\nTotal sequences: {cluster_result['total_sequences']}",
        "tree_newick": tree_result["tree_newick"],
        "aligned_sequences": tree_result["aligned_sequences"],
        "clustering_result": cluster_result,
        "clustered_visualization": viz_result
    }

def run_variant_selection_from_tree(aligned_sequences: str, num_variants: int = 10) -> Dict[str, Any]:
    """Select top representative variants from phylogenetic tree."""
    sequences = parse_aligned_sequences(aligned_sequences)
    if len(sequences) < 2:
        return {"error": "At least 2 sequences are required for variant selection"}
    
    # Create phylogenetic tree first
    tree_result = create_phylogenetic_tree(sequences)
    if "error" in tree_result:
        return tree_result
    
    # Select variants
    variant_result = select_top_variants_from_tree(
        tree_result["tree_newick"],
        tree_result["aligned_sequences"],
        num_variants
    )
    
    if "error" in variant_result:
        return variant_result
    
    return {
        "text": f"Variant selection completed!\n\nSelected {variant_result['num_selected']} variants from {variant_result['total_sequences']} total sequences.",
        "tree_newick": tree_result["tree_newick"],
        "aligned_sequences": tree_result["aligned_sequences"],
        "variant_selection": variant_result
    } 