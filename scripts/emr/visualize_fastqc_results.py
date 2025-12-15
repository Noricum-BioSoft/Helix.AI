#!/usr/bin/env python3
"""
Visualize FASTQ Quality Analysis Results
Creates interactive visualizations from the EMR FASTQ analysis results.

Usage:
    python visualize_fastqc_results.py --input s3://bucket/path/results.json
    python visualize_fastqc_results.py --input results.json
"""

import argparse
import json
import sys
import os
from pathlib import Path

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import plotly.express as px
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Warning: Plotly not available. Install with: pip install plotly")

try:
    import boto3
    BOTO3_AVAILABLE = True
except ImportError:
    BOTO3_AVAILABLE = False


def download_from_s3(s3_path: str, local_path: str) -> str:
    """Download file from S3 to local path."""
    if not BOTO3_AVAILABLE:
        raise ImportError("boto3 is required for S3 downloads. Install with: pip install boto3")
    
    # Parse S3 path
    if not s3_path.startswith("s3://"):
        raise ValueError(f"Invalid S3 path: {s3_path}")
    
    path_parts = s3_path[5:].split("/", 1)
    bucket = path_parts[0]
    key = path_parts[1] if len(path_parts) > 1 else ""
    
    # Download
    s3_client = boto3.client('s3')
    s3_client.download_file(bucket, key, local_path)
    return local_path


def load_results(input_path: str) -> dict:
    """Load results from local file or S3."""
    if input_path.startswith("s3://"):
        print(f"Downloading results from S3: {input_path}")
        local_path = "/tmp/fastqc_results.json"
        download_from_s3(input_path, local_path)
        input_path = local_path
    
    with open(input_path, 'r') as f:
        return json.load(f)


def create_visualizations(results: dict, output_dir: str = "fastqc_visualizations"):
    """Create interactive visualizations from results."""
    if not PLOTLY_AVAILABLE:
        print("Error: Plotly is required for visualizations")
        print("Install with: pip install plotly")
        return None
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Read Count Comparison (R1 vs R2)',
            'Length Distribution',
            'GC Content Comparison',
            'Summary Statistics'
        ),
        specs=[[{"type": "bar"}, {"type": "bar"}],
               [{"type": "bar"}, {"type": "table"}]]
    )
    
    # 1. Read Count Comparison
    r1_count = results['basic_statistics']['r1']['total_sequences']
    r2_count = results['basic_statistics']['r2']['total_sequences']
    
    fig.add_trace(
        go.Bar(x=['R1', 'R2'], y=[r1_count, r2_count],
              marker_color=['#1f77b4', '#ff7f0e'],
              text=[f"{r1_count:,}", f"{r2_count:,}"],
              textposition='auto',
              name='Read Count'),
        row=1, col=1
    )
    
    # 2. Length Distribution
    r1_lengths = results['length_distribution']['r1']
    r2_lengths = results['length_distribution']['r2']
    
    # Convert to lists for plotting
    lengths = sorted([int(k) for k in r1_lengths.keys()])
    r1_counts = [r1_lengths.get(str(l), 0) for l in lengths]
    r2_counts = [r2_lengths.get(str(l), 0) for l in lengths]
    
    fig.add_trace(
        go.Bar(x=lengths, y=r1_counts, name='R1', marker_color='#1f77b4'),
        row=1, col=2
    )
    fig.add_trace(
        go.Bar(x=lengths, y=r2_counts, name='R2', marker_color='#ff7f0e'),
        row=1, col=2
    )
    
    # 3. GC Content Comparison
    r1_gc_mean = results['gc_content']['r1']['mean']
    r2_gc_mean = results['gc_content']['r2']['mean']
    r1_gc_std = results['gc_content']['r1']['std']
    r2_gc_std = results['gc_content']['r2']['std']
    
    fig.add_trace(
        go.Bar(x=['R1', 'R2'], y=[r1_gc_mean, r2_gc_mean],
              error_y=dict(type='data', array=[r1_gc_std, r2_gc_std]),
              marker_color=['#2ca02c', '#d62728'],
              text=[f"{r1_gc_mean:.2f}%", f"{r2_gc_mean:.2f}%"],
              textposition='auto',
              name='GC Content'),
        row=2, col=1
    )
    
    # 4. Summary Statistics Table
    r1_stats = results['basic_statistics']['r1']
    r2_stats = results['basic_statistics']['r2']
    
    fig.add_trace(
        go.Table(
            header=dict(values=['Metric', 'R1', 'R2'],
                       fill_color='paleturquoise',
                       align='left'),
            cells=dict(values=[
                ['Total Sequences', 'Total Bases', 'Average Length', 'GC Content Mean', 'GC Content Std'],
                [f"{r1_stats['total_sequences']:,}",
                 f"{r1_stats['total_bases']:,}",
                 f"{r1_stats['average_length']:.1f} bp",
                 f"{r1_gc_mean:.2f}%",
                 f"{r1_gc_std:.2f}%"],
                [f"{r2_stats['total_sequences']:,}",
                 f"{r2_stats['total_bases']:,}",
                 f"{r2_stats['average_length']:.1f} bp",
                 f"{r2_gc_mean:.2f}%",
                 f"{r2_gc_std:.2f}%"]
            ],
            fill_color='lavender',
            align='left')
        ),
        row=2, col=2
    )
    
    # Update layout
    fig.update_xaxes(title_text="Read Type", row=1, col=1)
    fig.update_yaxes(title_text="Number of Reads", row=1, col=1)
    fig.update_xaxes(title_text="Length (bp)", row=1, col=2)
    fig.update_yaxes(title_text="Count", row=1, col=2)
    fig.update_xaxes(title_text="Read Type", row=2, col=1)
    fig.update_yaxes(title_text="GC Content (%)", row=2, col=1)
    
    fig.update_layout(
        height=1000,
        title_text="FASTQ Quality Analysis Results",
        showlegend=True
    )
    
    # Save as HTML
    html_path = os.path.join(output_dir, "fastqc_results.html")
    fig.write_html(html_path)
    print(f"✅ Interactive visualization saved to: {html_path}")
    
    # Also create individual charts
    create_individual_charts(results, output_dir)
    
    return html_path


def create_individual_charts(results: dict, output_dir: str):
    """Create individual chart files."""
    # Read count comparison
    r1_count = results['basic_statistics']['r1']['total_sequences']
    r2_count = results['basic_statistics']['r2']['total_sequences']
    
    fig1 = go.Figure(data=[
        go.Bar(x=['R1', 'R2'], y=[r1_count, r2_count],
               marker_color=['#1f77b4', '#ff7f0e'],
               text=[f"{r1_count:,}", f"{r2_count:,}"],
               textposition='auto')
    ])
    fig1.update_layout(
        title='Read Count Comparison',
        xaxis_title='Read Type',
        yaxis_title='Number of Reads',
        height=400
    )
    fig1.write_html(os.path.join(output_dir, "read_count_comparison.html"))
    
    # GC content comparison with error bars
    r1_gc_mean = results['gc_content']['r1']['mean']
    r2_gc_mean = results['gc_content']['r2']['mean']
    r1_gc_std = results['gc_content']['r1']['std']
    r2_gc_std = results['gc_content']['r2']['std']
    
    fig2 = go.Figure(data=[
        go.Bar(x=['R1', 'R2'], y=[r1_gc_mean, r2_gc_mean],
               error_y=dict(type='data', array=[r1_gc_std, r2_gc_std]),
               marker_color=['#2ca02c', '#d62728'],
               text=[f"{r1_gc_mean:.2f}%", f"{r2_gc_mean:.2f}%"],
               textposition='auto')
    ])
    fig2.update_layout(
        title='GC Content Comparison',
        xaxis_title='Read Type',
        yaxis_title='GC Content (%)',
        height=400
    )
    fig2.write_html(os.path.join(output_dir, "gc_content_comparison.html"))
    
    print(f"✅ Individual charts saved to: {output_dir}")


def main():
    parser = argparse.ArgumentParser(description='Visualize FASTQ Quality Analysis Results')
    parser.add_argument('--input', required=True,
                       help='Path to results.json (local file or s3:// path)')
    parser.add_argument('--output-dir', default='fastqc_visualizations',
                       help='Output directory for visualizations (default: fastqc_visualizations)')
    
    args = parser.parse_args()
    
    try:
        # Load results
        print(f"Loading results from: {args.input}")
        results = load_results(args.input)
        
        # Create visualizations
        print("Creating visualizations...")
        html_path = create_visualizations(results, args.output_dir)
        
        if html_path:
            print(f"\n✅ Visualization complete!")
            print(f"Open {html_path} in your browser to view the results.")
        else:
            print("\n❌ Failed to create visualizations")
            return 1
        
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(main())







