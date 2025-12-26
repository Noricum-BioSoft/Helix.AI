#!/bin/bash
# Quick script to download and visualize FASTQ results from S3

set -e

if [ $# -lt 1 ]; then
    echo "Usage: $0 <s3-path-to-results.json>"
    echo "Example: $0 s3://noricum-ngs-data/fastqc-results/GRCh38.p12.MafHi/results.json"
    exit 1
fi

S3_PATH="$1"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

echo "=========================================="
echo "FASTQ Results Visualizer"
echo "=========================================="
echo "S3 Path: $S3_PATH"
echo ""

# Download results
LOCAL_RESULTS="/tmp/fastqc_results_$$.json"
echo "Downloading results from S3..."
aws s3 cp "$S3_PATH" "$LOCAL_RESULTS" || {
    echo "Error: Failed to download results from $S3_PATH"
    exit 1
}
echo "✅ Results downloaded"

# Create visualizations
echo "Creating visualizations..."
cd "$PROJECT_ROOT"
python3 "$SCRIPT_DIR/visualize_fastqc_results.py" \
    --input "$LOCAL_RESULTS" \
    --output-dir "fastqc_viz"

# Clean up
rm -f "$LOCAL_RESULTS"

echo ""
echo "=========================================="
echo "✅ Visualization complete!"
echo "=========================================="
echo "Open fastqc_viz/fastqc_results.html in your browser"
echo ""

# Try to open in browser (macOS)
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "Opening in browser..."
    open "fastqc_viz/fastqc_results.html"
fi










