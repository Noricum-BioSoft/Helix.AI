#!/bin/bash

# Monitor Docker Build Progress
# This script shows the current build stage and estimated progress

TERMINAL_FILE="/Users/eoberortner/.cursor/projects/Users-eoberortner-git-Helix-AI/terminals/393883.txt"

echo "🐳 Docker Build Monitor - Full Biotools Image"
echo "=============================================="
echo ""
echo "Building: helix-biotools:latest"
echo "Estimated time: ~20 minutes"
echo ""
echo "Build stages (22 total):"
echo "  1-3:   Base setup (Python, dependencies)"
echo "  4-8:   Core tools (FastQC, Muscle, Clustalo, Trimmomatic, SAMtools)"
echo "  9-14:  Aligners (BWA, Bowtie2, STAR, HISAT2, Subread)"
echo "  15:    R + Bioconductor (DESeq2, edgeR, limma, Seurat)"
echo "  16-18: Python requirements"
echo "  19-22: Copy application files"
echo ""

while true; do
    if [ ! -f "$TERMINAL_FILE" ]; then
        echo "⏳ Waiting for build to start..."
        sleep 2
        continue
    fi
    
    # Get current stage
    CURRENT_STAGE=$(tail -100 "$TERMINAL_FILE" | grep -oE '\[[0-9]+/22\]' | tail -1)
    
    # Get recent activity
    RECENT=$(tail -5 "$TERMINAL_FILE" | grep -v "^---" | grep -v "^$")
    
    # Check if build completed or failed
    if tail -20 "$TERMINAL_FILE" | grep -q "Successfully tagged helix-biotools:latest"; then
        echo ""
        echo "✅ BUILD COMPLETE!"
        echo ""
        docker images | grep helix-biotools
        exit 0
    elif tail -20 "$TERMINAL_FILE" | grep -q "ERROR:"; then
        echo ""
        echo "❌ BUILD FAILED!"
        echo ""
        tail -30 "$TERMINAL_FILE"
        exit 1
    fi
    
    clear
    echo "🐳 Docker Build Monitor - Full Biotools Image"
    echo "=============================================="
    echo ""
    echo "Current Stage: $CURRENT_STAGE"
    echo ""
    echo "Recent activity:"
    echo "$RECENT"
    echo ""
    echo "⏳ Building... (press Ctrl+C to stop monitoring, build continues in background)"
    
    sleep 3
done
