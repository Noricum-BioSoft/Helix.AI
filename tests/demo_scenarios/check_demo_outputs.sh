#!/bin/bash
# Check Demo Outputs on S3

echo "======================================================================"
echo "  Checking Demo 3 Outputs"
echo "======================================================================"
echo ""

echo "1. Checking FastQC outputs (Step 1 - completed):"
echo "   Location: s3://noricum-ngs-data/results/complete-workflow/fastqc/"
aws s3 ls s3://noricum-ngs-data/results/complete-workflow/fastqc/
echo ""

echo "2. Checking alignment outputs (Step 2 - should be empty):"
echo "   Location: s3://noricum-ngs-data/results/complete-workflow/aligned/"
aws s3 ls s3://noricum-ngs-data/results/complete-workflow/aligned/ 2>&1 || echo "   [Empty or doesn't exist - expected, step 2 didn't run]"
echo ""

echo "3. Checking quantification outputs (Step 3 - should be empty):"
echo "   Location: s3://noricum-ngs-data/results/complete-workflow/counts/"
aws s3 ls s3://noricum-ngs-data/results/complete-workflow/counts/ 2>&1 || echo "   [Empty or doesn't exist - expected, step 3 didn't run]"
echo ""

echo "======================================================================"
echo "  Alternative: Check default FastQC location"
echo "======================================================================"
echo ""
echo "FastQC may have placed outputs in default location:"
echo "   s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/fastqc-results/"
aws s3 ls s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/test/fastqc-results/
echo ""

echo "======================================================================"
echo "  How to download FastQC reports"
echo "======================================================================"
echo ""
echo "If outputs exist, download them with:"
echo "  aws s3 cp s3://noricum-ngs-data/results/complete-workflow/fastqc/test_mate_R1_fastqc.html ."
echo "  open test_mate_R1_fastqc.html"
echo ""
