#!/bin/bash
#
# Deploy fixed FastQC analysis script to S3
#
# This script uploads the bugfix to the EMR scripts location in S3

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
SCRIPT_FILE="$SCRIPT_DIR/fastqc_analysis_v2.py"
S3_BUCKET="noricum-ngs-data"
S3_PATH="s3://$S3_BUCKET/emr-scripts/fastqc_analysis_v2.py"

echo "=============================================="
echo "Deploying Fixed FastQC Script to S3"
echo "=============================================="
echo ""

# Check if script exists
if [ ! -f "$SCRIPT_FILE" ]; then
    echo "❌ Error: Script file not found: $SCRIPT_FILE"
    exit 1
fi

echo "✅ Script found: $SCRIPT_FILE"
echo ""

# Verify Python syntax
echo "Verifying Python syntax..."
if ! python -m py_compile "$SCRIPT_FILE" 2>/dev/null; then
    echo "❌ Error: Script has syntax errors"
    exit 1
fi
echo "✅ Syntax valid"
echo ""

# Check AWS credentials
echo "Checking AWS credentials..."
if ! aws sts get-caller-identity >/dev/null 2>&1; then
    echo "❌ Error: AWS credentials not configured"
    echo "   Run: aws configure"
    exit 1
fi
echo "✅ AWS credentials configured"
echo ""

# Backup existing script
echo "Backing up existing script..."
BACKUP_PATH="s3://$S3_BUCKET/emr-scripts/backups/fastqc_analysis_v2_$(date +%Y%m%d_%H%M%S).py"
if aws s3 cp "$S3_PATH" "$BACKUP_PATH" 2>/dev/null; then
    echo "✅ Backup created: $BACKUP_PATH"
else
    echo "⚠️  No existing script to backup (this might be the first deployment)"
fi
echo ""

# Upload new script
echo "Uploading fixed script to S3..."
echo "  Source: $SCRIPT_FILE"
echo "  Destination: $S3_PATH"
echo ""

if aws s3 cp "$SCRIPT_FILE" "$S3_PATH"; then
    echo ""
    echo "✅ Script uploaded successfully!"
    echo ""
    echo "=============================================="
    echo "Deployment Complete"
    echo "=============================================="
    echo ""
    echo "Next steps:"
    echo "1. Test with large dataset:"
    echo "   python tests/demo_scenarios/run_fastqc_full_execution.py --large"
    echo ""
    echo "2. Monitor execution in EMR console"
    echo ""
    echo "3. Check that cleanup happens at the end (not beginning)"
    echo ""
    echo "4. Verify no FileNotFoundException errors"
    echo ""
else
    echo ""
    echo "❌ Upload failed"
    exit 1
fi
