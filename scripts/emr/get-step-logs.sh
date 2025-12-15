#!/bin/bash
# Get detailed logs for a failed EMR step

set -e

CLUSTER_ID="${1:-}"
STEP_ID="${2:-}"
REGION="${AWS_REGION:-us-east-1}"
S3_BUCKET="${S3_DATASET_BUCKET:-noricum-ngs-data}"

if [ -z "$CLUSTER_ID" ] || [ -z "$STEP_ID" ]; then
    echo "Usage: $0 <cluster-id> <step-id>"
    echo ""
    echo "Example:"
    echo "  $0 j-12QYDE51Q9LDP s-0531929FDXFHEGGF63R"
    exit 1
fi

echo "=========================================="
echo "Fetching EMR Step Logs"
echo "=========================================="
echo "Cluster ID: $CLUSTER_ID"
echo "Step ID: $STEP_ID"
echo ""

# Create temp directory for logs
LOG_DIR="/tmp/emr-logs/$CLUSTER_ID/$STEP_ID"
mkdir -p "$LOG_DIR"

# Get step status
echo "Step Status:"
aws emr describe-step \
    --cluster-id "$CLUSTER_ID" \
    --step-id "$STEP_ID" \
    --region "$REGION" \
    --query 'Step.Status' \
    --output json | jq '.'

echo ""
echo "Step Configuration:"
aws emr describe-step \
    --cluster-id "$CLUSTER_ID" \
    --step-id "$STEP_ID" \
    --region "$REGION" \
    --query 'Step.Config' \
    --output json | jq '.'

echo ""
echo "Downloading logs from S3..."

# Download all available logs
aws s3 sync \
    "s3://${S3_BUCKET}/emr-logs/$CLUSTER_ID/steps/$STEP_ID/" \
    "$LOG_DIR/" \
    --region "$REGION" || echo "Warning: Some logs may not be available yet"

# Decompress and display logs
echo ""
echo "=========================================="
echo "Controller Log:"
echo "=========================================="
if [ -f "$LOG_DIR/controller.gz" ]; then
    gunzip -c "$LOG_DIR/controller.gz" | tail -100
else
    echo "Controller log not found"
fi

echo ""
echo "=========================================="
echo "STDOUT:"
echo "=========================================="
if [ -f "$LOG_DIR/stdout.gz" ]; then
    gunzip -c "$LOG_DIR/stdout.gz"
elif [ -f "$LOG_DIR/stdout" ]; then
    cat "$LOG_DIR/stdout"
else
    echo "STDOUT log not found"
fi

echo ""
echo "=========================================="
echo "STDERR:"
echo "=========================================="
if [ -f "$LOG_DIR/stderr.gz" ]; then
    gunzip -c "$LOG_DIR/stderr.gz"
elif [ -f "$LOG_DIR/stderr" ]; then
    cat "$LOG_DIR/stderr"
else
    echo "STDERR log not found"
    echo ""
    echo "Checking for YARN application logs..."
    
    # Try to get YARN application ID from step
    APP_ID=$(aws emr describe-step \
        --cluster-id "$CLUSTER_ID" \
        --step-id "$STEP_ID" \
        --region "$REGION" \
        --query 'Step.Status.Timeline' \
        --output json | jq -r '.StartDateTime' | tr -d ':-' | cut -c1-13) || true
    
    if [ -n "$APP_ID" ]; then
        echo "Looking for YARN application logs with timestamp: $APP_ID"
        echo "You may need to SSH into the master node to view YARN logs:"
        MASTER_DNS=$(aws emr describe-cluster \
            --cluster-id "$CLUSTER_ID" \
            --region "$REGION" \
            --query 'Cluster.MasterPublicDnsName' \
            --output text)
        echo "  ssh hadoop@$MASTER_DNS"
        echo "  yarn logs -applicationId <app-id>"
    fi
fi

echo ""
echo "=========================================="
echo "All log files:"
echo "=========================================="
ls -lh "$LOG_DIR/"

echo ""
echo "Logs saved to: $LOG_DIR"







