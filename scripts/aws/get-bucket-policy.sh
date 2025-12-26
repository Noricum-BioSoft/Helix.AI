#!/bin/bash
# Get and display S3 bucket policy in a readable format
# Usage: ./get-bucket-policy.sh [bucket-name]

set -e

BUCKET_NAME="${1:-noricum-ngs-data}"

if [ -z "$BUCKET_NAME" ]; then
    echo "Usage: $0 [bucket-name]"
    echo "Example: $0 noricum-ngs-data"
    exit 1
fi

echo "Getting bucket policy for: $BUCKET_NAME"
echo ""

# Get bucket policy and format with jq
aws s3api get-bucket-policy \
    --bucket "$BUCKET_NAME" \
    --query 'Policy' \
    --output text | jq '.'


