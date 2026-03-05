#!/bin/bash
# Upload GRCh38.p12.MafHi dataset to S3
# This script uploads the paired-end reads to the noricum-ngs-data bucket
# Optimized for large file uploads with retry logic

BUCKET="${S3_DATASET_BUCKET:-noricum-ngs-data}"
DATASET_ID="${1:-GRCh38.p12.MafHi}"
SOURCE_DIR="${2:-data/paired-end-reads/GRCh38.p12.MafHi}"
REGION="${AWS_REGION:-us-east-1}"
MAX_RETRIES=3

if [ ! -d "$SOURCE_DIR" ]; then
    echo "Error: Source directory not found: $SOURCE_DIR"
    echo "Usage: $0 [dataset_id] [source_directory]"
    echo "Example: $0 GRCh38.p12.MafHi data/paired-end-reads/GRCh38.p12.MafHi"
    exit 1
fi

echo "=========================================="
echo "Uploading dataset to S3"
echo "=========================================="
echo "Dataset ID: $DATASET_ID"
echo "Source: $SOURCE_DIR"
echo "Destination: s3://$BUCKET/datasets/$DATASET_ID/"
echo ""

# Check if AWS CLI is available
if ! command -v aws &> /dev/null; then
    echo "Error: AWS CLI not found. Please install it first."
    exit 1
fi

# Check AWS credentials
if ! aws sts get-caller-identity &> /dev/null; then
    echo "Error: AWS credentials not configured. Run 'aws configure' first."
    exit 1
fi

# Upload files with retry logic for large files
echo "Uploading files..."
echo "Note: Large files (5.8GB each) may take 10-30 minutes depending on connection speed."
echo "The upload will automatically retry on connection failures."
echo ""

# Configure AWS CLI for large file uploads
export AWS_MAX_ATTEMPTS=20
export AWS_RETRY_MODE=adaptive

# Note: For optimal large file uploads, configure multipart settings in ~/.aws/config:
# [default]
# s3 = 
#     multipart_threshold = 64MB
#     multipart_chunksize = 8MB      # Smaller chunks = more stable
#     max_concurrent_requests = 1     # Single request = more reliable
#
# This helps prevent "Connection was closed" errors during multipart uploads.

RETRY_COUNT=0
MAX_RETRIES=3
UPLOAD_SUCCESS=false

while [ $RETRY_COUNT -lt $MAX_RETRIES ] && [ "$UPLOAD_SUCCESS" = false ]; do
    if [ $RETRY_COUNT -gt 0 ]; then
        echo ""
        echo "Retry attempt $RETRY_COUNT of $MAX_RETRIES..."
        echo "Waiting 10 seconds before retry..."
        sleep 10
    fi
    
    # Use sync which automatically resumes interrupted uploads
    
    if aws s3 sync "$SOURCE_DIR/" \
      "s3://$BUCKET/datasets/$DATASET_ID/" \
      --storage-class STANDARD_IA \
      --region "$REGION" \
      --exclude "*.DS_Store" \
      --exclude "*.git*" \
      --exclude ".git/*"; then
        UPLOAD_SUCCESS=true
        echo ""
        echo "✅ Upload completed successfully!"
    else
        RETRY_COUNT=$((RETRY_COUNT + 1))
        if [ $RETRY_COUNT -lt $MAX_RETRIES ]; then
            echo "⚠️  Upload encountered errors. Retrying..."
        fi
    fi
done

# Final check
if [ "$UPLOAD_SUCCESS" = false ]; then
    echo ""
    echo "❌ Upload failed after $MAX_RETRIES attempts."
    echo ""
    echo "Troubleshooting tips:"
    echo "1. Check your internet connection stability"
    echo "2. Verify AWS credentials: aws sts get-caller-identity"
    echo "3. Check bucket permissions: aws s3 ls s3://$BUCKET/"
    echo "4. Try uploading a single file manually to test:"
    echo "   aws s3 cp \"$SOURCE_DIR/mate_R1.fq\" \"s3://$BUCKET/datasets/$DATASET_ID/mate_R1.fq\" --region $REGION"
    echo ""
    echo "You can also try using AWS S3 Transfer Acceleration if available:"
    echo "  aws configure set s3.use_accelerate_endpoint true"
    echo ""
    echo "To resume the upload, simply run the same command again:"
    echo "  $0 $DATASET_ID $SOURCE_DIR"
    exit 1
fi

echo ""
echo "=========================================="
echo "✅ Upload complete!"
echo "=========================================="
echo ""
echo "Dataset files are now available at:"
echo "  s3://$BUCKET/datasets/$DATASET_ID/"
echo ""
echo "To link this dataset to a session, use the API:"
echo "  POST /session/{session_id}/link-dataset"
echo "  {"
echo "    \"dataset_id\": \"$DATASET_ID\","
echo "    \"dataset_path\": \"datasets/$DATASET_ID\""
echo "  }"
echo ""
echo "Or use the Dataset Linker UI in the frontend."

