#!/bin/bash
# Submit a universal (non-FastQC) Helix.AI job to EMR using the universal runner.

set -e

CLUSTER_ID="${EMR_CLUSTER_ID:-}"
REGION="${AWS_REGION:-us-east-1}"
S3_BUCKET="${S3_DATASET_BUCKET:-noricum-ngs-data}"
SCRIPT_BUCKET="${S3_SCRIPT_BUCKET:-noricum-ngs-data}"

PAYLOAD_S3_URI="${1:-}"
OUTPUT_S3_PREFIX="${2:-s3://${S3_BUCKET}/emr-results/}"

if [ -z "$CLUSTER_ID" ]; then
  echo "Error: EMR_CLUSTER_ID not set"
  exit 1
fi

if [ -z "$PAYLOAD_S3_URI" ]; then
  echo "Usage: $0 s3://bucket/payload.json [s3://bucket/output-prefix/]"
  exit 1
fi

echo "=========================================="
echo "Submitting Universal EMR Job"
echo "=========================================="
echo "Cluster ID: $CLUSTER_ID"
echo "Payload: $PAYLOAD_S3_URI"
echo "Output: $OUTPUT_S3_PREFIX"
echo ""

# Upload runner script to S3
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUNNER_FILE="${SCRIPT_DIR}/universal_emr_runner.py"
RUNNER_S3_PATH="s3://${SCRIPT_BUCKET}/emr-scripts/universal_emr_runner.py"

echo "Uploading universal runner to S3..."
aws s3 cp "$RUNNER_FILE" "$RUNNER_S3_PATH"
echo "✅ Runner uploaded to $RUNNER_S3_PATH"

# Wrapper script
WRAPPER_SCRIPT="/tmp/universal_wrapper_$$.sh"
LOG_FILE="/tmp/universal_wrapper_$$.log"
cat > "$WRAPPER_SCRIPT" <<WRAPPER_EOF
#!/bin/bash
set -e
set -x

exec > >(tee -a $LOG_FILE)
exec 2>&1

echo "=========================================="
echo "Universal EMR Wrapper Script"
echo "=========================================="
echo "EMR_STEP_ID: \${EMR_STEP_ID}"
echo "Log file: $LOG_FILE"
echo ""

aws s3 cp $RUNNER_S3_PATH /tmp/universal_emr_runner.py
aws s3 cp $PAYLOAD_S3_URI /tmp/payload.json
chmod +x /tmp/universal_emr_runner.py

python3 /tmp/universal_emr_runner.py \
  --payload-file /tmp/payload.json \
  --output-s3 "$OUTPUT_S3_PREFIX\${EMR_STEP_ID}/" \
  --job-id "\${EMR_STEP_ID}"

EXIT_CODE=\${PIPESTATUS[0]}

echo "Uploading wrapper log..."
aws s3 cp $LOG_FILE s3://${S3_BUCKET}/emr-logs/$CLUSTER_ID/steps/\${EMR_STEP_ID}/wrapper.log || true

exit \$EXIT_CODE
WRAPPER_EOF

chmod +x "$WRAPPER_SCRIPT"

WRAPPER_S3_PATH="s3://${SCRIPT_BUCKET}/emr-scripts/universal_wrapper_$$.sh"
echo "Uploading wrapper script to S3..."
aws s3 cp "$WRAPPER_SCRIPT" "$WRAPPER_S3_PATH"
rm -f "$WRAPPER_SCRIPT"
echo "✅ Wrapper uploaded to $WRAPPER_S3_PATH"

echo "Submitting EMR step..."
STEP_JSON=$(cat <<EOF
[
  {
    "Type": "CUSTOM_JAR",
    "Name": "Helix Universal Job",
    "ActionOnFailure": "CONTINUE",
    "Jar": "command-runner.jar",
    "Args": [
      "bash", "-c",
      "aws s3 cp $WRAPPER_S3_PATH /tmp/universal_wrapper.sh && chmod +x /tmp/universal_wrapper.sh && /tmp/universal_wrapper.sh"
    ]
  }
]
EOF
)

TMP_STEP_FILE=$(mktemp)
echo "$STEP_JSON" > "$TMP_STEP_FILE"

STEP_ID=$(aws emr add-steps \
  --cluster-id "$CLUSTER_ID" \
  --region "$REGION" \
  --steps file://"$TMP_STEP_FILE" \
  --query 'StepIds[0]' \
  --output text)

rm -f "$TMP_STEP_FILE"

if [ -z "$STEP_ID" ] || [ "$STEP_ID" == "None" ]; then
  echo "❌ Failed to submit job"
  exit 1
fi

echo "✅ Job submitted successfully!"
echo ""
echo "Step ID: $STEP_ID"
echo "Results prefix: ${OUTPUT_S3_PREFIX}${STEP_ID}/"





