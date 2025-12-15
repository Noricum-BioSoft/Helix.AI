#!/bin/bash
# Submit FASTQ quality analysis job to EMR cluster

set -e

# Configuration
CLUSTER_ID="${EMR_CLUSTER_ID:-}"
REGION="${AWS_REGION:-us-east-1}"
S3_BUCKET="${S3_DATASET_BUCKET:-noricum-ngs-data}"
SCRIPT_BUCKET="${S3_SCRIPT_BUCKET:-noricum-ngs-data}"

# Input files (default to GRCh38.p12.MafHi dataset)
R1_PATH="${1:-s3://${S3_BUCKET}/datasets/GRCh38.p12.MafHi/mate_R1.fq}"
R2_PATH="${2:-s3://${S3_BUCKET}/datasets/GRCh38.p12.MafHi/mate_R2.fq}"
OUTPUT_PATH="${3:-s3://${S3_BUCKET}/fastqc-results/GRCh38.p12.MafHi/}"

if [ -z "$CLUSTER_ID" ]; then
    echo "Error: EMR_CLUSTER_ID not set"
    echo ""
    echo "Usage:"
    echo "  export EMR_CLUSTER_ID=j-XXXXXXXXXXXXX"
    echo "  $0 [r1_path] [r2_path] [output_path]"
    echo ""
    echo "Or:"
    echo "  EMR_CLUSTER_ID=j-XXXXXXXXXXXXX $0"
    echo ""
    echo "Example:"
    echo "  EMR_CLUSTER_ID=j-12QYDE51Q9LDP $0"
    exit 1
fi

echo "=========================================="
echo "Submitting FASTQ Quality Analysis Job"
echo "=========================================="
echo "Cluster ID: $CLUSTER_ID"
echo "R1 Input: $R1_PATH"
echo "R2 Input: $R2_PATH"
echo "Output: $OUTPUT_PATH"
echo ""

# Check if cluster is ready
CLUSTER_STATE=$(aws emr describe-cluster \
    --cluster-id "$CLUSTER_ID" \
    --region "$REGION" \
    --query 'Cluster.Status.State' \
    --output text 2>/dev/null)

if [ "$CLUSTER_STATE" != "WAITING" ] && [ "$CLUSTER_STATE" != "RUNNING" ]; then
    echo "Warning: Cluster $CLUSTER_ID is not ready. Current state: $CLUSTER_STATE"
    echo ""
    
    # If cluster is terminated, try to find an active cluster
    if [[ "$CLUSTER_STATE" == *"TERMINATED"* ]] || [ -z "$CLUSTER_STATE" ]; then
        echo "Attempting to find an active EMR cluster..."
        
        # Try to find an active cluster
        # --active flag already filters for WAITING or RUNNING clusters
        ACTIVE_CLUSTER=$(aws emr list-clusters \
            --region "$REGION" \
            --active \
            --query 'Clusters[0].[Id,Name,Status.State]' \
            --output text 2>/dev/null)
        
        if [ -n "$ACTIVE_CLUSTER" ]; then
            NEW_CLUSTER_ID=$(echo "$ACTIVE_CLUSTER" | awk '{print $1}')
            NEW_CLUSTER_NAME=$(echo "$ACTIVE_CLUSTER" | awk '{print $2}')
            NEW_CLUSTER_STATE=$(echo "$ACTIVE_CLUSTER" | awk '{print $3}')
            
            echo "✅ Found active cluster: $NEW_CLUSTER_ID ($NEW_CLUSTER_NAME) - State: $NEW_CLUSTER_STATE"
            echo "Using active cluster instead..."
            CLUSTER_ID="$NEW_CLUSTER_ID"
            CLUSTER_STATE="$NEW_CLUSTER_STATE"
        else
            echo "❌ No active EMR clusters found."
            echo ""
            echo "To create a new cluster, run:"
            echo "  ./scripts/aws/setup-emr-cluster.sh"
            echo ""
            echo "Or set EMR_CLUSTER_ID to an active cluster ID:"
            echo "  export EMR_CLUSTER_ID=j-XXXXXXXXXXXXX"
            echo ""
            echo "To list all clusters (including terminated):"
            echo "  aws emr list-clusters --region $REGION"
            echo ""
            echo "To list only active clusters:"
            echo "  aws emr list-clusters --region $REGION --active"
            exit 1
        fi
    else
        echo "Cluster must be in WAITING or RUNNING state to submit jobs."
        echo ""
        echo "Current cluster state: $CLUSTER_STATE"
        echo ""
        echo "Please wait for the cluster to be ready, or use a different cluster:"
        echo "  export EMR_CLUSTER_ID=j-XXXXXXXXXXXXX"
        exit 1
    fi
fi

# Upload PySpark script to S3 (use v2 which handles large files better)
# Get script directory (where this script is located)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="s3://${SCRIPT_BUCKET}/emr-scripts/fastqc_analysis_v2.py"
SCRIPT_FILE="${SCRIPT_DIR}/fastqc_analysis_v2.py"

# Fallback to v1 if v2 doesn't exist
if [ ! -f "$SCRIPT_FILE" ]; then
    SCRIPT_FILE="${SCRIPT_DIR}/fastqc_analysis.py"
    SCRIPT_PATH="s3://${SCRIPT_BUCKET}/emr-scripts/fastqc_analysis.py"
fi

echo "Uploading PySpark script to S3..."
aws s3 cp "$SCRIPT_FILE" "$SCRIPT_PATH" || {
    echo "Error: Failed to upload script. Make sure $SCRIPT_FILE exists."
    exit 1
}
echo "✅ Script uploaded to $SCRIPT_PATH"

# Create wrapper script that downloads and runs the PySpark script
WRAPPER_SCRIPT="/tmp/fastqc_wrapper_$$.sh"
LOG_FILE="/tmp/fastqc_wrapper_$$.log"
cat > "$WRAPPER_SCRIPT" <<WRAPPER_EOF
#!/bin/bash
set -e
set -x

# Redirect all output to log file and stdout
exec > >(tee -a $LOG_FILE)
exec 2>&1

echo "=========================================="
echo "FASTQ Analysis Wrapper Script"
echo "=========================================="
echo "Log file: $LOG_FILE"
echo ""

# Download PySpark script from S3
echo "Downloading PySpark script from S3..."
aws s3 cp $SCRIPT_PATH /tmp/fastqc_script.py

# Make it executable
chmod +x /tmp/fastqc_script.py

# Verify script was downloaded
if [ ! -f /tmp/fastqc_script.py ]; then
    echo "ERROR: Failed to download script from $SCRIPT_PATH"
    exit 1
fi

echo "Script downloaded successfully"
echo "File size: \$(wc -l < /tmp/fastqc_script.py) lines"
echo ""

# Run Spark job with error handling
echo "Starting Spark job..."
echo "Command: spark-submit --master yarn --deploy-mode cluster ..."
echo ""

# Run spark-submit and capture exit code properly
# Using client mode temporarily for easier debugging (driver runs on master node)
# Switch back to cluster mode for production
# Add EMR filesystem JARs to classpath to fix EmrFileSystem ClassNotFoundException
# EMR filesystem JARs can be in multiple locations depending on EMR version
echo "Searching for EMR filesystem JARs..."
EMRFS_JARS=""
# Try common locations
for dir in /usr/share/aws/emr/emrfs/lib/ /usr/share/aws/emr/emrfs/ /usr/lib/hadoop/lib/ /usr/share/aws/emr/; do
    if [ -d "$dir" ]; then
        FOUND=$(find "$dir" -name "*emrfs*.jar" -o -name "*aws*.jar" 2>/dev/null | head -5 | tr '\n' ',' | sed 's/,$//')
        if [ -n "$FOUND" ]; then
            EMRFS_JARS="$FOUND"
            echo "Found EMR filesystem JARs in $dir: $EMRFS_JARS"
            break
        fi
    fi
done

# Also try to find hadoop-aws JARs which provide S3A support
if [ -z "$EMRFS_JARS" ]; then
    HADOOP_AWS=$(find /usr/lib/hadoop* /usr/share/aws /opt -name "*hadoop-aws*.jar" 2>/dev/null | head -3 | tr '\n' ',' | sed 's/,$//')
    if [ -n "$HADOOP_AWS" ]; then
        EMRFS_JARS="$HADOOP_AWS"
        echo "Found Hadoop AWS JARs: $EMRFS_JARS"
    fi
fi

# Also try to find AWS SDK JARs needed for S3A
if [ -z "$EMRFS_JARS" ]; then
    AWS_SDK=$(find /usr/lib/hadoop* /usr/share/aws /opt -name "*aws-java-sdk*.jar" 2>/dev/null | head -2 | tr '\n' ',' | sed 's/,$//')
    if [ -n "$AWS_SDK" ]; then
        EMRFS_JARS="$AWS_SDK"
        echo "Found AWS SDK JARs: $EMRFS_JARS"
    fi
fi

if [ -n "$EMRFS_JARS" ]; then
    JARS_ARG="--jars $EMRFS_JARS"
else
    echo "Warning: EMR filesystem JARs not found in standard locations"
    echo "Will try to use default EMR S3 configuration"
    JARS_ARG=""
fi

spark-submit \\
  --master yarn \\
  --deploy-mode client \\
  --conf spark.pyspark.python=/usr/bin/python3 \\
  --conf spark.driver.memory=4g \\
  --conf spark.executor.memory=8g \\
  $JARS_ARG \\
  /tmp/fastqc_script.py \\
  --input-r1 "$R1_PATH" \\
  --input-r2 "$R2_PATH" \\
  --output "$OUTPUT_PATH" 2>&1 | tee -a $LOG_FILE

EXIT_CODE=\${PIPESTATUS[0]}

echo ""
if [ \$EXIT_CODE -eq 0 ]; then
    echo "✅ Spark job completed successfully"
else
    echo "❌ Spark job failed with exit code: \$EXIT_CODE"
    echo ""
    echo "Note: In cluster mode, driver logs are on the worker node."
    echo "To view detailed error logs, check YARN application logs:"
    echo "  yarn logs -applicationId <app-id>"
    echo ""
    echo "Or check the Spark UI at the tracking URL shown above."
fi

# Upload log file to S3 for debugging (always upload, success or failure)
echo "Uploading log file to S3..."
aws s3 cp $LOG_FILE s3://${S3_BUCKET}/emr-logs/$CLUSTER_ID/steps/\${EMR_STEP_ID}/wrapper.log || true

exit \$EXIT_CODE
WRAPPER_EOF

chmod +x "$WRAPPER_SCRIPT"

# Upload wrapper script to S3
WRAPPER_S3_PATH="s3://${SCRIPT_BUCKET}/emr-scripts/fastqc_wrapper_$$.sh"
echo "Uploading wrapper script to S3..."
aws s3 cp "$WRAPPER_SCRIPT" "$WRAPPER_S3_PATH" || {
    echo "Error: Failed to upload wrapper script"
    rm -f "$WRAPPER_SCRIPT"
    exit 1
}
rm -f "$WRAPPER_SCRIPT"
echo "✅ Wrapper script uploaded to $WRAPPER_S3_PATH"

# Submit Spark job
echo ""
echo "Submitting FASTQ analysis job to EMR cluster..."

# Create step JSON - download and execute wrapper script
STEP_JSON=$(cat <<EOF
[
  {
    "Type": "CUSTOM_JAR",
    "Name": "FASTQ Quality Analysis",
    "ActionOnFailure": "CONTINUE",
    "Jar": "command-runner.jar",
    "Args": [
      "bash", "-c",
      "aws s3 cp $WRAPPER_S3_PATH /tmp/fastqc_wrapper.sh && chmod +x /tmp/fastqc_wrapper.sh && /tmp/fastqc_wrapper.sh"
    ]
  }
]
EOF
)

# Use file-based input to avoid shell escaping issues
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
echo ""
echo "To monitor job status:"
echo "  aws emr describe-step --cluster-id $CLUSTER_ID --step-id $STEP_ID --region $REGION --query 'Step.Status.State' --output text"
echo ""
echo "Or use the helper:"
echo "  source scripts/aws/emr-commands.sh"
echo "  emr_step_status $CLUSTER_ID $STEP_ID"
echo ""
echo "To view logs:"
echo "  Check S3: s3://${S3_BUCKET}/emr-logs/$CLUSTER_ID/steps/$STEP_ID/"
echo ""
echo "Results will be available at: $OUTPUT_PATH"

