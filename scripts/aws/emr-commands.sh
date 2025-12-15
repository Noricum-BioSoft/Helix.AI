#!/bin/bash
# Common EMR commands for Helix.AI FASTQ processing
# Source this file or use functions directly

# Configuration
REGION="${AWS_REGION:-us-east-1}"

# Function to get cluster status
emr_status() {
    local CLUSTER_ID=$1
    if [ -z "$CLUSTER_ID" ]; then
        echo "Usage: emr_status <cluster-id>"
        return 1
    fi
    aws emr describe-cluster \
        --cluster-id "$CLUSTER_ID" \
        --region "$REGION" \
        --query 'Cluster.Status.State' \
        --output text
}

# Function to list all EMR clusters
emr_list() {
    aws emr list-clusters \
        --region "$REGION" \
        --active \
        --query 'Clusters[*].[Id,Name,Status.State,CreationDateTime]' \
        --output table
}

# Function to get master node DNS
emr_master_dns() {
    local CLUSTER_ID=$1
    if [ -z "$CLUSTER_ID" ]; then
        echo "Usage: emr_master_dns <cluster-id>"
        return 1
    fi
    aws emr describe-cluster \
        --cluster-id "$CLUSTER_ID" \
        --region "$REGION" \
        --query 'Cluster.MasterPublicDnsName' \
        --output text
}

# Function to terminate cluster
emr_terminate() {
    local CLUSTER_ID=$1
    if [ -z "$CLUSTER_ID" ]; then
        echo "Usage: emr_terminate <cluster-id>"
        return 1
    fi
    echo "Terminating cluster $CLUSTER_ID..."
    aws emr terminate-clusters \
        --cluster-ids "$CLUSTER_ID" \
        --region "$REGION"
    echo "Termination initiated. Cluster will shut down in a few minutes."
}

# Function to submit Spark job
emr_submit_job() {
    local CLUSTER_ID=$1
    local SCRIPT_PATH=$2
    local ARGS=$3
    
    if [ -z "$CLUSTER_ID" ] || [ -z "$SCRIPT_PATH" ]; then
        echo "Usage: emr_submit_job <cluster-id> <s3-script-path> [args]"
        echo "Example: emr_submit_job j-123456789 s3://bucket/scripts/fastqc.py 's3://bucket/input/*.fq'"
        return 1
    fi
    
    aws emr add-steps \
        --cluster-id "$CLUSTER_ID" \
        --region "$REGION" \
        --steps Type=SPARK,Name="FASTQ Quality Analysis",ActionOnFailure=CONTINUE,Args="[--deploy-mode,cluster,--py-files,$SCRIPT_PATH,$ARGS]" \
        --query 'StepIds[0]' \
        --output text
}

# Function to check step status
emr_step_status() {
    local CLUSTER_ID=$1
    local STEP_ID=$2
    
    if [ -z "$CLUSTER_ID" ] || [ -z "$STEP_ID" ]; then
        echo "Usage: emr_step_status <cluster-id> <step-id>"
        return 1
    fi
    
    aws emr describe-step \
        --cluster-id "$CLUSTER_ID" \
        --step-id "$STEP_ID" \
        --region "$REGION" \
        --query 'Step.Status.State' \
        --output text
}

# Function to SSH into master node
emr_ssh() {
    local CLUSTER_ID=$1
    local KEY_FILE=$2
    
    if [ -z "$CLUSTER_ID" ]; then
        echo "Usage: emr_ssh <cluster-id> [key-file]"
        return 1
    fi
    
    local MASTER_DNS=$(emr_master_dns "$CLUSTER_ID")
    if [ -z "$MASTER_DNS" ]; then
        echo "Error: Could not get master DNS. Is cluster running?"
        return 1
    fi
    
    if [ -z "$KEY_FILE" ]; then
        echo "SSH command: ssh hadoop@$MASTER_DNS"
        echo "If you have a key file, use: ssh -i <key-file> hadoop@$MASTER_DNS"
    else
        ssh -i "$KEY_FILE" hadoop@"$MASTER_DNS"
    fi
}

# Function to create cluster with custom config
emr_create() {
    local CLUSTER_NAME=$1
    local INSTANCE_TYPE=${2:-m5.xlarge}
    local INSTANCE_COUNT=${3:-3}
    
    if [ -z "$CLUSTER_NAME" ]; then
        echo "Usage: emr_create <cluster-name> [instance-type] [instance-count]"
        echo "Example: emr_create helix-fastqc m5.xlarge 3"
        return 1
    fi
    
    echo "Creating EMR cluster: $CLUSTER_NAME"
    echo "Instance Type: $INSTANCE_TYPE"
    echo "Instance Count: $INSTANCE_COUNT"
    echo ""
    echo "Run ./scripts/aws/setup-emr-cluster.sh for full setup, or use AWS Console."
}

# Main menu
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    echo "EMR Helper Commands for Helix.AI"
    echo ""
    echo "Available functions:"
    echo "  emr_list              - List all active clusters"
    echo "  emr_status <id>       - Get cluster status"
    echo "  emr_master_dns <id>   - Get master node DNS"
    echo "  emr_terminate <id>    - Terminate cluster"
    echo "  emr_submit_job <id> <script> [args] - Submit Spark job"
    echo "  emr_step_status <id> <step> - Check job status"
    echo "  emr_ssh <id> [key]    - SSH into master node"
    echo ""
    echo "To use these functions, source this file:"
    echo "  source scripts/aws/emr-commands.sh"
fi







