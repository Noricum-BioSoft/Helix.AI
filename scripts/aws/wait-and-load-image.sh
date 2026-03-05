#!/bin/bash
# Wait for SSM to become available, then load image from S3

set -euo pipefail

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }

# Load config
CONFIG_FILE="${1:-deploy.config}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/${CONFIG_FILE}"

STACK_NAME="${STACK_NAME:-HelixAIStack}"
S3_BUCKET="helix-images-${AWS_ACCOUNT_ID}-${AWS_REGION}"

# Get EC2 instance ID
EC2_INSTANCE_ID=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2InstanceId'].OutputValue" \
    --output text 2>/dev/null || echo "")

if [ -z "${EC2_INSTANCE_ID}" ]; then
    echo "❌ EC2 instance ID not found"
    exit 1
fi

print_info "Waiting for SSM to become available for instance: ${EC2_INSTANCE_ID}"
print_info "This may take 5-10 minutes after instance launch..."
print_info ""

MAX_WAIT=600  # 10 minutes
WAITED=0
CHECK_INTERVAL=30

while [ $WAITED -lt $MAX_WAIT ]; do
    # Try to send a simple SSM command
    if aws ssm send-command \
        --instance-ids "${EC2_INSTANCE_ID}" \
        --region "${AWS_REGION}" \
        --document-name "AWS-RunShellScript" \
        --parameters 'commands=["echo SSM is working"]' \
        --query "Command.CommandId" \
        --output text >/dev/null 2>&1; then
        
        print_success "SSM is now available!"
        print_info ""
        
        # Now load the image
        print_info "Loading Docker image from S3..."
        COMMAND_ID=$(aws ssm send-command \
            --instance-ids "${EC2_INSTANCE_ID}" \
            --region "${AWS_REGION}" \
            --document-name "AWS-RunShellScript" \
            --parameters "commands=[
                'cd /opt/helix-ai',
                'aws s3 cp s3://${S3_BUCKET}/helix-backend-latest.tar.gz /tmp/helix-backend.tar.gz --region ${AWS_REGION}',
                'gunzip -f /tmp/helix-backend.tar.gz',
                'docker load -i /tmp/helix-backend.tar',
                'docker tag helix-backend:latest helix-ai-backend:latest || true',
                'docker stop helix-ai-backend || true',
                'docker rm helix-ai-backend || true',
                'docker-compose up -d',
                'rm -f /tmp/helix-backend.tar*',
                'echo Image loaded and container started!',
                'docker ps | grep helix-ai-backend'
            ]" \
            --query "Command.CommandId" \
            --output text)
        
        print_success "Load command sent! Command ID: ${COMMAND_ID}"
        print_info ""
        print_info "Check status with:"
        print_info "  aws ssm get-command-invocation --command-id ${COMMAND_ID} --instance-id ${EC2_INSTANCE_ID} --region ${AWS_REGION}"
        print_info ""
        print_info "Or wait 30 seconds and check:"
        sleep 30
        aws ssm get-command-invocation \
            --command-id "${COMMAND_ID}" \
            --instance-id "${EC2_INSTANCE_ID}" \
            --region "${AWS_REGION}" \
            --query '[Status,StandardOutputContent,StandardErrorContent]' \
            --output text
        
        exit 0
    fi
    
    print_warning "SSM not ready yet... (waited ${WAITED}s / ${MAX_WAIT}s)"
    sleep $CHECK_INTERVAL
    WAITED=$((WAITED + CHECK_INTERVAL))
done

print_warning "SSM did not become available within ${MAX_WAIT} seconds"
print_info ""
print_info "Alternative options:"
print_info "  1. Restart the instance to trigger updated user data (with S3 fallback)"
print_info "  2. Check IAM role permissions for SSM"
print_info "  3. Wait longer and try again"
exit 1
