#!/bin/bash
# Register EC2 instance with ALB target group
# This script should be run after the stack is deployed

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Load configuration
CONFIG_FILE="${1:-deploy.config}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ ! -f "${SCRIPT_DIR}/${CONFIG_FILE}" ]; then
    print_error "Configuration file not found: ${SCRIPT_DIR}/${CONFIG_FILE}"
    exit 1
fi

source "${SCRIPT_DIR}/${CONFIG_FILE}"

STACK_NAME="${STACK_NAME:-HelixAIStack}"

print_info "Registering EC2 instance with ALB target group..."
print_info "Stack: ${STACK_NAME}"
print_info "Region: ${AWS_REGION}"
print_info ""

# Get stack outputs
EC2_INSTANCE_ID=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2InstanceId'].OutputValue" \
    --output text 2>/dev/null || echo "")

EC2_TARGET_GROUP_ARN=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2TargetGroupARN'].OutputValue" \
    --output text 2>/dev/null || echo "")

if [ -z "${EC2_INSTANCE_ID}" ]; then
    print_error "EC2 instance ID not found in stack outputs"
    exit 1
fi

if [ -z "${EC2_TARGET_GROUP_ARN}" ]; then
    print_error "EC2 target group ARN not found in stack outputs"
    exit 1
fi

print_info "EC2 Instance ID: ${EC2_INSTANCE_ID}"
print_info "Target Group ARN: ${EC2_TARGET_GROUP_ARN}"
print_info ""

# Register the instance
print_info "Registering instance with target group..."
aws elbv2 register-targets \
    --target-group-arn "${EC2_TARGET_GROUP_ARN}" \
    --targets "Id=${EC2_INSTANCE_ID}" \
    --region "${AWS_REGION}" > /dev/null

print_success "Instance registered with target group!"

# Check target health
print_info "Checking target health..."
sleep 5
aws elbv2 describe-target-health \
    --target-group-arn "${EC2_TARGET_GROUP_ARN}" \
    --region "${AWS_REGION}" \
    --query "TargetHealthDescriptions[0].TargetHealth" \
    --output text

print_info ""
print_success "✅ EC2 instance registered with ALB target group!"
