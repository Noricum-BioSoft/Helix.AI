#!/usr/bin/env bash
set -euo pipefail

# Helper script to check CloudFormation stack status and wait for completion

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STACK_NAME="${1:-HelixAIStack}"
AWS_REGION="${AWS_REGION:-$(aws configure get region 2>/dev/null || echo 'us-west-1')}"

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_info "Checking stack status: ${STACK_NAME}"
print_info "Region: ${AWS_REGION}"
echo ""

# Get stack status
STATUS=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].StackStatus' \
    --output text 2>/dev/null || echo "NOT_FOUND")

if [ "${STATUS}" == "NOT_FOUND" ]; then
    print_error "Stack '${STACK_NAME}' not found in region ${AWS_REGION}"
    exit 1
fi

echo "Stack Status: ${STATUS}"
echo ""

# Check if stack is in progress
if [[ "${STATUS}" == *"_IN_PROGRESS" ]]; then
    print_warning "Stack deployment is still in progress!"
    echo ""
    print_info "Recent stack events:"
    aws cloudformation describe-stack-events \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" \
        --query 'StackEvents[0:5].[Timestamp,ResourceStatus,ResourceType,LogicalResourceId]' \
        --output table
    
    echo ""
    print_info "To monitor progress in real-time, run:"
    echo "  watch -n 5 'aws cloudformation describe-stack-events --stack-name ${STACK_NAME} --region ${AWS_REGION} --query \"StackEvents[0].[ResourceStatus,ResourceType,LogicalResourceId]\" --output table'"
    echo ""
    print_info "Or wait for completion with:"
    echo "  aws cloudformation wait stack-create-complete --stack-name ${STACK_NAME} --region ${AWS_REGION}"
    
elif [[ "${STATUS}" == *"COMPLETE"* ]] && [[ "${STATUS}" != *"ROLLBACK"* ]]; then
    print_success "Stack deployment completed successfully!"
    echo ""
    print_info "Stack outputs:"
    aws cloudformation describe-stacks \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" \
        --query 'Stacks[0].Outputs' \
        --output table
    
    echo ""
    print_info "Getting ALB URL..."
    ALB_URL=$(aws cloudformation describe-stacks \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" \
        --query 'Stacks[0].Outputs[?OutputKey==`ALBDNSName`].OutputValue' \
        --output text 2>/dev/null || echo "")
    
    if [ -n "${ALB_URL}" ] && [ "${ALB_URL}" != "None" ]; then
        echo ""
        print_success "Backend API URL: http://${ALB_URL}"
        echo ""
        echo "Add this to your deploy.config:"
        echo "  VITE_API_BASE_URL=http://${ALB_URL}"
    fi
else
    print_error "Stack is in an error state: ${STATUS}"
    echo ""
    print_info "Recent error events:"
    aws cloudformation describe-stack-events \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" \
        --query 'StackEvents[?ResourceStatus==`CREATE_FAILED` || ResourceStatus==`UPDATE_FAILED` || ResourceStatus==`DELETE_FAILED`].[Timestamp,ResourceStatus,ResourceType,LogicalResourceId,ResourceStatusReason]' \
        --output table
fi




