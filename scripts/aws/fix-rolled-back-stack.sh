#!/usr/bin/env bash
set -euo pipefail

# Script to fix a rolled back CloudFormation stack
# This will delete the failed stack and prepare for redeployment

STACK_NAME="${1:-HelixAIStack}"
AWS_REGION="${AWS_REGION:-us-west-1}"

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

# Check current status
STATUS=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].StackStatus' \
    --output text 2>/dev/null || echo "NOT_FOUND")

if [ "${STATUS}" == "NOT_FOUND" ]; then
    print_error "Stack '${STACK_NAME}' not found!"
    exit 1
fi

echo ""
print_info "Current stack status: ${STATUS}"

if [ "${STATUS}" != "ROLLBACK_COMPLETE" ] && [ "${STATUS}" != "CREATE_FAILED" ] && [ "${STATUS}" != "ROLLBACK_IN_PROGRESS" ]; then
    print_warning "Stack is not in a failed/rollback state (Status: ${STATUS})"
    print_info "Nothing to fix. Stack status is: ${STATUS}"
    exit 0
fi

echo ""
print_warning "Stack is in failed state. You need to DELETE it before redeploying."
echo ""
print_info "Option 1: Delete the stack now (this script)"
print_info "Option 2: Delete manually: aws cloudformation delete-stack --stack-name ${STACK_NAME} --region ${AWS_REGION}"
echo ""

read -p "Delete the stack now? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    print_info "Cancelled. To delete manually, run:"
    echo "  aws cloudformation delete-stack --stack-name ${STACK_NAME} --region ${AWS_REGION}"
    exit 0
fi

echo ""
print_info "Deleting stack: ${STACK_NAME}..."

# Delete stack
if aws cloudformation delete-stack --stack-name "${STACK_NAME}" --region "${AWS_REGION}"; then
    print_success "Stack deletion initiated!"
    echo ""
    print_info "Waiting for stack deletion to complete..."
    echo "  (This may take 5-10 minutes)"
    echo ""
    
    if aws cloudformation wait stack-delete-complete \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" 2>&1; then
        print_success "Stack deleted successfully!"
        echo ""
        print_info "Next steps:"
        echo "  1. Ensure Docker image is correct in ECR"
        echo "  2. Redeploy the stack:"
        echo "     cd infrastructure"
        echo "     cdk deploy"
        echo ""
        print_warning "Note: Some resources (like ECR repository and S3 bucket) are set to RETAIN"
        print_warning "      so they will persist. The stack will recreate other resources."
    else
        print_error "Stack deletion failed or timed out"
        echo ""
        print_info "Check deletion status:"
        echo "  aws cloudformation describe-stacks --stack-name ${STACK_NAME} --region ${AWS_REGION}"
        exit 1
    fi
else
    print_error "Failed to initiate stack deletion"
    exit 1
fi
















