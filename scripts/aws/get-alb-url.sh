#!/usr/bin/env bash
set -euo pipefail

# Helper script to get the ALB DNS name from CloudFormation stack outputs
# This is useful for finding the backend API URL to use in deploy.config

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

# Default stack name (can be overridden)
STACK_NAME="${1:-HelixAIStack}"
AWS_REGION="${AWS_REGION:-$(aws configure get region 2>/dev/null || echo 'us-east-1')}"

print_info "Fetching ALB DNS name from CloudFormation stack: ${STACK_NAME}"
print_info "AWS Region: ${AWS_REGION}"

# Check if stack exists
if ! aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].StackName' \
    --output text &>/dev/null; then
    print_warning "Stack '${STACK_NAME}' not found!"
    echo ""
    echo "Available stacks:"
    aws cloudformation list-stacks \
        --stack-status-filter CREATE_COMPLETE UPDATE_COMPLETE \
        --query 'StackSummaries[?contains(StackName, `Helix`) || contains(StackName, `helix`)].StackName' \
        --output table \
        --region "${AWS_REGION}" 2>/dev/null || echo "No stacks found"
    exit 1
fi

# Check stack status first
STACK_STATUS=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].StackStatus' \
    --output text 2>/dev/null || echo "UNKNOWN")

if [[ "${STACK_STATUS}" == *"_IN_PROGRESS" ]]; then
    print_warning "Stack deployment is still in progress (Status: ${STACK_STATUS})"
    echo ""
    echo "Stack outputs are only available after deployment completes."
    echo "This typically takes 10-15 minutes for the initial deployment."
    echo ""
    print_info "To check stack status, run:"
    echo "  ./check-stack-status.sh ${STACK_NAME}"
    echo ""
    print_info "Or wait for completion with:"
    echo "  aws cloudformation wait stack-create-complete --stack-name ${STACK_NAME} --region ${AWS_REGION}"
    echo ""
    exit 0
elif [[ "${STACK_STATUS}" == *"ROLLBACK"* ]] || [[ "${STACK_STATUS}" == *"FAILED"* ]]; then
    print_error "Stack is in an error state: ${STACK_STATUS}"
    echo ""
    echo "Check stack events for details:"
    echo "  aws cloudformation describe-stack-events --stack-name ${STACK_NAME} --region ${AWS_REGION} --query 'StackEvents[?ResourceStatus==\`CREATE_FAILED\`]' --output table"
    exit 1
fi

# Get ALB DNS name from stack outputs
ALB_URL=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].Outputs[?OutputKey==`ALBDNSName`].OutputValue' \
    --output text 2>/dev/null)

if [ -z "${ALB_URL}" ] || [ "${ALB_URL}" == "None" ]; then
    print_warning "ALB DNS name not found in stack outputs!"
    echo ""
    echo "All stack outputs:"
    aws cloudformation describe-stacks \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" \
        --query 'Stacks[0].Outputs' \
        --output table
    exit 1
fi

# Construct full URL
BACKEND_URL="http://${ALB_URL}"

echo ""
print_success "Backend API URL found!"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Backend API URL: ${BACKEND_URL}"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Test if backend is accessible
print_info "Testing backend health endpoint..."
if curl -s --max-time 5 "${BACKEND_URL}/health" > /dev/null 2>&1; then
    print_success "Backend is accessible and responding!"
    echo ""
    echo "Health check response:"
    curl -s "${BACKEND_URL}/health" | python3 -m json.tool 2>/dev/null || curl -s "${BACKEND_URL}/health"
else
    print_warning "Backend is not responding yet. This is normal if:"
    echo "  - ECS service is still starting up (wait 2-5 minutes)"
    echo "  - No tasks are running yet"
    echo "  - Health checks haven't passed"
    echo ""
    echo "You can still use this URL in deploy.config, but the backend"
    echo "needs to be healthy before the frontend can connect to it."
fi

echo ""
print_info "To use this URL in deploy.config:"
echo ""
echo "  VITE_API_BASE_URL=${BACKEND_URL}"
echo ""
echo "Or copy this line:"
echo ""
echo "${BACKEND_URL}" | pbcopy 2>/dev/null && echo "✓ Copied to clipboard!" || echo "${BACKEND_URL}"

