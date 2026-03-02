#!/bin/bash
# Switch ALB listener to use EC2 target group instead of ECS

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

print_info "Switching ALB to use EC2 target group..."
print_info "Stack: ${STACK_NAME}"
print_info "Region: ${AWS_REGION}"
print_info ""

# Get stack outputs
EC2_TARGET_GROUP_ARN=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2TargetGroupARN'].OutputValue" \
    --output text 2>/dev/null || echo "")

if [ -z "${EC2_TARGET_GROUP_ARN}" ]; then
    print_error "EC2 target group ARN not found in stack outputs"
    exit 1
fi

# Get ALB ARN from stack resources
ALB_ARN=$(aws cloudformation describe-stack-resources \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "StackResources[?ResourceType=='AWS::ElasticLoadBalancingV2::LoadBalancer'].PhysicalResourceId" \
    --output text 2>/dev/null | head -1 || echo "")

if [ -z "${ALB_ARN}" ]; then
    print_error "ALB ARN not found in stack resources"
    exit 1
fi

# Get HTTP listener ARN (port 80)
HTTP_LISTENER_ARN=$(aws elbv2 describe-listeners \
    --load-balancer-arn "${ALB_ARN}" \
    --region "${AWS_REGION}" \
    --query "Listeners[?Port==\`80\`].ListenerArn" \
    --output text 2>/dev/null | head -1 || echo "")

if [ -z "${HTTP_LISTENER_ARN}" ]; then
    print_error "HTTP listener (port 80) not found"
    exit 1
fi

print_info "Found:"
print_info "  ALB ARN: ${ALB_ARN}"
print_info "  HTTP Listener ARN: ${HTTP_LISTENER_ARN}"
print_info "  EC2 Target Group ARN: ${EC2_TARGET_GROUP_ARN}"
print_info ""

# Update the listener to use EC2 target group
print_info "Updating HTTP listener to use EC2 target group..."
aws elbv2 modify-listener \
    --listener-arn "${HTTP_LISTENER_ARN}" \
    --default-actions "Type=forward,TargetGroupArn=${EC2_TARGET_GROUP_ARN}" \
    --region "${AWS_REGION}" > /dev/null

print_success "HTTP listener updated to use EC2 target group"

# Update HTTPS listener if it exists (port 443)
HTTPS_LISTENER_ARN=$(aws elbv2 describe-listeners \
    --load-balancer-arn "${ALB_ARN}" \
    --region "${AWS_REGION}" \
    --query "Listeners[?Port==\`443\`].ListenerArn" \
    --output text 2>/dev/null | head -1 || echo "")

if [ -n "${HTTPS_LISTENER_ARN}" ]; then
    print_info "Updating HTTPS listener to use EC2 target group..."
    aws elbv2 modify-listener \
        --listener-arn "${HTTPS_LISTENER_ARN}" \
        --default-actions "Type=forward,TargetGroupArn=${EC2_TARGET_GROUP_ARN}" \
        --region "${AWS_REGION}" > /dev/null
    
    print_success "HTTPS listener updated to use EC2 target group"
fi

print_info ""
print_success "✅ ALB switched to use EC2 backend!"
print_info ""
print_info "The ALB will now route traffic to the EC2 instance."
print_info "It may take a few minutes for the health checks to pass."
