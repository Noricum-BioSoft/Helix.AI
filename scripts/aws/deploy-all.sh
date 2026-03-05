#!/bin/bash
# Complete deployment script for Helix.AI to AWS
# This script deploys everything: infrastructure, backend, and frontend

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

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Load configuration
CONFIG_FILE="${1:-deploy.config}"

if [ ! -f "${SCRIPT_DIR}/${CONFIG_FILE}" ]; then
    print_error "Configuration file not found: ${SCRIPT_DIR}/${CONFIG_FILE}"
    print_info "Creating example configuration file..."
    "${SCRIPT_DIR}/setup-deployment.sh" || exit 1
fi

source "${SCRIPT_DIR}/${CONFIG_FILE}"

print_info "=========================================="
print_info "  Helix.AI Complete AWS Deployment"
print_info "=========================================="
print_info ""
print_info "This script will:"
print_info "  1. Deploy infrastructure (CDK stack)"
print_info "  2. Build and push backend Docker image"
print_info "  3. Deploy frontend to S3/CloudFront"
print_info "  4. Configure ALB to use EC2 backend"
print_info ""
read -p "Continue? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    print_info "Deployment cancelled"
    exit 0
fi

print_info ""

# Step 1: Deploy Infrastructure
print_info "=========================================="
print_info "Step 1: Deploying Infrastructure (CDK)"
print_info "=========================================="

cd "${PROJECT_ROOT}/infrastructure"

if [ ! -d "node_modules" ]; then
    print_info "Installing CDK dependencies..."
    npm install
fi

print_info "Synthesizing CDK stack..."
cdk synth

print_info "Deploying CDK stack..."
cdk deploy --require-approval never

print_success "Infrastructure deployed!"
print_info ""

# Step 2: Deploy Backend to EC2
print_info "=========================================="
print_info "Step 2: Deploying Backend to EC2"
print_info "=========================================="

cd "${SCRIPT_DIR}"
"${SCRIPT_DIR}/deploy-to-ec2.sh" "${CONFIG_FILE}"

print_info ""

# Step 3: Switch ALB to use EC2
print_info "=========================================="
print_info "Step 3: Switching ALB to EC2 Backend"
print_info "=========================================="

read -p "Switch ALB to use EC2 backend? (y/N): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    "${SCRIPT_DIR}/switch-alb-to-ec2.sh" "${CONFIG_FILE}"
else
    print_info "Skipping ALB switch. You can do it later with: ./switch-alb-to-ec2.sh"
fi

print_info ""

# Step 4: Deploy Frontend
print_info "=========================================="
print_info "Step 4: Deploying Frontend"
print_info "=========================================="

if [ "${SKIP_FRONTEND:-false}" != "true" ]; then
    "${SCRIPT_DIR}/deploy.sh" "${CONFIG_FILE}"
else
    print_info "Skipping frontend deployment (SKIP_FRONTEND=true)"
fi

print_info ""
print_success "=========================================="
print_success "🎉 Complete Deployment Finished!"
print_success "=========================================="
print_info ""
print_info "Next steps:"
print_info "  1. Update your frontend to point to the backend URL"
print_info "  2. Test the application at the CloudFront URL"
print_info "  3. Monitor CloudWatch logs for any issues"
print_info ""
print_info "To view stack outputs:"
print_info "  aws cloudformation describe-stacks --stack-name ${STACK_NAME:-HelixAIStack} --query 'Stacks[0].Outputs'"
