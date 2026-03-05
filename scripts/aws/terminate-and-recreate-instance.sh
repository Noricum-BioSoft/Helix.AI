#!/bin/bash
# Terminate and let CloudFormation recreate the EC2 instance with new IAM role and user data

set -euo pipefail

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
print_warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
print_error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Load config
CONFIG_FILE="${1:-deploy.config}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/${CONFIG_FILE}"

STACK_NAME="${STACK_NAME:-HelixAIStack}"

print_info "This will:"
print_info "  1. Update the CDK stack (fixes IAM role and user data)"
print_info "  2. Terminate the current EC2 instance"
print_info "  3. The stack will auto-recreate it with the fixes"
print_info ""

read -p "Continue? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    print_info "Cancelled"
    exit 0
fi

print_info ""

# Step 1: Update stack
print_info "Step 1: Updating CDK stack with fixes..."
cd "${SCRIPT_DIR}/../../infrastructure"
cdk deploy --require-approval never

print_success "Stack updated!"
print_info ""

# Step 2: Get instance ID
EC2_INSTANCE_ID=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2InstanceId'].OutputValue" \
    --output text 2>/dev/null || echo "")

if [ -z "${EC2_INSTANCE_ID}" ]; then
    print_error "EC2 instance ID not found"
    exit 1
fi

print_info "Step 2: Terminating instance ${EC2_INSTANCE_ID}..."
print_warning "This will destroy the current instance!"
print_info ""

read -p "Terminate instance ${EC2_INSTANCE_ID}? (y/N): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    print_info "Instance termination cancelled"
    exit 0
fi

aws ec2 terminate-instances --instance-ids "${EC2_INSTANCE_ID}" --region "${AWS_REGION}" > /dev/null

print_success "Instance termination requested"
print_info ""
print_info "The stack will automatically create a new instance with:"
print_info "  - Fixed user data (dnf instead of apt-get)"
print_info "  - SSM permissions (so SSM will work)"
print_info "  - S3 fallback in user data (will download from S3 if ECR fails)"
print_info ""
print_info "Monitor new instance creation with:"
print_info "  aws ec2 describe-instances --filters 'Name=tag:aws:cloudformation:stack-name,Values=${STACK_NAME}' --region ${AWS_REGION} --query 'Reservations[0].Instances[0].[InstanceId,State.Name]' --output text"
print_info ""
print_info "Once the new instance is running, the image should load automatically from S3!"
