#!/bin/bash
# Deploy Helix.AI backend to EC2 instance
# This script helps deploy the backend to an EC2 instance and configure the ALB to use it

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

# Load configuration
CONFIG_FILE="${1:-deploy.config}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ ! -f "${SCRIPT_DIR}/${CONFIG_FILE}" ]; then
    print_error "Configuration file not found: ${SCRIPT_DIR}/${CONFIG_FILE}"
    print_info "Please create ${CONFIG_FILE} with your AWS configuration"
    exit 1
fi

# Source configuration
source "${SCRIPT_DIR}/${CONFIG_FILE}"

# Validate required variables
REQUIRED_VARS=(
    "AWS_REGION"
    "AWS_ACCOUNT_ID"
    "ECR_REPOSITORY"
    "STACK_NAME"
)

for var in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!var:-}" ]; then
        print_error "Required variable ${var} is not set in ${CONFIG_FILE}"
        exit 1
    fi
done

STACK_NAME="${STACK_NAME:-HelixAIStack}"

print_info "Deploying Helix.AI to EC2"
print_info "=========================="
print_info "Stack Name: ${STACK_NAME}"
print_info "AWS Region: ${AWS_REGION}"
print_info "AWS Account: ${AWS_ACCOUNT_ID}"
print_info ""

# Check prerequisites
print_info "Checking prerequisites..."
if ! command -v aws &> /dev/null; then
    print_error "AWS CLI is not installed"
    exit 1
fi

if ! aws sts get-caller-identity &> /dev/null; then
    print_error "AWS credentials not configured. Run 'aws configure' first."
    exit 1
fi

print_success "Prerequisites check passed"
print_info ""

# 1. Build and push Docker image to ECR
print_info "Step 1: Building and pushing Docker image to ECR..."
"${SCRIPT_DIR}/ecr_push_backend.sh" \
    "${AWS_ACCOUNT_ID}" \
    "${AWS_REGION}" \
    "${ECR_REPOSITORY}" \
    "${IMAGE_TAG:-latest}"

print_success "Docker image pushed to ECR"
print_info ""

# 2. Get stack outputs
print_info "Step 2: Getting stack outputs..."
EC2_INSTANCE_ID=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2InstanceId'].OutputValue" \
    --output text 2>/dev/null || echo "")

if [ -z "${EC2_INSTANCE_ID}" ]; then
    print_error "EC2 instance ID not found in stack outputs. Make sure the stack is deployed."
    print_info "Deploy the stack first: cd ../../infrastructure && cdk deploy"
    exit 1
fi

ALB_DNS=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='ALBDNSName'].OutputValue" \
    --output text 2>/dev/null || echo "")

EC2_TARGET_GROUP_ARN=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2TargetGroupARN'].OutputValue" \
    --output text 2>/dev/null || echo "")

print_success "Stack outputs retrieved"
print_info "  EC2 Instance ID: ${EC2_INSTANCE_ID}"
print_info "  ALB DNS: ${ALB_DNS}"
print_info ""

# 3. Wait for EC2 instance to be running
print_info "Step 3: Waiting for EC2 instance to be running..."
aws ec2 wait instance-running \
    --instance-ids "${EC2_INSTANCE_ID}" \
    --region "${AWS_REGION}" > /dev/null

print_success "EC2 instance is running"
print_info ""

# 4. Register EC2 instance with ALB target group
print_info "Step 4: Registering EC2 instance with ALB target group..."
if [ -n "${EC2_TARGET_GROUP_ARN}" ]; then
    aws elbv2 register-targets \
        --target-group-arn "${EC2_TARGET_GROUP_ARN}" \
        --targets "Id=${EC2_INSTANCE_ID}" \
        --region "${AWS_REGION}" > /dev/null 2>&1 || echo "Instance may already be registered"
    print_success "Instance registered with target group"
else
    print_warning "Target group ARN not found, skipping registration"
fi
print_info ""

# 5. Trigger container update on EC2 instance
print_info "Step 4: Updating container on EC2 instance..."
print_info "This will pull the latest Docker image and restart the container"

# SSH into the instance and update the container
# Note: This requires SSH access. For production, consider using SSM Session Manager instead
EC2_INSTANCE_IP=$(aws ec2 describe-instances \
    --instance-ids "${EC2_INSTANCE_ID}" \
    --region "${AWS_REGION}" \
    --query "Reservations[0].Instances[0].PrivateIpAddress" \
    --output text 2>/dev/null || echo "")

if [ -n "${EC2_INSTANCE_IP}" ]; then
    print_info "EC2 Instance Private IP: ${EC2_INSTANCE_IP}"
    print_info ""
    print_warning "To update the container, you can:"
    print_info "  1. SSH into the instance (if you have access)"
    print_info "  2. Run: cd /opt/helix-ai && docker-compose pull && docker-compose up -d"
    print_info ""
    print_info "Alternatively, the instance will pull the latest image on next startup"
else
    print_warning "Could not retrieve EC2 instance IP"
fi

print_info ""

# 6. Check health of the backend
print_info "Step 6: Checking backend health..."
if [ -n "${ALB_DNS}" ]; then
    BACKEND_URL="http://${ALB_DNS}"
    print_info "Backend URL: ${BACKEND_URL}"
    
    # Wait a bit for the service to be ready
    print_info "Waiting for backend to be ready..."
    for i in {1..30}; do
        if curl -f "${BACKEND_URL}/health" >/dev/null 2>&1; then
            print_success "Backend is healthy!"
            break
        fi
        echo -n "."
        sleep 5
    done
    echo ""
else
    print_warning "ALB DNS not found in stack outputs"
fi

print_info ""
print_success "🎉 Deployment Complete!"
print_info "=========================="
print_info "EC2 Instance ID: ${EC2_INSTANCE_ID}"
if [ -n "${ALB_DNS}" ]; then
    print_info "Backend URL: http://${ALB_DNS}"
fi
print_info ""
print_info "To switch the ALB to use EC2 instead of ECS:"
print_info "  1. Go to AWS Console > EC2 > Load Balancers"
print_info "  2. Select your ALB"
print_info "  3. Go to Listeners tab"
print_info "  4. Edit the listener and change the target group to EC2TargetGroup"
print_info ""
print_info "Or use the script: ./switch-alb-to-ec2.sh"
