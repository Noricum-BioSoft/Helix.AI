#!/usr/bin/env bash
set -euo pipefail

# Master Deployment Script for Helix.AI to AWS
# This script automates the complete deployment process including:
# - Backend Docker image build and push to ECR
# - Frontend build and sync to S3
# - CloudFront invalidation
# - ECS service update (if configured)

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
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [ ! -f "${SCRIPT_DIR}/${CONFIG_FILE}" ]; then
    print_error "Configuration file not found: ${SCRIPT_DIR}/${CONFIG_FILE}"
    print_info "Creating example configuration file..."
    cat > "${SCRIPT_DIR}/deploy.config.example" << 'EOF'
# Helix.AI AWS Deployment Configuration
# Copy this file to deploy.config and update with your values

# AWS Configuration
AWS_REGION=us-east-1
AWS_ACCOUNT_ID=your-account-id-here

# ECR Configuration
ECR_REPOSITORY=helix-backend
IMAGE_TAG=latest

# S3 Configuration
S3_BUCKET=helix-frontend-bucket

# CloudFront Configuration (optional)
CLOUDFRONT_DISTRIBUTION_ID=

# ECS Configuration (optional - for automatic service updates)
ECS_CLUSTER_NAME=
ECS_SERVICE_NAME=
ECS_TASK_DEFINITION_FAMILY=

# Frontend Configuration
VITE_API_BASE_URL=https://your-backend-alb-url.region.elb.amazonaws.com

# Deployment Options
SKIP_BACKEND=false
SKIP_FRONTEND=false
SKIP_CLOUDFRONT_INVALIDATION=false
SKIP_ECS_UPDATE=false
EOF
    print_warning "Please create ${CONFIG_FILE} with your AWS configuration"
    exit 1
fi

# Source configuration
source "${SCRIPT_DIR}/${CONFIG_FILE}"

# Validate required variables
REQUIRED_VARS=(
    "AWS_REGION"
    "AWS_ACCOUNT_ID"
    "ECR_REPOSITORY"
    "S3_BUCKET"
)

for var in "${REQUIRED_VARS[@]}"; do
    if [ -z "${!var:-}" ]; then
        print_error "Required variable ${var} is not set in ${CONFIG_FILE}"
        exit 1
    fi
done

# Build ECR URI
ECR_URI="${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com"

print_info "Starting Helix.AI AWS Deployment"
print_info "=================================="
print_info "AWS Region: ${AWS_REGION}"
print_info "AWS Account: ${AWS_ACCOUNT_ID}"
print_info "ECR Repository: ${ECR_REPOSITORY}"
print_info "S3 Bucket: ${S3_BUCKET}"
print_info "Image Tag: ${IMAGE_TAG:-latest}"

# Check prerequisites
print_info "Checking prerequisites..."

if ! command -v aws &> /dev/null; then
    print_error "AWS CLI is not installed"
    exit 1
fi

if ! command -v docker &> /dev/null; then
    print_error "Docker is not installed"
    exit 1
fi

if ! command -v node &> /dev/null; then
    print_error "Node.js is not installed"
    exit 1
fi

if ! aws sts get-caller-identity &> /dev/null; then
    print_error "AWS credentials not configured. Run 'aws configure' first."
    exit 1
fi

print_success "Prerequisites check passed"

# Deploy Backend
if [ "${SKIP_BACKEND:-false}" != "true" ]; then
    print_info "=================================="
    print_info "Deploying Backend to ECR"
    print_info "=================================="
    
    "${SCRIPT_DIR}/ecr_push_backend.sh" \
        "${AWS_ACCOUNT_ID}" \
        "${AWS_REGION}" \
        "${ECR_REPOSITORY}" \
        "${IMAGE_TAG:-latest}"
    
    print_success "Backend deployment completed"
    
    # Update ECS service if configured
    if [ "${SKIP_ECS_UPDATE:-false}" != "true" ] && [ -n "${ECS_CLUSTER_NAME:-}" ] && [ -n "${ECS_SERVICE_NAME:-}" ]; then
        print_info "=================================="
        print_info "Updating ECS Service"
        print_info "=================================="
        
        # Get current task definition
        if [ -n "${ECS_TASK_DEFINITION_FAMILY:-}" ]; then
            TASK_DEF=$(aws ecs describe-task-definition \
                --task-definition "${ECS_TASK_DEFINITION_FAMILY}" \
                --query 'taskDefinition' \
                --region "${AWS_REGION}" \
                --output json)
            
            # Create new task definition with updated image
            NEW_TASK_DEF=$(echo "${TASK_DEF}" | jq --arg IMAGE "${ECR_URI}/${ECR_REPOSITORY}:${IMAGE_TAG:-latest}" \
                '.containerDefinitions[0].image = $IMAGE | del(.taskDefinitionArn) | del(.revision) | del(.status) | del(.requiresAttributes) | del(.compatibilities) | del(.registeredAt) | del(.registeredBy)')
            
            # Register new task definition
            NEW_TASK_DEF_ARN=$(aws ecs register-task-definition \
                --cli-input-json "${NEW_TASK_DEF}" \
                --region "${AWS_REGION}" \
                --query 'taskDefinition.taskDefinitionArn' \
                --output text)
            
            print_success "New task definition registered: ${NEW_TASK_DEF_ARN}"
            
            # Update ECS service
            aws ecs update-service \
                --cluster "${ECS_CLUSTER_NAME}" \
                --service "${ECS_SERVICE_NAME}" \
                --task-definition "${NEW_TASK_DEF_ARN}" \
                --force-new-deployment \
                --region "${AWS_REGION}" > /dev/null
            
            print_success "ECS service update initiated. New tasks will be deployed shortly."
        else
            print_warning "ECS_TASK_DEFINITION_FAMILY not set, skipping ECS update"
        fi
    fi
else
    print_info "Skipping backend deployment (SKIP_BACKEND=true)"
fi

# Deploy Frontend
if [ "${SKIP_FRONTEND:-false}" != "true" ]; then
    print_info "=================================="
    print_info "Deploying Frontend to S3"
    print_info "=================================="
    
    # Set VITE_API_BASE_URL if provided
    if [ -n "${VITE_API_BASE_URL:-}" ]; then
        export VITE_API_BASE_URL="${VITE_API_BASE_URL}"
        print_info "Using VITE_API_BASE_URL=${VITE_API_BASE_URL}"
    fi
    
    "${SCRIPT_DIR}/s3_sync_frontend.sh" \
        "${AWS_REGION}" \
        "${S3_BUCKET}" \
        "${CLOUDFRONT_DISTRIBUTION_ID:-}"
    
    print_success "Frontend deployment completed"
else
    print_info "Skipping frontend deployment (SKIP_FRONTEND=true)"
fi

print_info "=================================="
print_success "🎉 Deployment Complete!"
print_info "=================================="
print_info "Backend Image: ${ECR_URI}/${ECR_REPOSITORY}:${IMAGE_TAG:-latest}"
print_info "Frontend URL: https://${S3_BUCKET}.s3-website-${AWS_REGION}.amazonaws.com"
if [ -n "${CLOUDFRONT_DISTRIBUTION_ID:-}" ]; then
    CLOUDFRONT_DOMAIN=$(aws cloudfront get-distribution \
        --id "${CLOUDFRONT_DISTRIBUTION_ID}" \
        --query 'Distribution.DomainName' \
        --region "${AWS_REGION}" \
        --output text 2>/dev/null || echo "Check CloudFront console")
    print_info "CloudFront URL: https://${CLOUDFRONT_DOMAIN}"
fi

