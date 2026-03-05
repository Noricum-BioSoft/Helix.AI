#!/usr/bin/env bash
set -euo pipefail

# Script to restart the backend ECS service with updated code
# This will rebuild the Docker image and force a new deployment

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

# Configuration
STACK_NAME="${1:-HelixAIStack}"
AWS_REGION="${AWS_REGION:-us-west-1}"
IMAGE_TAG="${IMAGE_TAG:-latest}"

print_info "Restarting backend ECS service"
print_info "Stack: ${STACK_NAME}"
print_info "Region: ${AWS_REGION}"
echo ""

# Get stack outputs
print_info "Fetching ECS configuration from CloudFormation stack..."

CLUSTER_NAME=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].Outputs[?OutputKey==`ECSClusterName`].OutputValue' \
    --output text 2>/dev/null || echo "")

SERVICE_NAME=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].Outputs[?OutputKey==`ECSServiceName`].OutputValue' \
    --output text 2>/dev/null || echo "")

ECR_REPO=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].Outputs[?OutputKey==`ECRRepositoryURI`].OutputValue' \
    --output text 2>/dev/null || echo "")

AWS_ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text --region "${AWS_REGION}")

if [ -z "${CLUSTER_NAME}" ] || [ "${CLUSTER_NAME}" == "None" ]; then
    print_error "Could not find ECS cluster name in stack outputs"
    exit 1
fi

if [ -z "${SERVICE_NAME}" ] || [ "${SERVICE_NAME}" == "None" ]; then
    print_error "Could not find ECS service name in stack outputs"
    exit 1
fi

if [ -z "${ECR_REPO}" ] || [ "${ECR_REPO}" == "None" ]; then
    # Try to construct it from account ID
    ECR_REPO="${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com/helix-ai-backend"
    print_warning "ECR repository not found in outputs, using: ${ECR_REPO}"
fi

print_success "Found ECS configuration:"
echo "  Cluster: ${CLUSTER_NAME}"
echo "  Service: ${SERVICE_NAME}"
echo "  ECR Repo: ${ECR_REPO}"
echo ""

# Check if we're in the project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [ ! -f "${PROJECT_ROOT}/backend/Dockerfile" ]; then
    print_error "Could not find backend/Dockerfile. Are you in the project root?"
    exit 1
fi

# Step 1: Build and push Docker image
print_info "=================================="
print_info "Step 1: Building and pushing Docker image"
print_info "=================================="

cd "${PROJECT_ROOT}"

# Login to ECR
print_info "Logging in to ECR..."
aws ecr get-login-password --region "${AWS_REGION}" | \
    docker login --username AWS --password-stdin "${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com"

# Build for linux/amd64 (required for ECS Fargate)
print_info "Building Docker image for linux/amd64 platform..."
IMAGE_NAME="${ECR_REPO}:${IMAGE_TAG}"
docker build --platform linux/amd64 -f backend/Dockerfile -t "${IMAGE_NAME}" .

# Push to ECR with retry logic
print_info "Pushing to ECR (this may take a few minutes)..."
MAX_RETRIES=3
RETRY_COUNT=0
PUSH_SUCCESS=false

while [ $RETRY_COUNT -lt $MAX_RETRIES ]; do
    if docker push "${IMAGE_NAME}"; then
        PUSH_SUCCESS=true
        break
    else
        RETRY_COUNT=$((RETRY_COUNT + 1))
        if [ $RETRY_COUNT -lt $MAX_RETRIES ]; then
            print_warning "Push failed (attempt ${RETRY_COUNT}/${MAX_RETRIES}). Retrying in 5 seconds..."
            sleep 5
            # Re-login to ECR in case token expired
            print_info "Refreshing ECR login..."
            aws ecr get-login-password --region "${AWS_REGION}" | \
                docker login --username AWS --password-stdin "${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com"
        fi
    fi
done

if [ "$PUSH_SUCCESS" = false ]; then
    print_error "Failed to push Docker image after ${MAX_RETRIES} attempts"
    print_info "You can try pushing manually with:"
    echo "  docker push ${IMAGE_NAME}"
    exit 1
fi

print_success "Docker image pushed: ${IMAGE_NAME}"
echo ""

# Step 2: Force new deployment
print_info "=================================="
print_info "Step 2: Forcing new ECS deployment"
print_info "=================================="

print_info "Updating ECS service to force new deployment..."
aws ecs update-service \
    --cluster "${CLUSTER_NAME}" \
    --service "${SERVICE_NAME}" \
    --force-new-deployment \
    --region "${AWS_REGION}" > /dev/null

print_success "ECS service update initiated!"
echo ""
print_info "The service is now deploying new tasks with the updated code."
print_info "This may take a few minutes. You can monitor the deployment with:"
echo ""
echo "  aws ecs describe-services \\"
echo "    --cluster ${CLUSTER_NAME} \\"
echo "    --services ${SERVICE_NAME} \\"
echo "    --region ${AWS_REGION} \\"
echo "    --query 'services[0].deployments[*].[status,desiredCount,runningCount]' \\"
echo "    --output table"
echo ""

print_info "Or check the deployment status:"
echo "  ./scripts/aws/check-deployment-status.sh"
echo ""

print_success "Backend restart initiated! ✓"

