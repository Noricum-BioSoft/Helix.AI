#!/bin/bash
# Retry script for ECR push with better error handling

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
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [ ! -f "${SCRIPT_DIR}/${CONFIG_FILE}" ]; then
    print_error "Configuration file not found: ${SCRIPT_DIR}/${CONFIG_FILE}"
    exit 1
fi

source "${SCRIPT_DIR}/${CONFIG_FILE}"

print_info "Retrying ECR push with retry logic..."
print_info "Repository: ${ECR_REPOSITORY}"
print_info "Region: ${AWS_REGION}"
print_info ""

# Function to push with retries
push_with_retry() {
    local max_attempts=3
    local attempt=1
    
    while [ $attempt -le $max_attempts ]; do
        print_info "Push attempt $attempt of $max_attempts..."
        
        if docker push "${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com/${ECR_REPOSITORY}:${IMAGE_TAG:-latest}" 2>&1 | tee /tmp/docker-push.log; then
            print_success "Push succeeded!"
            return 0
        else
            if [ $attempt -lt $max_attempts ]; then
                print_warning "Push failed. Waiting 5 seconds before retry..."
                sleep 5
            fi
            attempt=$((attempt + 1))
        fi
    done
    
    print_error "Push failed after $max_attempts attempts"
    return 1
}

# Check if image exists locally
IMAGE_NAME="${ECR_REPOSITORY}:${IMAGE_TAG:-latest}"
if ! docker images "${IMAGE_NAME}" | grep -q "${ECR_REPOSITORY}"; then
    print_warning "Image ${IMAGE_NAME} not found locally. Building first..."
    cd "${PROJECT_ROOT}"
    "${SCRIPT_DIR}/ecr_push_backend.sh" "${AWS_ACCOUNT_ID}" "${AWS_REGION}" "${ECR_REPOSITORY}" "${IMAGE_TAG:-latest}"
else
    # Just retry the push
    print_info "Image found locally. Retrying push..."
    
    # Login to ECR
    ECR_URI="${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com"
    print_info "Logging in to ECR..."
    aws ecr get-login-password --region "${AWS_REGION}" | docker login --username AWS --password-stdin "${ECR_URI}"
    
    # Tag if needed
    docker tag "${IMAGE_NAME}" "${ECR_URI}/${ECR_REPOSITORY}:${IMAGE_TAG:-latest}"
    
    # Push with retry
    push_with_retry
fi
