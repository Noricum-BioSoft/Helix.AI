#!/usr/bin/env bash
set -euo pipefail

# Script to delete retained resources (ECR repository and S3 bucket)
# These need to be deleted before redeploying the stack

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

AWS_REGION="${AWS_REGION:-us-west-1}"
PROJECT_NAME="helix-ai"
ACCOUNT_ID=$(aws sts get-caller-identity --query Account --output text)

print_warning "This will DELETE the following retained resources:"
echo "  • ECR Repository: ${PROJECT_NAME}-backend"
echo "  • S3 Bucket: ${PROJECT_NAME}-frontend-${ACCOUNT_ID}-${AWS_REGION}"
echo ""
print_warning "ALL Docker images in ECR will be DELETED!"
print_warning "ALL files in S3 bucket will be DELETED!"
echo ""
read -p "Are you sure you want to continue? (yes/NO): " -r
echo

if [[ ! $REPLY =~ ^[Yy][Ee][Ss]$ ]]; then
    print_info "Cancelled."
    exit 0
fi

# Delete ECR repository
print_info "Deleting ECR repository..."
ECR_REPO="${PROJECT_NAME}-backend"
if aws ecr describe-repositories --repository-names "$ECR_REPO" --region "$AWS_REGION" &>/dev/null; then
    print_warning "Deleting all images in ECR repository..."
    aws ecr list-images --repository-name "$ECR_REPO" --region "$AWS_REGION" --query 'imageIds[*]' --output json | \
        jq -r '.[] | "\(.imageDigest)"' | \
        while read digest; do
            aws ecr batch-delete-image --repository-name "$ECR_REPO" --image-ids imageDigest="$digest" --region "$AWS_REGION" &>/dev/null || true
        done
    
    print_info "Deleting ECR repository..."
    aws ecr delete-repository --repository-name "$ECR_REPO" --region "$AWS_REGION" --force
    print_success "ECR repository deleted"
else
    print_info "ECR repository doesn't exist, skipping"
fi

# Delete S3 bucket
print_info "Deleting S3 bucket..."
S3_BUCKET="${PROJECT_NAME}-frontend-${ACCOUNT_ID}-${AWS_REGION}"
if aws s3 ls "s3://${S3_BUCKET}" &>/dev/null 2>&1; then
    print_warning "Emptying S3 bucket..."
    aws s3 rm "s3://${S3_BUCKET}" --recursive --region "$AWS_REGION" || true
    
    print_info "Deleting S3 bucket..."
    aws s3 rb "s3://${S3_BUCKET}" --region "$AWS_REGION" --force
    print_success "S3 bucket deleted"
else
    print_info "S3 bucket doesn't exist, skipping"
fi

echo ""
print_success "All retained resources deleted!"
echo ""
print_info "You can now redeploy the stack:"
echo "  cd infrastructure"
echo "  cdk deploy"
















