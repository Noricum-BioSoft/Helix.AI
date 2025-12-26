#!/usr/bin/env bash
set -euo pipefail

# Script to wait for CloudFormation stack to complete and get ALB URL

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

print_info "Waiting for CloudFormation stack: $STACK_NAME"
print_info "Region: $AWS_REGION"
echo ""

# Check current status
CURRENT_STATUS=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].StackStatus' \
    --output text 2>/dev/null || echo "NOT_FOUND")

if [ "${CURRENT_STATUS}" == "NOT_FOUND" ]; then
    print_error "Stack '${STACK_NAME}' not found!"
    exit 1
fi

if [ "${CURRENT_STATUS}" == "CREATE_COMPLETE" ] || [ "${CURRENT_STATUS}" == "UPDATE_COMPLETE" ]; then
    print_success "Stack is already complete (Status: ${CURRENT_STATUS})"
else
    print_info "Current status: ${CURRENT_STATUS}"
    print_info "Waiting for stack to complete (this may take 10-15 minutes)..."
    echo ""
    
    # Wait for stack to complete
    if aws cloudformation wait stack-create-complete \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" 2>&1; then
        print_success "Stack deployment completed successfully!"
    else
        print_error "Stack deployment failed or timed out"
        echo ""
        print_info "Checking final stack status..."
        FINAL_STATUS=$(aws cloudformation describe-stacks \
            --stack-name "${STACK_NAME}" \
            --region "${AWS_REGION}" \
            --query 'Stacks[0].StackStatus' \
            --output text 2>/dev/null)
        print_error "Final status: ${FINAL_STATUS}"
        
        if [ "${FINAL_STATUS}" == "CREATE_FAILED" ] || [ "${FINAL_STATUS}" == "ROLLBACK_COMPLETE" ]; then
            echo ""
            print_info "Recent stack events:"
            aws cloudformation describe-stack-events \
                --stack-name "${STACK_NAME}" \
                --region "${AWS_REGION}" \
                --max-items 5 \
                --query 'StackEvents[?ResourceStatus==`CREATE_FAILED` || ResourceStatus==`ROLLBACK_IN_PROGRESS`].[Timestamp,ResourceStatus,ResourceType,LogicalResourceId,ResourceStatusReason]' \
                --output table
        fi
        exit 1
    fi
fi

echo ""
print_info "Fetching ALB DNS name from stack outputs..."

# Get ALB URL
ALB_URL=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].Outputs[?OutputKey==`ALBDNSName`].OutputValue' \
    --output text 2>/dev/null)

if [ -z "${ALB_URL}" ] || [ "${ALB_URL}" == "None" ]; then
    print_error "ALB DNS name not found in stack outputs!"
    echo ""
    print_info "Available stack outputs:"
    aws cloudformation describe-stacks \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" \
        --query 'Stacks[0].Outputs[*].[OutputKey,OutputValue]' \
        --output table
    exit 1
fi

# Get other outputs
CLOUDFRONT_URL=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].Outputs[?OutputKey==`CloudFrontDistributionDomainName`].OutputValue' \
    --output text 2>/dev/null || echo "")

S3_BUCKET=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].Outputs[?OutputKey==`FrontendBucketName`].OutputValue' \
    --output text 2>/dev/null || echo "")

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Stack Deployment Complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
print_success "Backend ALB URL:"
echo "  http://${ALB_URL}"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Add this to your GitHub repository secrets:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  VITE_API_BASE_URL = http://${ALB_URL}"
echo ""

if [ -n "${CLOUDFRONT_URL}" ] && [ "${CLOUDFRONT_URL}" != "None" ]; then
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  Frontend URLs:"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo ""
    print_success "CloudFront URL:"
    echo "  https://${CLOUDFRONT_URL}"
    echo ""
    if [ -n "${S3_BUCKET}" ] && [ "${S3_BUCKET}" != "None" ]; then
        echo "  S3 Bucket: ${S3_BUCKET}"
    fi
    echo ""
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
print_info "Next steps:"
echo "  1. Add VITE_API_BASE_URL to GitHub Secrets"
echo "  2. Test backend health: curl http://${ALB_URL}/health"
echo "  3. Run GitHub Actions workflow to deploy frontend"
echo ""

# Test ALB health endpoint
print_info "Testing backend health endpoint..."
if curl -f -s -o /dev/null -w "%{http_code}" "http://${ALB_URL}/health" --max-time 10 | grep -q "200"; then
    print_success "Backend is healthy! ✓"
else
    print_warning "Backend health check failed (may still be starting)"
    echo "  ECS service may need a few more minutes to fully start"
fi

echo ""
print_info "To check ECS service status:"
echo "  cd scripts/aws && ./check-stack-status.sh ${STACK_NAME}"




