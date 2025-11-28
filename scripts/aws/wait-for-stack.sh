#!/usr/bin/env bash
set -euo pipefail

# Helper script to wait for CloudFormation stack to complete and show outputs

STACK_NAME="${1:-HelixAIStack}"
AWS_REGION="${AWS_REGION:-$(aws configure get region 2>/dev/null || echo 'us-west-1')}"
TIMEOUT="${2:-1800}"  # 30 minutes default timeout

print_info() {
    echo -e "\033[0;34m[INFO]\033[0m $1"
}

print_success() {
    echo -e "\033[0;32m[SUCCESS]\033[0m $1"
}

print_warning() {
    echo -e "\033[1;33m[WARNING]\033[0m $1"
}

print_error() {
    echo -e "\033[0;31m[ERROR]\033[0m $1"
}

print_info "Waiting for stack '${STACK_NAME}' to complete..."
print_info "Region: ${AWS_REGION}"
print_info "Timeout: ${TIMEOUT} seconds"
echo ""

START_TIME=$(date +%s)
ELAPSED=0
LAST_STATUS=""

while [ $ELAPSED -lt $TIMEOUT ]; do
    STATUS=$(aws cloudformation describe-stacks \
        --stack-name "${STACK_NAME}" \
        --region "${AWS_REGION}" \
        --query 'Stacks[0].StackStatus' \
        --output text 2>/dev/null || echo "NOT_FOUND")
    
    if [ "$STATUS" != "$LAST_STATUS" ]; then
        echo "[$(date +%H:%M:%S)] Stack Status: ${STATUS}"
        LAST_STATUS="$STATUS"
    fi
    
    case "$STATUS" in
        CREATE_COMPLETE|UPDATE_COMPLETE)
            print_success "Stack deployment completed successfully!"
            echo ""
            break
            ;;
        CREATE_FAILED|UPDATE_FAILED|ROLLBACK_COMPLETE|DELETE_COMPLETE)
            print_error "Stack deployment failed or was deleted!"
            echo ""
            echo "Recent stack events:"
            aws cloudformation describe-stack-events \
                --stack-name "${STACK_NAME}" \
                --region "${AWS_REGION}" \
                --max-items 10 \
                --query 'StackEvents[?ResourceStatus==`CREATE_FAILED` || ResourceStatus==`UPDATE_FAILED` || ResourceStatus==`DELETE_FAILED`].{Time:Timestamp,Resource:LogicalResourceId,Status:ResourceStatus,Reason:ResourceStatusReason}' \
                --output table
            exit 1
            ;;
        NOT_FOUND)
            print_error "Stack '${STACK_NAME}' not found!"
            exit 1
            ;;
        *)
            # Still in progress
            sleep 10
            ELAPSED=$(($(date +%s) - START_TIME))
            printf "\r[%02d:%02d] Waiting... " $((ELAPSED/60)) $((ELAPSED%60))
            ;;
    esac
done

if [ $ELAPSED -ge $TIMEOUT ]; then
    print_error "Timeout waiting for stack to complete!"
    exit 1
fi

echo ""
print_info "Stack Outputs:"
echo ""
aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'Stacks[0].Outputs' \
    --output table

echo ""
print_success "Stack is ready! You can now get the ALB URL using: ./get-alb-url.sh"

