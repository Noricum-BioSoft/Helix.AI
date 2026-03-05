#!/bin/bash
# Wait for any in-progress CloudFormation updates, then deploy
# This prevents "stack is currently being updated" errors

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

# Set required environment variables
export AWS_PROFILE=copilot-admin
export AWS_REGION=us-west-1

STACK_NAME="helix-ai-production-backend"

echo "⏳ Checking CloudFormation stack status..."

# Check if stack is updating
STATUS=$(aws cloudformation describe-stacks \
    --stack-name "$STACK_NAME" \
    --region "$AWS_REGION" \
    --query 'Stacks[0].StackStatus' \
    --output text 2>/dev/null || echo "NOT_FOUND")

if [ "$STATUS" = "NOT_FOUND" ]; then
    echo "✅ Stack not found, proceeding with deployment..."
elif [ "$STATUS" = "UPDATE_IN_PROGRESS" ] || [ "$STATUS" = "CREATE_IN_PROGRESS" ]; then
    echo "⏳ Stack is currently updating (status: $STATUS)"
    echo "   Waiting for update to complete..."
    
    # Wait for stack to finish updating
    aws cloudformation wait stack-update-complete \
        --stack-name "$STACK_NAME" \
        --region "$AWS_REGION" 2>/dev/null || \
    aws cloudformation wait stack-create-complete \
        --stack-name "$STACK_NAME" \
        --region "$AWS_REGION" 2>/dev/null || true
    
    # Check final status
    FINAL_STATUS=$(aws cloudformation describe-stacks \
        --stack-name "$STACK_NAME" \
        --region "$AWS_REGION" \
        --query 'Stacks[0].StackStatus' \
        --output text)
    
    if [ "$FINAL_STATUS" = "UPDATE_COMPLETE" ] || [ "$FINAL_STATUS" = "CREATE_COMPLETE" ]; then
        echo "✅ Stack update completed successfully"
    elif [ "$FINAL_STATUS" = "UPDATE_ROLLBACK_COMPLETE" ] || [ "$FINAL_STATUS" = "ROLLBACK_COMPLETE" ]; then
        echo "⚠️  Stack update rolled back. Check CloudFormation console for details."
        echo "   You may need to fix issues before redeploying."
        exit 1
    else
        echo "⚠️  Stack in unexpected state: $FINAL_STATUS"
        echo "   Proceeding with deployment anyway..."
    fi
elif [ "$STATUS" = "UPDATE_ROLLBACK_IN_PROGRESS" ] || [ "$STATUS" = "ROLLBACK_IN_PROGRESS" ]; then
    echo "⏳ Stack is rolling back (status: $STATUS)"
    echo "   Waiting for rollback to complete..."
    
    # Wait for rollback to complete
    aws cloudformation wait stack-rollback-complete \
        --stack-name "$STACK_NAME" \
        --region "$AWS_REGION" 2>/dev/null || \
    aws cloudformation wait stack-update-complete \
        --stack-name "$STACK_NAME" \
        --region "$AWS_REGION" 2>/dev/null || true
    
    FINAL_STATUS=$(aws cloudformation describe-stacks \
        --stack-name "$STACK_NAME" \
        --region "$AWS_REGION" \
        --query 'Stacks[0].StackStatus' \
        --output text)
    
    if [ "$FINAL_STATUS" = "UPDATE_ROLLBACK_COMPLETE" ] || [ "$FINAL_STATUS" = "ROLLBACK_COMPLETE" ]; then
        echo "⚠️  Stack rollback completed. Previous deployment failed."
        echo "   Please investigate the cause before redeploying."
        echo "   Check CloudFormation events and ECS task logs for details."
        exit 1
    else
        echo "⚠️  Rollback finished with status: $FINAL_STATUS"
        exit 1
    fi
elif [ "$STATUS" = "UPDATE_ROLLBACK_COMPLETE" ] || [ "$STATUS" = "ROLLBACK_COMPLETE" ]; then
    echo "⚠️  Stack is in rollback state (status: $STATUS)"
    echo "   Previous deployment failed. Please investigate before redeploying."
    exit 1
else
    echo "✅ Stack is ready (status: $STATUS)"
fi

echo ""
echo "🚀 Deploying backend service to production environment..."

# Deploy
copilot svc deploy --name backend --env production

echo ""
echo "✅ Deployment initiated. Monitor with:"
echo "   copilot svc status --name backend --env production"
echo "   copilot svc logs --name backend --env production --follow"

