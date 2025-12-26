#!/bin/bash
# Fix rolled back stack by deleting and redeploying with desired_count=0

set -e

STACK_NAME="HelixAIStack"
REGION="us-west-1"

echo "🔧 Fixing Rolled Back Stack"
echo "============================"
echo ""

# Check current status
STATUS=$(aws cloudformation describe-stacks \
  --stack-name $STACK_NAME \
  --region $REGION \
  --query 'Stacks[0].StackStatus' \
  --output text 2>/dev/null || echo "NOT_FOUND")

echo "Current stack status: $STATUS"
echo ""

if [ "$STATUS" = "ROLLBACK_COMPLETE" ]; then
    echo "⚠️  Stack is in ROLLBACK_COMPLETE state"
    echo "   Need to delete it before redeploying"
    echo ""
    read -p "Delete stack and redeploy? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Cancelled"
        exit 1
    fi
    
    echo ""
    echo "🗑️  Deleting stack..."
    aws cloudformation delete-stack \
      --stack-name $STACK_NAME \
      --region $REGION
    
    echo "⏳ Waiting for deletion to complete..."
    aws cloudformation wait stack-delete-complete \
      --stack-name $STACK_NAME \
      --region $REGION
    
    echo "✅ Stack deleted"
    echo ""
fi

echo "📝 Modifying CDK code to set desired_count=0..."
echo "   (This allows stack to complete, then we can debug ECS service)"
echo ""
echo "Please modify infrastructure/helix_infrastructure/helix_stack.py:"
echo "  Change line 273 from: desired_count=1"
echo "  To:                   desired_count=0"
echo ""
read -p "Press Enter after making this change..."
echo ""

echo "🚀 Deploying stack..."
cd infrastructure
cdk deploy HelixAIStack --require-approval never

echo ""
echo "✅ Deployment complete!"
echo ""
echo "Next steps:"
echo "1. Check that stack deployed successfully"
echo "2. Investigate why ECS tasks were failing"
echo "3. Fix the issue"
echo "4. Update desired_count back to 1"
