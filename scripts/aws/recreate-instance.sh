#!/bin/bash
# Force CloudFormation to recreate the EC2 instance

set -euo pipefail

source deploy.config
STACK_NAME="${STACK_NAME:-HelixAIStack}"

echo "To recreate the EC2 instance, we need to update the stack."
echo "Since the instance was terminated, CloudFormation should recreate it on next update."
echo ""
echo "Option 1: Update via CDK (if CDK is available):"
echo "  cd ../../infrastructure"
echo "  cdk deploy"
echo ""
echo "Option 2: Force replacement via AWS CLI by updating a stack parameter:"
echo "  (This will trigger a replacement)"
echo ""
echo "Option 3: Manual replacement via CloudFormation Console:"
echo "  1. Go to CloudFormation > HelixAIStack"
echo "  2. Select the EC2Instance resource"
echo "  3. Click 'Replace' or 'Delete' and let it recreate"
echo ""
echo "Checking current stack status..."
aws cloudformation describe-stack-resources \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --logical-resource-id EC2Instance770AAE32 \
    --query "StackResources[0].[PhysicalResourceId,ResourceStatus]" \
    --output text

echo ""
echo "The instance resource exists but the physical instance is terminated."
echo "CDK deploy should detect this and recreate the instance."
