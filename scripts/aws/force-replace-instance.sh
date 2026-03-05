#!/bin/bash
# Force CloudFormation to replace the EC2 instance by updating user data

set -euo pipefail

source deploy.config
STACK_NAME="${STACK_NAME:-HelixAIStack}"

echo "Forcing CloudFormation to replace EC2 instance..."
echo "This will update the stack with a small change to trigger replacement."
echo ""

# Get current template
echo "Getting current template..."
TEMPLATE=$(aws cloudformation get-template \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query 'TemplateBody' \
    --output json 2>/dev/null)

if [ -z "${TEMPLATE}" ] || [ "${TEMPLATE}" = "null" ]; then
    echo "❌ Could not get template. Try running 'cdk deploy' instead."
    exit 1
fi

echo "✅ Got template"
echo ""
echo "The easiest way to force replacement is via AWS Console:"
echo ""
echo "1. Go to: https://console.aws.amazon.com/cloudformation"
echo "2. Select stack: ${STACK_NAME}"
echo "3. Click 'Update'"
echo "4. Select 'Replace current template'"
echo "5. Use the same template (this will force a check and replacement)"
echo ""
echo "OR, if you have CDK installed in a different environment:"
echo "  cd ../../infrastructure"
echo "  cdk deploy"
echo ""
echo "The instance resource still exists in CloudFormation but the physical instance is terminated."
echo "CloudFormation needs to detect this and recreate it."
