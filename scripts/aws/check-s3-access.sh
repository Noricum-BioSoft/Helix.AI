#!/bin/bash
# Check S3 bucket access from ECS task role perspective
# This script helps diagnose S3 permission issues

set -e

BUCKET_NAME="${1:-noricum-ngs-data}"
TASK_ROLE_ARN="${2:-}"

if [ -z "$TASK_ROLE_ARN" ]; then
    echo "Usage: $0 [bucket-name] [task-role-arn]"
    echo "Example: $0 noricum-ngs-data arn:aws:iam::794270057041:role/HelixAIStack-TaskRole30FC0FBB-kTaxcfQAiUut"
    exit 1
fi

echo "=========================================="
echo "S3 Access Diagnostics"
echo "=========================================="
echo "Bucket: $BUCKET_NAME"
echo "Task Role ARN: $TASK_ROLE_ARN"
echo ""

# Extract account ID and role name
ACCOUNT_ID=$(echo "$TASK_ROLE_ARN" | cut -d: -f5)
ROLE_NAME=$(echo "$TASK_ROLE_ARN" | cut -d/ -f2)

echo "1. Checking if bucket exists and is accessible..."
if aws s3 ls "s3://$BUCKET_NAME/" > /dev/null 2>&1; then
    echo "   ✅ Bucket is accessible from current credentials"
else
    echo "   ❌ Cannot access bucket from current credentials"
    echo "   Error details:"
    aws s3 ls "s3://$BUCKET_NAME/" 2>&1 || true
fi
echo ""

echo "2. Getting bucket region..."
BUCKET_REGION=$(aws s3api get-bucket-location --bucket "$BUCKET_NAME" --query 'LocationConstraint' --output text 2>/dev/null || echo "us-east-1")
if [ "$BUCKET_REGION" == "None" ] || [ -z "$BUCKET_REGION" ]; then
    BUCKET_REGION="us-east-1"
fi
echo "   Bucket region: $BUCKET_REGION"
echo ""

echo "3. Checking bucket policy..."
if aws s3api get-bucket-policy --bucket "$BUCKET_NAME" --region "$BUCKET_REGION" > /dev/null 2>&1; then
    echo "   ✅ Bucket has a policy"
    aws s3api get-bucket-policy --bucket "$BUCKET_NAME" --region "$BUCKET_REGION" --query 'Policy' --output text | python3 -m json.tool 2>/dev/null || \
        aws s3api get-bucket-policy --bucket "$BUCKET_NAME" --region "$BUCKET_REGION" --query 'Policy' --output text
else
    echo "   ℹ️  No bucket policy found (or access denied to view policy)"
fi
echo ""

echo "4. Checking IAM role policies..."
echo "   Checking inline policies..."
aws iam list-role-policies --role-name "$ROLE_NAME" --query 'PolicyNames' --output table 2>/dev/null || echo "   ❌ Cannot list inline policies"
echo ""
echo "   Checking attached managed policies..."
aws iam list-attached-role-policies --role-name "$ROLE_NAME" --query 'AttachedPolicies' --output table 2>/dev/null || echo "   ❌ Cannot list attached policies"
echo ""

echo "5. Testing S3 operations with task role (if you can assume it)..."
echo "   To test with the task role, you would need to:"
echo "   1. Assume the role using: aws sts assume-role --role-arn '$TASK_ROLE_ARN' --role-session-name test-session"
echo "   2. Use the temporary credentials to test S3 access"
echo ""

echo "6. Recommendations:"
echo "   - If bucket is in a different AWS account, you need:"
echo "     1. Bucket policy allowing the task role ARN"
echo "     2. IAM role trust relationship (if cross-account)"
echo "   - If bucket is in the same account, ensure:"
echo "     1. IAM role has s3:GetObject, s3:HeadObject, s3:PutObject permissions"
echo "     2. Bucket policy doesn't deny the role"
echo "     3. Bucket ACLs allow access (if used)"
echo ""


