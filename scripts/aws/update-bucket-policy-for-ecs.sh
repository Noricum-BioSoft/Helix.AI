#!/bin/bash
# Update bucket policy to allow ECS task role access
# This script adds the ECS task role to the bucket policy

set -e

BUCKET_NAME="${1:-noricum-ngs-data}"
TASK_ROLE_ARN="${2:-arn:aws:iam::794270057041:role/HelixAIStack-TaskRole30FC0FBB-kTaxcfQAiUut}"
REGION="${3:-us-east-1}"

if [ -z "$BUCKET_NAME" ]; then
    echo "Usage: $0 [bucket-name] [task-role-arn] [region]"
    echo "Example: $0 noricum-ngs-data arn:aws:iam::794270057041:role/HelixAIStack-TaskRole30FC0FBB-kTaxcfQAiUut us-east-1"
    exit 1
fi

echo "=========================================="
echo "Updating S3 Bucket Policy"
echo "=========================================="
echo "Bucket: $BUCKET_NAME"
echo "Task Role ARN: $TASK_ROLE_ARN"
echo "Region: $REGION"
echo ""

# Get current bucket policy (if it exists)
CURRENT_POLICY=$(aws s3api get-bucket-policy --bucket "$BUCKET_NAME" --region "$REGION" 2>/dev/null | jq -r '.Policy' || echo "{}")

if [ "$CURRENT_POLICY" == "{}" ] || [ -z "$CURRENT_POLICY" ]; then
    echo "No existing bucket policy found. Creating new policy..."
    # Create new policy
    # Note: ListBucket must be on bucket ARN only, object actions on bucket/* only
    NEW_POLICY=$(cat <<EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Sid": "AllowHelixAITaskRoleListBucket",
      "Effect": "Allow",
      "Principal": {
        "AWS": "$TASK_ROLE_ARN"
      },
      "Action": [
        "s3:ListBucket"
      ],
      "Resource": "arn:aws:s3:::$BUCKET_NAME"
    },
    {
      "Sid": "AllowHelixAITaskRoleObjects",
      "Effect": "Allow",
      "Principal": {
        "AWS": "$TASK_ROLE_ARN"
      },
      "Action": [
        "s3:GetObject",
        "s3:HeadObject",
        "s3:PutObject"
      ],
      "Resource": "arn:aws:s3:::$BUCKET_NAME/*"
    }
  ]
}
EOF
)
else
    echo "Found existing bucket policy. Merging new statements..."
    # Parse current policy and add new statements (ListBucket and object actions must be separate)
    BUCKET_ARN="arn:aws:s3:::$BUCKET_NAME"
    BUCKET_ARN_OBJECTS="arn:aws:s3:::$BUCKET_NAME/*"
    NEW_POLICY=$(echo "$CURRENT_POLICY" | jq --arg ROLE_ARN "$TASK_ROLE_ARN" --arg BUCKET_ARN "$BUCKET_ARN" --arg BUCKET_ARN_OBJ "$BUCKET_ARN_OBJECTS" '
      .Statement += [
        {
          "Sid": "AllowHelixAITaskRoleListBucket",
          "Effect": "Allow",
          "Principal": {
            "AWS": $ROLE_ARN
          },
          "Action": [
            "s3:ListBucket"
          ],
          "Resource": $BUCKET_ARN
        },
        {
          "Sid": "AllowHelixAITaskRoleObjects",
          "Effect": "Allow",
          "Principal": {
            "AWS": $ROLE_ARN
          },
          "Action": [
            "s3:GetObject",
            "s3:HeadObject",
            "s3:PutObject"
          ],
          "Resource": $BUCKET_ARN_OBJ
        }
      ]
    ')
fi

# Save policy to temporary file
POLICY_FILE=$(mktemp)
echo "$NEW_POLICY" > "$POLICY_FILE"
echo "$NEW_POLICY" | jq '.' > "$POLICY_FILE"  # Pretty print

echo ""
echo "New bucket policy:"
cat "$POLICY_FILE"
echo ""
echo "Applying policy to bucket..."

# Apply the policy
if aws s3api put-bucket-policy --bucket "$BUCKET_NAME" --policy "file://$POLICY_FILE" --region "$REGION"; then
    echo "✅ Successfully updated bucket policy!"
    echo ""
    echo "Verifying policy..."
    aws s3api get-bucket-policy --bucket "$BUCKET_NAME" --region "$REGION" --query 'Policy' --output text | jq '.'
else
    echo "❌ Failed to update bucket policy"
    echo "You may need:"
    echo "  1. s3:PutBucketPolicy permission"
    echo "  2. Access to the bucket owner account"
    exit 1
fi

# Cleanup
rm -f "$POLICY_FILE"

echo ""
echo "✅ Done! The ECS task role should now be able to access the bucket."

