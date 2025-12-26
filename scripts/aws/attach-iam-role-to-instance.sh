#!/bin/bash
# Attach an IAM role to an existing EC2 instance for S3 access

set -e

INSTANCE_ID="${1:-${HELIX_EC2_INSTANCE_ID}}"
IAM_ROLE_NAME="${2:-${HELIX_EC2_IAM_ROLE_NAME:-HelixAI-EC2-S3Role}}"
REGION="${AWS_REGION:-us-west-1}"

if [ -z "$INSTANCE_ID" ]; then
    echo "Usage: $0 <instance-id> [iam-role-name]"
    echo "   or set HELIX_EC2_INSTANCE_ID and HELIX_EC2_IAM_ROLE_NAME"
    exit 1
fi

echo "Attaching IAM role to EC2 instance..."
echo "Instance ID: $INSTANCE_ID"
echo "IAM Role: $IAM_ROLE_NAME"
echo "Region: $REGION"
echo ""

# Check if instance profile exists
if aws iam get-instance-profile --instance-profile-name "$IAM_ROLE_NAME" --region "$REGION" &>/dev/null; then
    echo "✅ Instance profile exists: $IAM_ROLE_NAME"
    INSTANCE_PROFILE_ARN=$(aws iam get-instance-profile --instance-profile-name "$IAM_ROLE_NAME" --query 'InstanceProfile.Arn' --output text)
    echo "   ARN: $INSTANCE_PROFILE_ARN"
else
    echo "❌ Instance profile not found: $IAM_ROLE_NAME"
    echo ""
    echo "To create the IAM role and instance profile, run:"
    echo ""
    echo "1. Create IAM role:"
    echo "   aws iam create-role --role-name $IAM_ROLE_NAME \\"
    echo "     --assume-role-policy-document '{\"Version\":\"2012-10-17\",\"Statement\":[{\"Effect\":\"Allow\",\"Principal\":{\"Service\":\"ec2.amazonaws.com\"},\"Action\":\"sts:AssumeRole\"}]}'"
    echo ""
    echo "2. Attach S3 access policy:"
    echo "   aws iam attach-role-policy --role-name $IAM_ROLE_NAME \\"
    echo "     --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess"
    echo ""
    echo "3. Create instance profile:"
    echo "   aws iam create-instance-profile --instance-profile-name $IAM_ROLE_NAME"
    echo "   aws iam add-role-to-instance-profile --instance-profile-name $IAM_ROLE_NAME --role-name $IAM_ROLE_NAME"
    echo ""
    exit 1
fi

# Attach instance profile to instance
echo ""
echo "Attaching instance profile to EC2 instance..."
aws ec2 associate-iam-instance-profile \
    --instance-id "$INSTANCE_ID" \
    --iam-instance-profile "Name=$IAM_ROLE_NAME" \
    --region "$REGION"

echo ""
echo "✅ IAM role attached successfully!"
echo ""
echo "The instance will now have S3 access via IAM role credentials."
echo "Note: It may take a few seconds for the credentials to be available."




