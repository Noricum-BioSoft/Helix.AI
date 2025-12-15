#!/bin/bash
# Complete setup for an existing EC2 instance: Install conda + attach IAM role

set -e

INSTANCE_ID="${1:-${HELIX_EC2_INSTANCE_ID}}"
KEY_FILE="${2:-${HELIX_EC2_KEY_FILE}}"
IAM_ROLE_NAME="${3:-${HELIX_EC2_IAM_ROLE_NAME:-HelixAI-EC2-S3Role}}"
REGION="${AWS_REGION:-us-west-1}"

if [ -z "$INSTANCE_ID" ] || [ -z "$KEY_FILE" ]; then
    echo "Usage: $0 <instance-id> <key-file> [iam-role-name]"
    echo "   or set HELIX_EC2_INSTANCE_ID, HELIX_EC2_KEY_FILE, and optionally HELIX_EC2_IAM_ROLE_NAME"
    exit 1
fi

echo "=========================================="
echo "Helix.AI EC2 Instance Setup"
echo "=========================================="
echo "Instance ID: $INSTANCE_ID"
echo "Key File: $KEY_FILE"
echo "IAM Role: $IAM_ROLE_NAME"
echo "Region: $REGION"
echo ""

# Step 1: Create IAM role if it doesn't exist
echo "Step 1: Setting up IAM role for S3 access..."
if aws iam get-instance-profile --instance-profile-name "$IAM_ROLE_NAME" &>/dev/null; then
    echo "✅ IAM instance profile already exists: $IAM_ROLE_NAME"
else
    echo "Creating IAM role and instance profile..."
    
    # Create role
    if aws iam get-role --role-name "$IAM_ROLE_NAME" &>/dev/null; then
        echo "✅ IAM role already exists: $IAM_ROLE_NAME"
    else
        echo "Creating IAM role..."
        aws iam create-role --role-name "$IAM_ROLE_NAME" \
            --assume-role-policy-document '{
                "Version": "2012-10-17",
                "Statement": [{
                    "Effect": "Allow",
                    "Principal": {"Service": "ec2.amazonaws.com"},
                    "Action": "sts:AssumeRole"
                }]
            }'
        echo "✅ IAM role created"
    fi
    
    # Attach S3 policy
    echo "Attaching S3 access policy..."
    aws iam attach-role-policy --role-name "$IAM_ROLE_NAME" \
        --policy-arn arn:aws:iam::aws:policy/AmazonS3FullAccess
    echo "✅ S3 policy attached"
    
    # Create instance profile
    echo "Creating instance profile..."
    aws iam create-instance-profile --instance-profile-name "$IAM_ROLE_NAME" 2>/dev/null || echo "Instance profile may already exist"
    aws iam add-role-to-instance-profile --instance-profile-name "$IAM_ROLE_NAME" --role-name "$IAM_ROLE_NAME" 2>/dev/null || echo "Role may already be in instance profile"
    echo "✅ Instance profile created"
fi

# Step 2: Attach IAM role to instance
echo ""
echo "Step 2: Attaching IAM role to EC2 instance..."
# Check if already attached
CURRENT_PROFILE=$(aws ec2 describe-instances \
    --instance-ids "$INSTANCE_ID" \
    --region "$REGION" \
    --query 'Reservations[0].Instances[0].IamInstanceProfile.Arn' \
    --output text 2>/dev/null || echo "None")

if [ "$CURRENT_PROFILE" != "None" ] && [ -n "$CURRENT_PROFILE" ]; then
    echo "✅ IAM role already attached: $CURRENT_PROFILE"
else
    aws ec2 associate-iam-instance-profile \
        --instance-id "$INSTANCE_ID" \
        --iam-instance-profile "Name=$IAM_ROLE_NAME" \
        --region "$REGION"
    echo "✅ IAM role attached to instance"
fi

# Step 3: Install conda and tools
echo ""
echo "Step 3: Installing conda and bioinformatics tools..."
echo "This will take 10-15 minutes..."
echo ""

# Run the conda installation script
./scripts/aws/install-conda-on-instance.sh "$INSTANCE_ID" "$KEY_FILE"

echo ""
echo "=========================================="
echo "✅ Setup Complete!"
echo "=========================================="
echo ""
echo "Your EC2 instance is now configured with:"
echo "  ✅ Conda with Python 3.8+"
echo "  ✅ Bioinformatics tools (BBTools, samtools, etc.)"
echo "  ✅ Python packages (boto3, biopython, etc.)"
echo "  ✅ IAM role for S3 access"
echo ""
echo "You can now use the instance with:"
echo "  export HELIX_EC2_INSTANCE_ID=$INSTANCE_ID"
echo "  export HELIX_USE_EC2=true"

