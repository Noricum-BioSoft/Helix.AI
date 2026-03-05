#!/bin/bash
# Run the EC2 instance test suite on a remote instance

set -e

INSTANCE_ID="${1:-${HELIX_EC2_INSTANCE_ID}}"
KEY_FILE="${2:-${HELIX_EC2_KEY_FILE}}"
REGION="${AWS_REGION:-us-west-1}"

if [ -z "$INSTANCE_ID" ] || [ -z "$KEY_FILE" ]; then
    echo "Usage: $0 <instance-id> <key-file>"
    echo "   or set HELIX_EC2_INSTANCE_ID and HELIX_EC2_KEY_FILE"
    exit 1
fi

# Get instance IP
PUBLIC_IP=$(aws ec2 describe-instances \
    --instance-ids "$INSTANCE_ID" \
    --region "$REGION" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' \
    --output text)

if [ "$PUBLIC_IP" == "None" ] || [ -z "$PUBLIC_IP" ]; then
    echo "❌ Could not get public IP for instance $INSTANCE_ID"
    exit 1
fi

echo "Running EC2 instance test suite..."
echo "Instance ID: $INSTANCE_ID"
echo "Public IP: $PUBLIC_IP"
echo ""

# Find Python on instance
echo "Finding Python on instance..."
PYTHON_CMD=$(ssh -i "$KEY_FILE" -o StrictHostKeyChecking=no ec2-user@$PUBLIC_IP \
    "which /opt/conda/bin/python3 || which python3 || which python" 2>/dev/null)

if [ -z "$PYTHON_CMD" ]; then
    echo "❌ Python not found on instance"
    exit 1
fi

echo "Using Python: $PYTHON_CMD"
echo ""

# Copy test script to instance
TEST_SCRIPT="/tmp/test_ec2_instance.py"
echo "Copying test script to instance..."
scp -i "$KEY_FILE" -o StrictHostKeyChecking=no \
    scripts/aws/test-ec2-instance.py \
    ec2-user@$PUBLIC_IP:$TEST_SCRIPT

# Run test script
echo "Running test suite..."
echo ""
ssh -i "$KEY_FILE" -o StrictHostKeyChecking=no ec2-user@$PUBLIC_IP \
    "cd /opt/helix-tools && AWS_REGION=$REGION AWS_DEFAULT_REGION=$REGION $PYTHON_CMD $TEST_SCRIPT"

# Cleanup
ssh -i "$KEY_FILE" -o StrictHostKeyChecking=no ec2-user@$PUBLIC_IP \
    "rm -f $TEST_SCRIPT" 2>/dev/null || true

