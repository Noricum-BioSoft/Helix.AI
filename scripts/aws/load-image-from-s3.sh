#!/bin/bash
# Manually load Docker image from S3 on EC2 instance
# Use this if SSM isn't working yet

set -euo pipefail

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m'

print_info() { echo -e "${BLUE}[INFO]${NC} $1"; }
print_success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }

# Load config
CONFIG_FILE="${1:-deploy.config}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/${CONFIG_FILE}"

STACK_NAME="${STACK_NAME:-HelixAIStack}"
S3_BUCKET="helix-images-${AWS_ACCOUNT_ID}-${AWS_REGION}"

# Get EC2 instance info
EC2_INSTANCE_ID=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2InstanceId'].OutputValue" \
    --output text 2>/dev/null || echo "")

if [ -z "${EC2_INSTANCE_ID}" ]; then
    echo "❌ EC2 instance ID not found"
    exit 1
fi

print_info "EC2 Instance ID: ${EC2_INSTANCE_ID}"
print_info "S3 Image: s3://${S3_BUCKET}/helix-backend-latest.tar.gz"
print_info ""

# Check if image exists in S3
if ! aws s3 ls "s3://${S3_BUCKET}/helix-backend-latest.tar.gz" --region "${AWS_REGION}" >/dev/null 2>&1; then
    echo "❌ Image not found in S3. Run push-via-s3.sh first."
    exit 1
fi

print_info "Image found in S3. Creating load script..."
print_info ""

# Create a simple script that can be run on EC2
cat > /tmp/load-helix-image.sh <<EOF
#!/bin/bash
set -e
cd /opt/helix-ai

echo "Downloading image from S3..."
aws s3 cp s3://${S3_BUCKET}/helix-backend-latest.tar.gz /tmp/helix-backend.tar.gz --region ${AWS_REGION}

echo "Decompressing..."
gunzip /tmp/helix-backend.tar.gz

echo "Loading Docker image..."
docker load -i /tmp/helix-backend.tar

echo "Tagging image..."
docker tag helix-backend:latest helix-ai-backend:latest || true

echo "Stopping old container..."
docker stop helix-ai-backend || true
docker rm helix-ai-backend || true

echo "Starting new container..."
docker-compose up -d

echo "Cleaning up..."
rm -f /tmp/helix-backend.tar*

echo "✅ Done! Container should be running."
docker ps | grep helix-ai-backend || echo "⚠️  Container not running, check logs: docker logs helix-ai-backend"
EOF

print_success "Load script created at /tmp/load-helix-image.sh"
print_info ""
print_info "To load the image on EC2, you can:"
print_info ""
print_info "Option 1: Use SSM Session Manager (if available):"
print_info "  aws ssm start-session --target ${EC2_INSTANCE_ID} --region ${AWS_REGION}"
print_info "  Then run: bash /tmp/load-helix-image.sh"
print_info ""
print_info "Option 2: SSH into the instance (if you have SSH key):"
print_info "  ssh ec2-user@<EC2_PUBLIC_IP>"
print_info "  Then copy and run the script above"
print_info ""
print_info "Option 3: Wait for instance to be ready, then try SSM command again:"
print_info "  aws ssm send-command \\"
print_info "    --instance-ids ${EC2_INSTANCE_ID} \\"
print_info "    --region ${AWS_REGION} \\"
print_info "    --document-name AWS-RunShellScript \\"
print_info "    --parameters 'commands=[\"aws s3 cp s3://${S3_BUCKET}/helix-backend-latest.tar.gz /tmp/\",\"gunzip /tmp/helix-backend.tar.gz\",\"docker load -i /tmp/helix-backend.tar\",\"cd /opt/helix-ai && docker-compose up -d\"]'"
print_info ""

# Try SSM one more time with simpler commands
print_info "Attempting SSM command with simpler approach..."
aws ssm send-command \
    --instance-ids "${EC2_INSTANCE_ID}" \
    --region "${AWS_REGION}" \
    --document-name "AWS-RunShellScript" \
    --parameters "commands=[
        'cd /opt/helix-ai',
        'aws s3 cp s3://${S3_BUCKET}/helix-backend-latest.tar.gz /tmp/helix-backend.tar.gz --region ${AWS_REGION}',
        'gunzip -f /tmp/helix-backend.tar.gz',
        'docker load -i /tmp/helix-backend.tar',
        'docker tag helix-backend:latest helix-ai-backend:latest || true',
        'docker stop helix-ai-backend || true',
        'docker rm helix-ai-backend || true',
        'docker-compose up -d',
        'rm -f /tmp/helix-backend.tar*',
        'echo Image loaded and container started!'
    ]" \
    --query "Command.CommandId" \
    --output text 2>&1 | tee /tmp/ssm-command.log

if [ $? -eq 0 ]; then
    COMMAND_ID=$(cat /tmp/ssm-command.log | tail -1)
    print_success "Command sent! Command ID: ${COMMAND_ID}"
    print_info "Check status with:"
    print_info "  aws ssm get-command-invocation --command-id ${COMMAND_ID} --instance-id ${EC2_INSTANCE_ID} --region ${AWS_REGION}"
else
    print_info "SSM command failed. Use one of the manual options above."
fi
