#!/bin/bash
# Workaround for Docker Desktop proxy issues: Save image, upload to S3, load on EC2
# This bypasses the problematic ECR push step

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Load configuration
CONFIG_FILE="${1:-deploy.config}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [ ! -f "${SCRIPT_DIR}/${CONFIG_FILE}" ]; then
    print_error "Configuration file not found: ${SCRIPT_DIR}/${CONFIG_FILE}"
    exit 1
fi

source "${SCRIPT_DIR}/${CONFIG_FILE}"

STACK_NAME="${STACK_NAME:-HelixAIStack}"
IMAGE_NAME="${ECR_REPOSITORY}:${IMAGE_TAG:-latest}"
S3_BUCKET="helix-images-${AWS_ACCOUNT_ID}-${AWS_REGION}"

print_info "Workaround: Uploading Docker image via S3 to bypass ECR push issues"
print_info "Image: ${IMAGE_NAME}"
print_info ""

# Build image if it doesn't exist
if ! docker images "${IMAGE_NAME}" | grep -q "${ECR_REPOSITORY}"; then
    print_info "Image not found locally. Building..."
    cd "${PROJECT_ROOT}"
    "${SCRIPT_DIR}/ecr_push_backend.sh" "${AWS_ACCOUNT_ID}" "${AWS_REGION}" "${ECR_REPOSITORY}" "${IMAGE_TAG:-latest}" 2>&1 | grep -v "Pushing\|Pushed" || true
    # Don't fail if push fails - we'll use S3 instead
fi

print_info "Saving Docker image to tar file..."
IMAGE_TAR="/tmp/helix-backend-${IMAGE_TAG:-latest}.tar"
docker save "${IMAGE_NAME}" -o "${IMAGE_TAR}"

# Compress to save upload time
print_info "Compressing image..."
gzip -f "${IMAGE_TAR}"
IMAGE_TAR_GZ="${IMAGE_TAR}.gz"

print_info "Creating S3 bucket for images..."
aws s3 mb "s3://${S3_BUCKET}" --region "${AWS_REGION}" 2>/dev/null || print_info "Bucket already exists"

print_info "Uploading image to S3 (this may take a few minutes)..."
aws s3 cp "${IMAGE_TAR_GZ}" "s3://${S3_BUCKET}/helix-backend-${IMAGE_TAG:-latest}.tar.gz" --region "${AWS_REGION}"

print_success "Image uploaded to S3: s3://${S3_BUCKET}/helix-backend-${IMAGE_TAG:-latest}.tar.gz"

# Get EC2 instance ID
EC2_INSTANCE_ID=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2InstanceId'].OutputValue" \
    --output text 2>/dev/null || echo "")

if [ -z "${EC2_INSTANCE_ID}" ]; then
    print_warning "EC2 instance ID not found. Image is in S3, you can manually load it later."
    exit 0
fi

print_info ""
print_info "Loading image on EC2 instance..."
print_info "This will download, decompress, and load the image"

# Create script to load image on EC2
# Use a temporary file for the script to avoid parameter parsing issues
TEMP_SCRIPT=$(mktemp)
cat > "${TEMP_SCRIPT}" <<EOF
#!/bin/bash
set -e
cd /opt/helix-ai

# Download image from S3
aws s3 cp s3://${S3_BUCKET}/helix-backend-${IMAGE_TAG:-latest}.tar.gz /tmp/helix-backend.tar.gz --region ${AWS_REGION}

# Decompress
gunzip /tmp/helix-backend.tar.gz

# Load image
docker load -i /tmp/helix-backend.tar

# Tag it for docker-compose
docker tag ${IMAGE_NAME} helix-ai-backend:latest || true

# Stop and remove old container
docker stop helix-ai-backend || true
docker rm helix-ai-backend || true

# Start new container
docker-compose up -d

# Cleanup
rm -f /tmp/helix-backend.tar*

echo "Image loaded and container restarted!"
EOF

# Upload script to S3 first, then execute it
print_info "Uploading load script to S3..."
aws s3 cp "${TEMP_SCRIPT}" "s3://${S3_BUCKET}/load-image.sh" --region "${AWS_REGION}"

# Execute via SSM - download script from S3 and run it
print_info "Executing load script on EC2..."
COMMAND_ID=$(aws ssm send-command \
    --instance-ids "${EC2_INSTANCE_ID}" \
    --region "${AWS_REGION}" \
    --document-name "AWS-RunShellScript" \
    --parameters "commands=[
        'aws s3 cp s3://${S3_BUCKET}/load-image.sh /tmp/load-image.sh --region ${AWS_REGION}',
        'chmod +x /tmp/load-image.sh',
        '/tmp/load-image.sh'
    ]" \
    --query "Command.CommandId" \
    --output text)

# Cleanup temp script
rm -f "${TEMP_SCRIPT}"

print_info ""
print_info "Command sent. Command ID: ${COMMAND_ID}"
print_info "Monitor progress with:"
print_info "  aws ssm get-command-invocation --command-id ${COMMAND_ID} --instance-id ${EC2_INSTANCE_ID} --region ${AWS_REGION}"

# Cleanup local files
rm -f "${IMAGE_TAR_GZ}"

print_info ""
print_success "✅ Image uploaded and loading script sent to EC2!"
