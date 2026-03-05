#!/bin/bash
# Alternative deployment: Build Docker image directly on EC2 instance
# This avoids Docker Desktop proxy issues when pushing large images

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

print_info "Building Docker image on EC2 instance..."
print_info "This avoids Docker Desktop proxy issues"
print_info ""

# Get EC2 instance ID
EC2_INSTANCE_ID=$(aws cloudformation describe-stacks \
    --stack-name "${STACK_NAME}" \
    --region "${AWS_REGION}" \
    --query "Stacks[0].Outputs[?OutputKey=='EC2InstanceId'].OutputValue" \
    --output text 2>/dev/null || echo "")

if [ -z "${EC2_INSTANCE_ID}" ]; then
    print_error "EC2 instance ID not found. Deploy infrastructure first."
    exit 1
fi

print_info "EC2 Instance ID: ${EC2_INSTANCE_ID}"
print_info ""

# Create a tarball of the project (excluding unnecessary files)
print_info "Creating project tarball..."
cd "${PROJECT_ROOT}"

# Create a temporary directory with only what we need
TEMP_DIR=$(mktemp -d)
trap "rm -rf ${TEMP_DIR}" EXIT

# Copy necessary directories
cp -r backend "${TEMP_DIR}/"
cp -r agents "${TEMP_DIR}/"
cp -r tools "${TEMP_DIR}/"
cp backend/Dockerfile "${TEMP_DIR}/"

# Create tarball
TARBALL="${TEMP_DIR}/helix-backend.tar.gz"
tar czf "${TARBALL}" -C "${TEMP_DIR}" backend agents tools Dockerfile

print_info "Uploading project files to EC2 instance..."
# Use SSM to copy files
aws ssm send-command \
    --instance-ids "${EC2_INSTANCE_ID}" \
    --region "${AWS_REGION}" \
    --document-name "AWS-RunShellScript" \
    --parameters "commands=[
        'mkdir -p /tmp/helix-build',
        'echo \"Tarball upload starting...\"'
    ]" \
    --output text > /dev/null

# Wait a moment
sleep 2

# Copy tarball using SCP or SSM (simpler: upload to S3 then download on EC2)
S3_BUCKET="helix-build-${AWS_ACCOUNT_ID}"
aws s3 mb "s3://${S3_BUCKET}" --region "${AWS_REGION}" 2>/dev/null || true
aws s3 cp "${TARBALL}" "s3://${S3_BUCKET}/helix-backend.tar.gz" --region "${AWS_REGION}"

print_info "Files uploaded to S3. Building on EC2..."
print_info ""

# Create build script on EC2
BUILD_SCRIPT=$(cat <<'EOF'
#!/bin/bash
set -e

cd /opt/helix-ai
rm -rf build || true
mkdir -p build
cd build

# Download from S3
aws s3 cp s3://helix-build-794270057041/helix-backend.tar.gz .
tar xzf helix-backend.tar.gz

# Build Docker image
docker build -f Dockerfile -t helix-ai-backend:latest .

# Stop old container if running
docker stop helix-ai-backend || true
docker rm helix-ai-backend || true

# Start new container
docker run -d \
    --name helix-ai-backend \
    --restart unless-stopped \
    -p 8001:8001 \
    --env-file /opt/helix-ai/.env \
    helix-ai-backend:latest

echo "Build and deployment complete!"
EOF
)

# Execute build on EC2 via SSM
print_info "Executing build on EC2 instance..."
aws ssm send-command \
    --instance-ids "${EC2_INSTANCE_ID}" \
    --region "${AWS_REGION}" \
    --document-name "AWS-RunShellScript" \
    --parameters "commands=[${BUILD_SCRIPT//$'\n'/ }]" \
    --output text

print_info ""
print_success "Build command sent to EC2 instance!"
print_info "The build is running in the background."
print_info "Check logs with:"
print_info "  aws ssm get-command-invocation --command-id <command-id> --instance-id ${EC2_INSTANCE_ID}"
