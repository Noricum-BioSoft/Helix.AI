#!/bin/bash
# Push to ECR with workaround for Docker Desktop proxy issues
# This script tries to push with increased timeouts and without proxy if possible

set -euo pipefail

# Load configuration
CONFIG_FILE="${1:-deploy.config}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

source "${SCRIPT_DIR}/${CONFIG_FILE}"

ECR_URI="${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com"
IMAGE_NAME="${ECR_REPOSITORY}:${IMAGE_TAG:-latest}"
FULL_IMAGE_NAME="${ECR_URI}/${ECR_REPOSITORY}:${IMAGE_TAG:-latest}"

echo "Logging in to ECR..."
aws ecr get-login-password --region "${AWS_REGION}" | docker login --username AWS --password-stdin "${ECR_URI}"

echo ""
echo "Trying to push with workarounds for Docker Desktop proxy issues..."
echo "If this fails, try:"
echo "  1. Restart Docker Desktop"
echo "  2. Wait a few minutes and try again"
echo "  3. Check your network connection"
echo ""

# Try pushing with DOCKER_BUILDKIT=0 to use legacy builder
export DOCKER_BUILDKIT=0

# Push with retry logic
MAX_ATTEMPTS=5
ATTEMPT=1

while [ $ATTEMPT -le $MAX_ATTEMPTS ]; do
    echo "Push attempt $ATTEMPT of $MAX_ATTEMPTS..."
    
    # Try with timeout command if available
    if command -v gtimeout &> /dev/null || command -v timeout &> /dev/null; then
        TIMEOUT_CMD="timeout 600"
    else
        TIMEOUT_CMD=""
    fi
    
    if $TIMEOUT_CMD docker push "${FULL_IMAGE_NAME}" 2>&1 | tee /tmp/docker-push.log; then
        echo "✅ Push succeeded!"
        exit 0
    else
        if grep -q "use of closed network connection\|broken pipe\|EOF" /tmp/docker-push.log; then
            echo "⚠️  Network/proxy error detected. This is often a Docker Desktop issue."
            if [ $ATTEMPT -lt $MAX_ATTEMPTS ]; then
                WAIT_TIME=$((ATTEMPT * 10))
                echo "Waiting ${WAIT_TIME} seconds before retry..."
                sleep $WAIT_TIME
            fi
        else
            echo "❌ Different error detected. Check the log above."
            exit 1
        fi
    fi
    
    ATTEMPT=$((ATTEMPT + 1))
done

echo "❌ Push failed after $MAX_ATTEMPTS attempts"
echo ""
echo "Recommended solutions:"
echo "  1. Restart Docker Desktop: Docker Desktop menu > Restart"
echo "  2. Try again in a few minutes (may be temporary network issue)"
echo "  3. Check Docker Desktop settings > Resources > Advanced"
echo "  4. Try disabling Docker Desktop proxy in settings"
echo ""
exit 1
