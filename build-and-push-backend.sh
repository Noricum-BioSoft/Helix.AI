#!/bin/bash
set -e

echo "🐳 Building and Pushing Backend Docker Image to ECR"
echo "===================================================="

# Configuration
AWS_ACCOUNT_ID="794270057041"
AWS_REGION="us-west-1"
ECR_REPO_NAME="helix-ai-backend"
IMAGE_TAG="latest"
MAX_RETRIES=3
RETRY_DELAY=5

ECR_URI="${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com"
IMAGE_NAME="${ECR_URI}/${ECR_REPO_NAME}:${IMAGE_TAG}"

# Check for Docker proxy configuration (common cause of ECR push failures)
echo "🔍 Checking Docker configuration..."
if docker info 2>/dev/null | grep -qi "proxy"; then
    echo "⚠️  WARNING: Docker proxy detected - this may cause ECR push failures"
    echo "   If push fails with 'use of closed network connection', try:"
    echo "   1. Docker Desktop > Settings > Resources > Network"
    echo "   2. Add '*.dkr.ecr.*.amazonaws.com' to 'No proxy for these hosts'"
    echo "   3. Or use: ./scripts/aws/build-and-push-from-ec2.sh (builds on EC2, avoids proxy)"
    echo ""
fi

# Function to login to ECR
ecr_login() {
    echo "🔐 Logging in to ECR..."
    aws ecr get-login-password --region ${AWS_REGION} | docker login --username AWS --password-stdin ${ECR_URI}
    if [ $? -eq 0 ]; then
        echo "✅ ECR login successful"
        return 0
    else
        echo "❌ ECR login failed"
        return 1
    fi
}

# Function to push with retries and better error handling
push_with_retry() {
    local attempt=1
    while [ $attempt -le $MAX_RETRIES ]; do
        echo ""
        echo "📤 Attempt $attempt/$MAX_RETRIES: Pushing to ECR..."
        
        # Re-authenticate before each retry (ECR tokens expire)
        if [ $attempt -gt 1 ]; then
            echo "🔄 Re-authenticating to ECR..."
            ecr_login || return 1
            # Increase delay for subsequent retries (network issues may need more time)
            RETRY_DELAY=$((RETRY_DELAY * 2))
        fi
        
        # Attempt push with better error handling
        # Use DOCKER_BUILDKIT=0 to avoid potential issues with buildkit and proxy
        local push_output
        push_output=$(DOCKER_BUILDKIT=0 docker push ${IMAGE_NAME} 2>&1) && {
            echo "✅ Successfully pushed ${IMAGE_NAME}"
            return 0
        } || {
            local exit_code=$?
            echo "❌ Push attempt $attempt failed (exit code: $exit_code)"
            echo "$push_output" | tail -5  # Show last few lines of error
            
            # Check for specific error types in output
            if echo "$push_output" | grep -qiE "timeout|connection.*closed|proxy|network"; then
                echo "⚠️  Network/proxy issue detected"
                echo "   This often happens with Docker Desktop proxy settings"
                echo "   Try disabling proxy for ECR or use EC2 build script"
            fi
            
            if [ $attempt -lt $MAX_RETRIES ]; then
                echo "⏳ Waiting ${RETRY_DELAY} seconds before retry..."
                sleep $RETRY_DELAY
            fi
            attempt=$((attempt + 1))
        fi
    done
    
    echo ""
    echo "❌ Failed to push after $MAX_RETRIES attempts"
    return 1
}

# Step 1: Build locally first (test it works)
# IMPORTANT: Build for linux/amd64 platform (required for ECS Fargate)
echo ""
echo "📦 Step 1: Building Docker image locally for linux/amd64 platform..."
docker build --platform linux/amd64 -f backend/Dockerfile -t ${ECR_REPO_NAME}:${IMAGE_TAG} .

# Step 2: Test the import works
echo ""
echo "🧪 Step 2: Testing import..."
if docker run --rm ${ECR_REPO_NAME}:${IMAGE_TAG} python -c "import backend.main_with_mcp; print('✅ Import successful!')" > /dev/null 2>&1; then
    echo "✅ Import test passed!"
else
    echo "❌ Import test failed! Fix the Docker image before pushing."
    exit 1
fi

# Step 3: Login to ECR
echo ""
ecr_login || exit 1

# Step 4: Tag for ECR
echo ""
echo "🏷️  Step 4: Tagging image for ECR..."
docker tag ${ECR_REPO_NAME}:${IMAGE_TAG} ${IMAGE_NAME}

# Step 5: Push to ECR with retries
echo ""
push_with_retry || {
    echo ""
    echo "❌ Push failed after multiple retries"
    echo ""
    echo "🔍 Common causes and solutions:"
    echo ""
    echo "🔧 Solution 1: Fix Docker Desktop Proxy (RECOMMENDED)"
    echo "   The error 'use of closed network connection' usually indicates proxy issues."
    echo "   Try disabling proxy for ECR endpoints:"
    echo "   1. Docker Desktop > Settings > Resources > Network"
    echo "   2. Add to 'No proxy for these hosts':"
    echo "      *.dkr.ecr.us-west-1.amazonaws.com"
    echo "      ${ECR_URI}"
    echo "   3. Apply & Restart"
    echo ""
    echo "🔧 Solution 2: Temporarily disable proxy"
    echo "   Unset proxy environment variables and retry:"
    echo "   unset HTTP_PROXY HTTPS_PROXY http_proxy https_proxy"
    echo "   docker push ${IMAGE_NAME}"
    echo ""
    echo "🔧 Solution 3: Use EC2 build (avoids local proxy issues)"
    echo "   Build and push from an EC2 instance:"
    echo "   ./scripts/aws/build-and-push-from-ec2.sh"
    echo ""
    echo "🔧 Solution 4: Check Docker Desktop resources"
    echo "   - Docker Desktop > Settings > Resources"
    echo "   - Increase memory to at least 4GB"
    echo "   - Increase disk space if needed"
    echo ""
    echo "🔧 Solution 5: Verify network connectivity"
    echo "   - Check AWS credentials: aws sts get-caller-identity"
    echo "   - Verify ECR repository: aws ecr describe-repositories --repository-names ${ECR_REPO_NAME} --region ${AWS_REGION}"
    echo "   - Test connectivity: curl -I https://${ECR_URI}"
    echo ""
    echo "🔧 Solution 6: Try pushing with increased verbosity"
    echo "   DOCKER_BUILDKIT=0 docker push ${IMAGE_NAME}"
    echo ""
    exit 1
}

echo ""
echo "✅ Successfully pushed ${IMAGE_NAME}"
echo ""
echo "🚀 You can now deploy the CDK stack - the ECS service will use this image."
