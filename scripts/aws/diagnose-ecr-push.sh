#!/bin/bash
# Diagnostic script for ECR push issues

echo "🔍 Diagnosing ECR Push Issues"
echo "=============================="
echo ""

# Configuration
AWS_ACCOUNT_ID="${AWS_ACCOUNT_ID:-794270057041}"
AWS_REGION="${AWS_REGION:-us-west-1}"
ECR_REPO_NAME="helix-ai-backend"

ECR_URI="${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com"

# Check 1: AWS CLI
echo "1️⃣ Checking AWS CLI..."
if command -v aws &> /dev/null; then
    AWS_VERSION=$(aws --version)
    echo "   ✅ AWS CLI installed: $AWS_VERSION"
else
    echo "   ❌ AWS CLI not found"
    exit 1
fi

# Check 2: AWS Credentials
echo ""
echo "2️⃣ Checking AWS credentials..."
if aws sts get-caller-identity &> /dev/null; then
    ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
    USER=$(aws sts get-caller-identity --query Arn --output text)
    echo "   ✅ AWS credentials valid"
    echo "   Account: $ACCOUNT"
    echo "   Identity: $USER"
else
    echo "   ❌ AWS credentials not configured or invalid"
    echo "   Run: aws configure"
    exit 1
fi

# Check 3: ECR Repository
echo ""
echo "3️⃣ Checking ECR repository..."
if aws ecr describe-repositories --repository-names ${ECR_REPO_NAME} --region ${AWS_REGION} &> /dev/null; then
    echo "   ✅ ECR repository exists: ${ECR_REPO_NAME}"
else
    echo "   ❌ ECR repository not found: ${ECR_REPO_NAME}"
    echo "   Creating repository..."
    aws ecr create-repository --repository-name ${ECR_REPO_NAME} --region ${AWS_REGION} || {
        echo "   ❌ Failed to create repository"
        exit 1
    }
    echo "   ✅ Repository created"
fi

# Check 4: Docker
echo ""
echo "4️⃣ Checking Docker..."
if command -v docker &> /dev/null; then
    DOCKER_VERSION=$(docker --version)
    echo "   ✅ Docker installed: $DOCKER_VERSION"
    
    if docker info &> /dev/null; then
        echo "   ✅ Docker daemon is running"
    else
        echo "   ❌ Docker daemon is not running"
        echo "   Start Docker Desktop or Docker daemon"
        exit 1
    fi
else
    echo "   ❌ Docker not found"
    exit 1
fi

# Check 5: Docker Platform
echo ""
echo "5️⃣ Checking Docker platform support..."
if docker buildx ls &> /dev/null; then
    echo "   ✅ Docker buildx available"
    docker buildx ls
else
    echo "   ⚠️  Docker buildx not available (may need to enable experimental features)"
fi

# Check 6: Network/Proxy
echo ""
echo "6️⃣ Checking network configuration..."
if [ -n "$HTTP_PROXY" ] || [ -n "$HTTPS_PROXY" ]; then
    echo "   ⚠️  Proxy detected:"
    [ -n "$HTTP_PROXY" ] && echo "   HTTP_PROXY=$HTTP_PROXY"
    [ -n "$HTTPS_PROXY" ] && echo "   HTTPS_PROXY=$HTTPS_PROXY"
    echo "   Note: Proxy may cause connection issues with ECR"
else
    echo "   ✅ No proxy configured"
fi

# Check 7: ECR Login
echo ""
echo "7️⃣ Testing ECR login..."
if aws ecr get-login-password --region ${AWS_REGION} | docker login --username AWS --password-stdin ${ECR_URI} &> /dev/null; then
    echo "   ✅ ECR login successful"
else
    echo "   ❌ ECR login failed"
    exit 1
fi

# Check 8: Image size
echo ""
echo "8️⃣ Checking for existing local images..."
if docker images ${ECR_REPO_NAME}:latest --format "{{.Size}}" | head -1 | grep -q .; then
    IMAGE_SIZE=$(docker images ${ECR_REPO_NAME}:latest --format "{{.Size}}" | head -1)
    echo "   ✅ Local image found: ${ECR_REPO_NAME}:latest ($IMAGE_SIZE)"
    echo "   Large images (>1GB) may timeout during push"
else
    echo "   ⚠️  No local image found - need to build first"
fi

# Check 9: Docker Desktop resources (if on Mac)
echo ""
if [[ "$OSTYPE" == "darwin"* ]]; then
    echo "9️⃣ Docker Desktop (macOS) recommendations:"
    echo "   - Ensure Docker Desktop has at least 4GB RAM allocated"
    echo "   - Check Settings > Resources > Advanced"
    echo "   - Disable VPN/proxy if possible during push"
    echo "   - Try: Docker Desktop > Settings > Resources > Network > Reset to factory defaults"
fi

# Summary
echo ""
echo "=============================="
echo "✅ Diagnostic checks complete"
echo ""
echo "If all checks passed but push still fails:"
echo "1. Try the improved build script with retries: ./build-and-push-backend.sh"
echo "2. Build from EC2 to avoid local network issues: ./scripts/aws/build-and-push-from-ec2.sh"
echo "3. Check Docker Desktop logs for network errors"
echo "4. Try pushing a smaller test image first to verify connectivity"





