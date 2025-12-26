#!/bin/bash
# Alternative: Build and push Docker image from EC2 to avoid local network/proxy issues
# This script can be run on an EC2 instance to build and push the image

set -e

echo "🐳 Building and Pushing Backend Docker Image to ECR (from EC2)"
echo "==============================================================="

# Configuration
AWS_ACCOUNT_ID="${AWS_ACCOUNT_ID:-794270057041}"
AWS_REGION="${AWS_REGION:-us-west-1}"
ECR_REPO_NAME="helix-ai-backend"
IMAGE_TAG="latest"

ECR_URI="${AWS_ACCOUNT_ID}.dkr.ecr.${AWS_REGION}.amazonaws.com"
IMAGE_NAME="${ECR_URI}/${ECR_REPO_NAME}:${IMAGE_TAG}"

# Check if we're on EC2
if ! curl -s --max-time 1 http://169.254.169.254/latest/meta-data/instance-id > /dev/null 2>&1; then
    echo "⚠️  This script is designed to run on an EC2 instance"
    echo "On EC2, IAM roles provide credentials automatically"
    echo "Proceeding anyway (assuming AWS credentials are configured)..."
fi

# Login to ECR (uses IAM role credentials on EC2)
echo ""
echo "🔐 Logging in to ECR..."
aws ecr get-login-password --region ${AWS_REGION} | docker login --username AWS --password-stdin ${ECR_URI}

# Build for linux/amd64 (required for ECS Fargate)
echo ""
echo "📦 Building Docker image for linux/amd64 platform..."
docker build --platform linux/amd64 -f backend/Dockerfile -t ${IMAGE_NAME} .

# Push to ECR
echo ""
echo "📤 Pushing to ECR..."
docker push ${IMAGE_NAME}

echo ""
echo "✅ Successfully pushed ${IMAGE_NAME}"
echo ""
echo "🚀 You can now force ECS to redeploy:"
echo "   aws ecs update-service --cluster helix-ai-cluster --service HelixAIStack-ServiceD69D759B-G23PygXszFhF --force-new-deployment --region us-west-1"


