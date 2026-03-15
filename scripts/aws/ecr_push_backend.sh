#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./scripts/aws/ecr_push_backend.sh <AWS_ACCOUNT_ID> <AWS_REGION> <ECR_REPO_NAME> [<IMAGE_TAG>]
#
# Example:
#   ./scripts/aws/ecr_push_backend.sh 123456789012 us-east-1 helix-backend v1
#
# Requirements:
# - AWS CLI v2 configured (aws configure)
# - Docker installed and logged in

ACCOUNT_ID="${1:-}"
REGION="${2:-}"
REPO_NAME="${3:-helix-backend}"
IMAGE_TAG="${4:-latest}"

if [[ -z "$ACCOUNT_ID" || -z "$REGION" ]]; then
  echo "Usage: $0 <AWS_ACCOUNT_ID> <AWS_REGION> <ECR_REPO_NAME> [<IMAGE_TAG>]"
  exit 1
fi

ECR_URI="${ACCOUNT_ID}.dkr.ecr.${REGION}.amazonaws.com"

echo "Logging in to ECR..."
aws ecr get-login-password --region "$REGION" | docker login --username AWS --password-stdin "${ECR_URI}"

echo "Ensuring ECR repository ${REPO_NAME} exists..."
aws ecr describe-repositories --repository-names "${REPO_NAME}" --region "$REGION" >/dev/null 2>&1 || \
  aws ecr create-repository --repository-name "${REPO_NAME}" --region "$REGION" >/dev/null

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

echo "Building backend image for linux/amd64..."
# Change to project root to ensure Docker build context is correct
cd "${PROJECT_ROOT}"

# ECS Fargate in this stack runs linux/amd64 tasks. Build and load the image
# explicitly for that platform to avoid manifest mismatch pull errors.
docker buildx build \
  --platform linux/amd64 \
  -f backend/Dockerfile \
  -t "${REPO_NAME}:${IMAGE_TAG}" \
  --load \
  .

echo "Tagging image..."
docker tag "${REPO_NAME}:${IMAGE_TAG}" "${ECR_URI}/${REPO_NAME}:${IMAGE_TAG}"

echo "Pushing to ECR..."
docker push "${ECR_URI}/${REPO_NAME}:${IMAGE_TAG}"

echo "Done. Image: ${ECR_URI}/${REPO_NAME}:${IMAGE_TAG}"


