#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   ./scripts/aws/s3_sync_frontend.sh <AWS_REGION> <S3_BUCKET_NAME> [<CLOUDFRONT_DIST_ID>]
#
# Example:
#   ./scripts/aws/s3_sync_frontend.sh us-east-1 helix-frontend-bucket E123ABC456DEF
#
# Requirements:
# - AWS CLI v2 configured (aws configure)
# - Node.js to build the frontend

REGION="${1:-}"
BUCKET="${2:-}"
CLOUDFRONT_ID="${3:-}"

if [[ -z "$REGION" || -z "$BUCKET" ]]; then
  echo "Usage: $0 <AWS_REGION> <S3_BUCKET_NAME> [<CLOUDFRONT_DIST_ID>]"
  exit 1
fi

echo "Building frontend..."
pushd frontend >/dev/null
npm ci
npm run build
popd >/dev/null

echo "Syncing to s3://${BUCKET} ..."
aws s3 sync frontend/dist "s3://${BUCKET}" --region "$REGION" --delete --acl public-read --cache-control "public,max-age=31536000,immutable"

# Ensure index.html is not cached aggressively
aws s3 cp "frontend/dist/index.html" "s3://${BUCKET}/index.html" --region "$REGION" --acl public-read --cache-control "no-cache, no-store, must-revalidate" --content-type "text/html"

if [[ -n "$CLOUDFRONT_ID" ]]; then
  echo "Creating CloudFront invalidation..."
  aws cloudfront create-invalidation --distribution-id "$CLOUDFRONT_ID" --paths "/*" >/dev/null
  echo "Invalidation requested."
fi

echo "Frontend deployed to S3."


