#!/bin/bash
# Deploy backend service using AWS Copilot
# This script ensures the correct environment variables are set

set -e

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

cd "$PROJECT_ROOT"

# Set required environment variables
export AWS_PROFILE=copilot-admin
export AWS_REGION=us-west-1

# Verify we're in the right directory
if [ ! -f "copilot/.workspace" ]; then
    echo "❌ Error: Not in Copilot workspace. Expected copilot/.workspace file."
    exit 1
fi

# Verify service manifest exists
if [ ! -f "copilot/backend/manifest.yml" ]; then
    echo "❌ Error: Service manifest not found at copilot/backend/manifest.yml"
    exit 1
fi

echo "🚀 Deploying backend service to production environment..."
echo "   AWS Profile: $AWS_PROFILE"
echo "   AWS Region: $AWS_REGION"
echo "   Working Directory: $(pwd)"
echo ""

# Deploy (don't use --app flag, let Copilot detect from workspace)
copilot svc deploy --name backend --env production

echo ""
echo "✅ Deployment initiated. Monitor with:"
echo "   copilot svc status --name backend --env production"
echo "   copilot svc logs --name backend --env production --follow"


