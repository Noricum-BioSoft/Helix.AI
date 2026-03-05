#!/bin/bash
# Simple script to run deployed backend tests

set -e

# Default backend URL (can be overridden via environment variable)
BACKEND_URL="${BACKEND_URL:-http://HelixA-ALBAE-D7BksiQIynZb-1051248867.us-west-1.elb.amazonaws.com}"

echo "🧪 Running Deployed Backend Tests"
echo "=================================="
echo "Backend URL: $BACKEND_URL"
echo ""

# Export for pytest
export BACKEND_URL

# Run pytest with the deployed backend tests
# Note: -m integration explicitly includes integration tests (pytest.ini excludes them by default)
python -m pytest tests/integration/test_deployed_backend.py -v --tb=short -m integration

echo ""
echo "✅ Tests complete!"

