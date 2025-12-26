#!/bin/bash
set -e

echo "🧪 Testing Docker Build and Run Locally"
echo "========================================"

# Build the image
echo "📦 Building Docker image..."
docker build -f backend/Dockerfile -t helix-ai-backend:test .

# Test that the container can start and imports work
echo ""
echo "🐳 Testing container startup..."
echo "Running: docker run --rm helix-ai-backend:test python -c 'import backend.main_with_mcp; print(\"✅ Import successful!\")'"

if docker run --rm helix-ai-backend:test python -c 'import backend.main_with_mcp; print("✅ Import successful!")' 2>&1; then
    echo ""
    echo "✅ Import test passed!"
else
    echo ""
    echo "❌ Import test failed!"
    exit 1
fi

# Test health endpoint (run in background, check, then stop)
echo ""
echo "🏥 Testing health endpoint..."
CONTAINER_ID=$(docker run -d -p 8002:8001 helix-ai-backend:test)

# Wait a bit for the server to start
sleep 5

if curl -f http://localhost:8002/health > /dev/null 2>&1; then
    echo "✅ Health check passed!"
else
    echo "❌ Health check failed!"
    docker logs $CONTAINER_ID
    docker stop $CONTAINER_ID > /dev/null 2>&1
    exit 1
fi

docker stop $CONTAINER_ID > /dev/null 2>&1
echo ""
echo "✅ All tests passed! Docker image is working correctly."
