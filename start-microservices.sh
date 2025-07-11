#!/bin/bash

# DataBloom.AI Microservices Startup Script

echo "🚀 Starting DataBloom.AI Microservices Architecture..."

# Check if Docker and Docker Compose are installed
if ! command -v docker &> /dev/null; then
    echo "❌ Docker is not installed. Please install Docker first."
    exit 1
fi

if ! docker compose version &> /dev/null; then
    echo "❌ Docker Compose is not available. Please install Docker Compose first."
    exit 1
fi

# Create shared directory structure
echo "📁 Creating shared directory structure..."
mkdir -p shared
cp requirements.txt api-gateway/
cp requirements.txt workflow-engine/
cp requirements.txt bioinformatics-services/alignment-service/
cp requirements.txt bioinformatics-services/mutation-service/
cp requirements.txt bioinformatics-services/selection-service/
cp requirements.txt bioinformatics-services/plasmid-service/
cp requirements.txt bioinformatics-services/synthesis-service/
cp requirements.txt nlp-service/
cp requirements.txt notification-service/

# Copy shared models and utils to each service
echo "📋 Copying shared components..."
cp shared/models.py api-gateway/shared/
cp shared/utils.py api-gateway/shared/
cp shared/models.py workflow-engine/shared/
cp shared/utils.py workflow-engine/shared/
cp shared/models.py bioinformatics-services/alignment-service/shared/
cp shared/utils.py bioinformatics-services/alignment-service/shared/

# Start infrastructure services first
echo "🏗️  Starting infrastructure services..."
docker compose up -d redis postgres minio

# Wait for infrastructure to be ready
echo "⏳ Waiting for infrastructure services to be ready..."
sleep 10

# Start all microservices
echo "🔧 Starting microservices..."
docker compose up -d

# Wait for services to start
echo "⏳ Waiting for services to start..."
sleep 15

# Check service health
echo "🏥 Checking service health..."
services=(
    "http://localhost:8000/health"  # API Gateway
    "http://localhost:8001/health"  # Workflow Engine
    "http://localhost:8002/health"  # Alignment Service
)

for service in "${services[@]}"; do
    echo "Checking $service..."
    if curl -f -s "$service" > /dev/null; then
        echo "✅ $service is healthy"
    else
        echo "❌ $service is not responding"
    fi
done

# Start frontend (if it exists)
if [ -d "frontend" ]; then
    echo "🌐 Starting frontend..."
    cd frontend
    npm install
    npm run dev &
    cd ..
fi

echo ""
echo "🎉 DataBloom.AI Microservices Architecture is running!"
echo ""
echo "📊 Service URLs:"
echo "   API Gateway:     http://localhost:8000"
echo "   Workflow Engine: http://localhost:8001"
echo "   Alignment:       http://localhost:8002"
echo "   Frontend:        http://localhost:5173"
echo ""
echo "📈 Monitoring:"
echo "   MinIO Console:   http://localhost:9001"
echo "   Redis:           localhost:6379"
echo "   PostgreSQL:      localhost:5432"
echo ""
echo "🔧 Management:"
echo "   View logs:       docker compose logs -f"
echo "   Stop services:   docker compose down"
echo "   Restart:         ./start-microservices.sh"
echo "" 