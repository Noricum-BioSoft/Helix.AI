#!/bin/bash

# Helix.AI Application Startup Script
# This script starts both the MCP server and the FastAPI server, along with the frontend

echo "🚀 Starting Helix.AI Application..."
echo "========================================"

# Function to check if a port is in use
check_port() {
    if lsof -Pi :$1 -sTCP:LISTEN -t >/dev/null ; then
        echo "⚠️  Port $1 is already in use"
        return 1
    else
        return 0
    fi
}

# Function to wait for a service to be ready
wait_for_service() {
    local url=$1
    local max_attempts=30
    local attempt=1
    
    echo "⏳ Waiting for service at $url..."
    
    while [ $attempt -le $max_attempts ]; do
        if curl -s "$url" > /dev/null 2>&1; then
            echo "✅ Service is ready!"
            return 0
        fi
        
        echo "   Attempt $attempt/$max_attempts..."
        sleep 2
        attempt=$((attempt + 1))
    done
    
    echo "❌ Service failed to start within expected time"
    return 1
}

# Function to check if a process is running
check_process() {
    local process_name=$1
    if pgrep -f "$process_name" > /dev/null; then
        return 0
    else
        return 1
    fi
}

# Check if required directories exist
if [ ! -d "backend" ]; then
    echo "❌ Backend directory not found"
    exit 1
fi

if [ ! -d "frontend" ]; then
    echo "❌ Frontend directory not found"
    exit 1
fi

# Check if Python is available
if ! command -v python &> /dev/null; then
    echo "❌ Python is not installed"
    exit 1
fi

# Check if Node.js is available
if ! command -v node &> /dev/null; then
    echo "❌ Node.js is not installed"
    exit 1
fi

# Check if npm is available
if ! command -v npm &> /dev/null; then
    echo "❌ npm is not installed"
    exit 1
fi

echo "✅ Prerequisites check passed"

# Check ports
check_port 8001 || exit 1
check_port 5173 || exit 1

# Install backend dependencies
echo "📦 Installing backend dependencies..."
cd backend

# Install essential packages for the application
echo "   Installing essential Python packages..."
python -m pip install fastapi uvicorn langchain langgraph pydantic requests pandas

# Try to install requirements.txt if it exists
if [ -f "requirements.txt" ]; then
    echo "   Installing from requirements.txt..."
    python -m pip install -r requirements.txt 2>/dev/null || echo "   ⚠️  Some packages from requirements.txt failed to install, continuing with essential packages"
else
    echo "   ⚠️  requirements.txt not found, using essential packages only"
fi

cd ..

# Install frontend dependencies
echo "📦 Installing frontend dependencies..."
cd frontend
if [ -f "package.json" ]; then
    npm install
else
    echo "   ⚠️  package.json not found, skipping frontend dependency installation"
fi
cd ..

# Start FastAPI server (main server with plasmid visualization)
echo "🔧 Starting FastAPI server with plasmid visualization..."
cd backend

# Set up Python path for tools
export PYTHONPATH="../tools:$PYTHONPATH"

# Test if plasmid visualizer can be imported
echo "   Testing plasmid visualizer import..."
python -c "import plasmid_visualizer; print('✅ Plasmid visualizer imported successfully')" 2>/dev/null || {
    echo "   ❌ Failed to import plasmid visualizer"
    echo "   Current directory: $(pwd)"
    echo "   PYTHONPATH: $PYTHONPATH"
    exit 1
}

echo "   Starting FastAPI server..."
python main_with_mcp.py &
BACKEND_PID=$!
cd ..

# Wait for FastAPI server to be ready
if wait_for_service "http://localhost:8001/health"; then
    echo "✅ FastAPI server is running on http://localhost:8001"
else
    echo "❌ FastAPI server failed to start"
    kill $BACKEND_PID 2>/dev/null
    exit 1
fi

# Start frontend server
echo "🎨 Starting frontend development server..."
cd frontend
npm run dev &
FRONTEND_PID=$!
cd ..

# Wait for frontend to be ready
if wait_for_service "http://localhost:5173"; then
    echo "✅ Frontend server is running on http://localhost:5173"
else
    echo "❌ Frontend server failed to start"
    kill $BACKEND_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    exit 1
fi

echo ""
echo "🎉 Helix.AI Application is now running!"
echo "=========================================="
echo "🌐 Frontend: http://localhost:5173"
echo "🔧 FastAPI:  http://localhost:8001"
echo "📚 API Docs: http://localhost:8001/docs"
echo "🧬 Plasmid Visualizer: Available in the UI"
echo ""
echo "Press Ctrl+C to stop all servers"

# Function to cleanup on exit
cleanup() {
    echo ""
    echo "🛑 Stopping servers..."
    kill $BACKEND_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    echo "✅ Servers stopped"
    exit 0
}

# Set up signal handlers
trap cleanup SIGINT SIGTERM

# Wait for user to stop the application
wait 