#!/bin/bash

# Helix.AI Unified Startup Script
# Clean, single-command startup for development and production

set -e  # Exit on any error

echo "🧬 Starting Helix.AI Unified System..."
echo "======================================"

# Configuration
FRONTEND_PORT=5173
BACKEND_PORT=8001
REDIS_PORT=6379

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Function to check if a port is available
check_port() {
    local port=$1
    if lsof -Pi :$port -sTCP:LISTEN -t >/dev/null 2>&1; then
        print_warning "Port $port is already in use"
        return 1
    fi
    return 0
}

# Function to wait for service to be ready
wait_for_service() {
    local url=$1
    local service_name=$2
    local max_attempts=30
    local attempt=1
    
    print_status "Waiting for $service_name..."
    
    while [ $attempt -le $max_attempts ]; do
        if curl -s "$url" > /dev/null 2>&1; then
            print_success "$service_name is ready!"
            return 0
        fi
        
        echo "   Attempt $attempt/$max_attempts..."
        sleep 2
        attempt=$((attempt + 1))
    done
    
    print_error "$service_name failed to start within expected time"
    return 1
}

# Check prerequisites
print_status "Checking prerequisites..."

# Check if required directories exist
if [ ! -d "backend" ] || [ ! -d "frontend" ]; then
    print_error "Backend or frontend directory not found"
    exit 1
fi

# Check Python
if ! command -v python &> /dev/null; then
    print_error "Python is not installed"
    exit 1
fi

# Check Node.js
if ! command -v node &> /dev/null; then
    print_error "Node.js is not installed"
    exit 1
fi

# Check npm
if ! command -v npm &> /dev/null; then
    print_error "npm is not installed"
    exit 1
fi

print_success "Prerequisites check passed"

# Check ports
print_status "Checking ports..."
check_port $BACKEND_PORT || exit 1
check_port $FRONTEND_PORT || exit 1

# Optional: Check Redis (for enhanced session management)
USE_REDIS=false
if command -v redis-server &> /dev/null; then
    if check_port $REDIS_PORT; then
        print_status "Starting Redis for enhanced session management..."
        redis-server --port $REDIS_PORT --daemonize yes
        USE_REDIS=true
        print_success "Redis started on port $REDIS_PORT"
    else
        print_warning "Redis port $REDIS_PORT is in use, using file-based sessions"
    fi
else
    print_warning "Redis not found, using file-based sessions"
fi

# Install backend dependencies
print_status "Installing backend dependencies..."
cd backend

# Install from requirements.txt
if [ -f "requirements.txt" ]; then
    python -m pip install -r requirements.txt
else
    print_warning "requirements.txt not found, installing essential packages..."
    python -m pip install fastapi uvicorn langchain langgraph pydantic requests pandas seaborn plotly
fi

# Test plasmid visualizer import
print_status "Testing plasmid visualizer..."
export PYTHONPATH="../tools:$PYTHONPATH"
python -c "import plasmid_visualizer; print('✅ Plasmid visualizer ready')" || {
    print_error "Failed to import plasmid visualizer"
    exit 1
}

cd ..

# Install frontend dependencies
print_status "Installing frontend dependencies..."
cd frontend
npm install
cd ..

# Start backend server
print_status "Starting backend server..."
cd backend
export PYTHONPATH="../tools:$PYTHONPATH"

# Set environment variables
if [ "$USE_REDIS" = true ]; then
    export REDIS_URL="redis://localhost:$REDIS_PORT"
    print_status "Using Redis for session management"
else
    print_status "Using file-based session management"
fi

# Start the enhanced MCP backend
python main_with_mcp.py &
BACKEND_PID=$!
cd ..

# Wait for backend
if wait_for_service "http://localhost:$BACKEND_PORT/health" "Backend server"; then
    print_success "Backend server is running on http://localhost:$BACKEND_PORT"
else
    print_error "Backend server failed to start"
    kill $BACKEND_PID 2>/dev/null
    exit 1
fi

# Start frontend server
print_status "Starting frontend server..."
cd frontend
npm run dev &
FRONTEND_PID=$!
cd ..

# Wait for frontend
if wait_for_service "http://localhost:$FRONTEND_PORT" "Frontend server"; then
    print_success "Frontend server is running on http://localhost:$FRONTEND_PORT"
else
    print_error "Frontend server failed to start"
    kill $BACKEND_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    exit 1
fi

echo ""
echo "🎉 Helix.AI Unified System is ready!"
echo "====================================="
echo "🌐 Frontend: http://localhost:$FRONTEND_PORT"
echo "🔧 Backend:  http://localhost:$BACKEND_PORT"
echo "📚 API Docs: http://localhost:$BACKEND_PORT/docs"
echo "🧬 Plasmid Visualizer: Available in the UI"
echo "📊 Session Management: $([ "$USE_REDIS" = true ] && echo "Redis" || echo "File-based")"
echo ""
echo "Press Ctrl+C to stop all servers"

# Cleanup function
cleanup() {
    echo ""
    print_status "Stopping servers..."
    kill $BACKEND_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    
    if [ "$USE_REDIS" = true ]; then
        print_status "Stopping Redis..."
        redis-cli -p $REDIS_PORT shutdown 2>/dev/null || true
    fi
    
    print_success "All servers stopped"
    exit 0
}

# Set up signal handlers
trap cleanup SIGINT SIGTERM

# Wait for user to stop
wait
