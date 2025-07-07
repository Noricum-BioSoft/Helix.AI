#!/bin/bash

# DataBloom.AI Application Startup Script
# This script starts both the MCP server and the FastAPI server, along with the frontend

echo "🚀 Starting DataBloom.AI Application..."
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
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed"
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
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
else
    echo "⚠️  requirements.txt not found, skipping backend dependency installation"
fi
cd ..

# Install frontend dependencies
echo "📦 Installing frontend dependencies..."
cd frontend
if [ -f "package.json" ]; then
    npm install
else
    echo "⚠️  package.json not found, skipping frontend dependency installation"
fi
cd ..

# Start MCP server
echo "🔧 Starting MCP server..."
cd backend

if [ -f "simple_mcp_server.py" ]; then
    echo "   Starting simplified MCP server..."
    echo "   Current directory: $(pwd)"
    
    # Check if conda is available and activate databloom environment
    if command -v conda &> /dev/null; then
        echo "   Activating databloom conda environment..."
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate databloom
        echo "   Python version: $(python --version)"
        echo "   Python executable: $(which python)"
    else
        echo "   ⚠️  Conda not found, using system Python"
        echo "   Python version: $(python --version)"
    fi
    
    echo "   Checking for required modules..."
    
    # Check if required modules are available
    PYTHONPATH="$PYTHONPATH:$(pwd)/../tools" python -c "import tools.mutations" 2>/dev/null && echo "   ✅ tools.mutations module found" || echo "   ❌ tools.mutations module not found"
    PYTHONPATH="$PYTHONPATH:$(pwd)/../tools" python -c "import tools.alignment" 2>/dev/null && echo "   ✅ tools.alignment module found" || echo "   ❌ tools.alignment module not found"
    PYTHONPATH="$PYTHONPATH:$(pwd)/../tools" python -c "import tools.data_science" 2>/dev/null && echo "   ✅ tools.data_science module found" || echo "   ❌ tools.data_science module not found"
    
    echo "   Starting MCP server..."
    python simple_mcp_server.py &
    MCP_PID=$!
    sleep 3  # Give MCP server time to start

    # No longer check for process name, rely on MCP server log output
    # If you want, you can tail the log or print a message here
    # echo "   MCP server started with PID $MCP_PID (check logs for details)"

else
    echo "❌ simple_mcp_server.py not found"
    echo "   Current directory contents:"
    ls -la
    exit 1
fi
cd ..

# Start FastAPI server
echo "🔧 Starting FastAPI server..."
cd backend
python main_with_mcp.py &
BACKEND_PID=$!
cd ..

# Wait for FastAPI server to be ready
if wait_for_service "http://localhost:8001/health"; then
    echo "✅ FastAPI server is running on http://localhost:8001"
else
    echo "❌ FastAPI server failed to start"
    kill $MCP_PID 2>/dev/null
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
    kill $MCP_PID 2>/dev/null
    kill $BACKEND_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    exit 1
fi

echo ""
echo "🎉 DataBloom.AI Application is now running!"
echo "=========================================="
echo "🌐 Frontend: http://localhost:5173"
echo "🔧 FastAPI:  http://localhost:8001"
echo "📚 API Docs: http://localhost:8001/docs"
echo "🔬 MCP Server: Running (PID: $MCP_PID)"
echo ""
echo "Press Ctrl+C to stop all servers"

# Function to cleanup on exit
cleanup() {
    echo ""
    echo "🛑 Stopping servers..."
    kill $MCP_PID 2>/dev/null
    kill $BACKEND_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    echo "✅ Servers stopped"
    exit 0
}

# Set up signal handlers
trap cleanup SIGINT SIGTERM

# Wait for user to stop the application
wait 