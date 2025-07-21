#!/bin/bash

# Helix.AI Development Startup Script
# Simplified script for development and testing

echo "🚀 Starting Helix.AI Development Environment..."
echo "=================================================="

# Check if required directories exist
if [ ! -d "backend" ] || [ ! -d "frontend" ]; then
    echo "❌ Backend or frontend directory not found"
    exit 1
fi

# Check ports
echo "🔍 Checking ports..."
if lsof -Pi :8001 -sTCP:LISTEN -t >/dev/null ; then
    echo "⚠️  Port 8001 is already in use"
    exit 1
fi

if lsof -Pi :5173 -sTCP:LISTEN -t >/dev/null ; then
    echo "⚠️  Port 5173 is already in use"
    exit 1
fi

echo "✅ Ports available"

# Install backend dependencies
echo "📦 Installing backend dependencies..."
cd backend

# Install essential packages
echo "   Installing essential packages..."
python -m pip install fastapi uvicorn langchain langgraph pydantic requests pandas seaborn

# Test plasmid visualizer import
echo "   Testing plasmid visualizer..."
export PYTHONPATH="../tools:$PYTHONPATH"
python -c "import plasmid_visualizer; print('✅ Plasmid visualizer ready')" || {
    echo "❌ Failed to import plasmid visualizer"
    exit 1
}

cd ..

# Install frontend dependencies
echo "📦 Installing frontend dependencies..."
cd frontend
npm install
cd ..

# Start backend server
echo "🔧 Starting backend server..."
cd backend
export PYTHONPATH="../tools:$PYTHONPATH"
python main_with_mcp.py &
BACKEND_PID=$!
cd ..

# Wait for backend
echo "⏳ Waiting for backend server..."
sleep 5
if curl -s http://localhost:8001/health > /dev/null; then
    echo "✅ Backend server is running on http://localhost:8001"
else
    echo "❌ Backend server failed to start"
    kill $BACKEND_PID 2>/dev/null
    exit 1
fi

# Start frontend server
echo "🎨 Starting frontend server..."
cd frontend
npm run dev &
FRONTEND_PID=$!
cd ..

# Wait for frontend
echo "⏳ Waiting for frontend server..."
sleep 5
if curl -s http://localhost:5173 > /dev/null; then
    echo "✅ Frontend server is running on http://localhost:5173"
else
    echo "❌ Frontend server failed to start"
    kill $BACKEND_PID 2>/dev/null
    kill $FRONTEND_PID 2>/dev/null
    exit 1
fi

echo ""
echo "🎉 Helix.AI Development Environment is ready!"
echo "================================================"
echo "🌐 Frontend: http://localhost:5173"
echo "🔧 Backend:  http://localhost:8001"
echo "📚 API Docs: http://localhost:8001/docs"
echo "🧬 Plasmid Visualizer: Available in the UI"
echo ""
echo "Press Ctrl+C to stop all servers"

# Cleanup function
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

# Wait for user to stop
wait 