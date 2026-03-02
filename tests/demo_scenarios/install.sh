#!/bin/bash
# Quick installation script with uv support
# Automatically detects and uses uv if available, falls back to pip

set -e  # Exit on error

cd "$(dirname "$0")"

echo "╔══════════════════════════════════════════════════════════════════════╗"
echo "║       Helix.AI Demo/Eval System - Dependency Installation           ║"
echo "╚══════════════════════════════════════════════════════════════════════╝"
echo ""

# Detect Python
if command -v python3 &> /dev/null; then
    PYTHON=python3
elif command -v python &> /dev/null; then
    PYTHON=python
else
    echo "❌ Python not found. Please install Python 3.10 or later."
    exit 1
fi

echo "🐍 Python: $($PYTHON --version)"
echo ""

# Check if uv is installed
if command -v uv &> /dev/null; then
    echo "🚀 Using uv (fast installer)"
    UV_AVAILABLE=true
    echo "   uv version: $(uv --version)"
else
    echo "📦 Using pip (standard installer)"
    UV_AVAILABLE=false
    echo "   💡 Tip: Install uv for 10-100x faster installations!"
    echo "   curl -LsSf https://astral.sh/uv/install.sh | sh"
fi
echo ""

# Parse arguments
INSTALL_TYPE="all"
FORCE_PIP=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --mock-only)
            INSTALL_TYPE="mock"
            shift
            ;;
        --backend-only)
            INSTALL_TYPE="backend"
            shift
            ;;
        --dev)
            INSTALL_TYPE="dev"
            shift
            ;;
        --pip)
            FORCE_PIP=true
            shift
            ;;
        --help)
            echo "Usage: ./install.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --mock-only      Install only framework dependencies (mock mode)"
            echo "  --backend-only   Install only backend dependencies (real execution)"
            echo "  --dev            Install all dependencies including dev tools"
            echo "  --pip            Force use of pip even if uv is available"
            echo "  --help           Show this help message"
            echo ""
            echo "Examples:"
            echo "  ./install.sh                    # Install all dependencies (default)"
            echo "  ./install.sh --mock-only        # Install only for mock testing"
            echo "  ./install.sh --dev              # Install everything including dev tools"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Run './install.sh --help' for usage information"
            exit 1
            ;;
    esac
done

echo "📋 Installation type: $INSTALL_TYPE"
echo ""

# Install dependencies
if [[ "$UV_AVAILABLE" == "true" && "$FORCE_PIP" == "false" ]]; then
    echo "⏳ Installing dependencies with uv..."
    echo ""
    
    case $INSTALL_TYPE in
        mock)
            uv pip install -e ".[test]"
            ;;
        backend)
            uv pip install -e ".[backend]"
            ;;
        dev)
            uv pip install -e ".[all,dev]"
            ;;
        all)
            uv pip install -e ".[all]"
            ;;
    esac
else
    echo "⏳ Installing dependencies with pip..."
    echo ""
    
    case $INSTALL_TYPE in
        mock)
            $PYTHON -m pip install -e "."
            ;;
        backend)
            $PYTHON -m pip install -e ".[backend]"
            ;;
        dev)
            $PYTHON -m pip install -e ".[all,dev]"
            ;;
        all)
            $PYTHON -m pip install -e ".[all]"
            ;;
    esac
fi

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Verify installation
echo "✅ Verifying installation..."
echo ""

$PYTHON check_dependencies.py

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""

# Show next steps
case $INSTALL_TYPE in
    mock)
        echo "🎉 Mock mode dependencies installed!"
        echo ""
        echo "You can now run:"
        echo "  • pytest test_scenarios.py -v"
        echo "  • python demo_cli.py"
        echo ""
        echo "To install backend dependencies for real execution:"
        echo "  ./install.sh --backend-only"
        ;;
    backend)
        echo "🎉 Backend dependencies installed!"
        echo ""
        echo "You can now run:"
        echo "  • python test_real_execution.py --scenario real_fastqc_small"
        ;;
    *)
        echo "🎉 All dependencies installed!"
        echo ""
        echo "You can now run:"
        echo "  • pytest test_scenarios.py -v          # Run all tests"
        echo "  • python demo_cli.py                   # Interactive demo"
        echo "  • python test_real_execution.py --compare  # Infrastructure test"
        echo ""
        echo "For real execution (requires AWS + API keys):"
        echo "  • python test_real_execution.py --scenario real_fastqc_small"
        ;;
esac

echo ""
echo "📖 For more information, see:"
echo "   • README.md"
echo "   • HOW_TO_RUN.md"
echo "   • INSTALL_DEPENDENCIES.md"
echo ""
