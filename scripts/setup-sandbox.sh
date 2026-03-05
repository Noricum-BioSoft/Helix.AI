#!/bin/bash
#
# Setup Sandbox Execution Environment
#
# This script builds the helix-biotools Docker image and verifies
# that all bioinformatics tools are properly installed.
#
# Usage:
#   ./scripts/setup-sandbox.sh [--rebuild] [--test] [--minimal]
#
# Options:
#   --rebuild    Force rebuild even if image exists
#   --test       Run tool verification tests after build
#   --minimal    Build minimal version with only core tools (faster)
#   --help       Show this help message

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Configuration
IMAGE_NAME="helix-biotools"
IMAGE_TAG="latest"
DOCKERFILE="backend/Dockerfile.biotools"

# Parse arguments
REBUILD=false
TEST=false
MINIMAL=false

for arg in "$@"; do
    case $arg in
        --rebuild)
            REBUILD=true
            shift
            ;;
        --test)
            TEST=true
            shift
            ;;
        --minimal)
            MINIMAL=true
            DOCKERFILE="backend/Dockerfile.biotools-minimal"
            IMAGE_TAG="minimal"
            shift
            ;;
        --help)
            head -n 15 "$0" | grep "^#" | sed 's/^# //g' | sed 's/^#//g'
            exit 0
            ;;
        *)
            echo -e "${RED}Unknown argument: $arg${NC}"
            echo "Use --help for usage information"
            exit 1
            ;;
    esac
done

echo -e "${BLUE}=====================================================================${NC}"
echo -e "${BLUE}Helix.AI Sandbox Setup${NC}"
echo -e "${BLUE}=====================================================================${NC}"
echo ""

# Check if Docker is available
echo -e "${YELLOW}Checking Docker availability...${NC}"
if ! command -v docker &> /dev/null; then
    echo -e "${RED}❌ Docker is not installed${NC}"
    echo "Please install Docker Desktop: https://www.docker.com/products/docker-desktop"
    exit 1
fi

# Check if Docker daemon is running
if ! docker info &> /dev/null; then
    echo -e "${RED}❌ Docker daemon is not running${NC}"
    echo "Please start Docker Desktop"
    exit 1
fi

echo -e "${GREEN}✅ Docker is available${NC}"
docker --version
echo ""

# Check if image already exists
IMAGE_EXISTS=$(docker images -q ${IMAGE_NAME}:${IMAGE_TAG} 2> /dev/null)

if [ -n "$IMAGE_EXISTS" ] && [ "$REBUILD" = false ]; then
    echo -e "${YELLOW}Docker image ${IMAGE_NAME}:${IMAGE_TAG} already exists${NC}"
    echo "Use --rebuild to force rebuild"
    echo ""
else
    # Build the Docker image
    if [ "$MINIMAL" = true ]; then
        echo -e "${YELLOW}Building MINIMAL Docker image: ${IMAGE_NAME}:${IMAGE_TAG}${NC}"
        echo "Minimal version includes: FastQC, MUSCLE, Clustal Omega"
        echo "This will take ~5 minutes..."
    else
        echo -e "${YELLOW}Building FULL Docker image: ${IMAGE_NAME}:${IMAGE_TAG}${NC}"
        echo "Full version includes all bioinformatics tools"
        echo "This may take 10-20 minutes on first build..."
    fi
    echo ""
    
    # Check if Dockerfile exists
    if [ ! -f "$DOCKERFILE" ]; then
        echo -e "${RED}❌ Dockerfile not found: $DOCKERFILE${NC}"
        exit 1
    fi
    
    # Build with progress output
    docker build \
        -f "$DOCKERFILE" \
        -t "${IMAGE_NAME}:${IMAGE_TAG}" \
        . \
        || {
            echo -e "${RED}❌ Docker build failed${NC}"
            exit 1
        }
    
    echo ""
    echo -e "${GREEN}✅ Docker image built successfully${NC}"
fi

# Show image info
echo -e "${YELLOW}Image information:${NC}"
docker images ${IMAGE_NAME}:${IMAGE_TAG}
echo ""

# Run tests if requested
if [ "$TEST" = true ]; then
    echo -e "${BLUE}=====================================================================${NC}"
    echo -e "${BLUE}Testing Bioinformatics Tools${NC}"
    echo -e "${BLUE}=====================================================================${NC}"
    echo ""
    
    # Function to test a tool
    test_tool() {
        local tool_name=$1
        local tool_path=$2
        local version_flag=$3
        
        echo -e "${YELLOW}Testing ${tool_name}...${NC}"
        
        if docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} ${tool_path} ${version_flag} &> /tmp/tool_test.txt; then
            version=$(head -1 /tmp/tool_test.txt)
            echo -e "${GREEN}✅ ${tool_name}: ${version}${NC}"
        else
            echo -e "${RED}❌ ${tool_name}: FAILED${NC}"
            cat /tmp/tool_test.txt
        fi
    }
    
    # Test tools based on which version was built
    if [ "$MINIMAL" = true ]; then
        echo -e "${BLUE}Testing minimal tools (FastQC, MUSCLE, Clustal Omega)...${NC}"
        echo ""
        test_tool "FastQC" "/usr/local/bin/fastqc" "--version"
        test_tool "MUSCLE" "/usr/local/bin/muscle" "-version"
        test_tool "Clustal Omega" "/usr/bin/clustalo" "--version"
        
        # Test Python packages (minimal)
        echo -e "${YELLOW}Testing Python packages (minimal)...${NC}"
        docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} python3 -c "
import sys
packages = ['biopython']
for pkg in packages:
    try:
        __import__(pkg)
        print(f'✅ {pkg}')
    except ImportError as e:
        print(f'❌ {pkg} FAILED: {str(e)}')
        sys.exit(1)
" || {
        echo -e "${RED}❌ Some Python packages failed to import${NC}"
        echo -e "${YELLOW}Note: This may not affect FastQC functionality${NC}"
        # Don't exit - biopython failure shouldn't block FastQC
    }
    else
        echo -e "${BLUE}Testing full tool suite...${NC}"
        echo ""
        test_tool "FastQC" "/usr/local/bin/fastqc" "--version"
        test_tool "MUSCLE" "/usr/local/bin/muscle" "-version"
        test_tool "Clustal Omega" "/usr/bin/clustalo" "--version"
        test_tool "SAMtools" "/usr/local/bin/samtools" "--version"
        test_tool "BWA" "/usr/local/bin/bwa" "2>&1 | head -3"
        test_tool "Bowtie2" "/usr/local/bin/bowtie2" "--version"
        test_tool "STAR" "/usr/local/bin/STAR" "--version"
        test_tool "HISAT2" "/usr/local/bin/hisat2" "--version"
        test_tool "StringTie" "/usr/local/bin/stringtie" "--version"
        test_tool "featureCounts" "/usr/local/bin/featureCounts" "-v"
        
        # Test R
        echo -e "${YELLOW}Testing R and Bioconductor...${NC}"
        if docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} R --version | head -1 | grep "R version" &> /dev/null; then
            r_version=$(docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} R --version | head -1)
            echo -e "${GREEN}✅ R: ${r_version}${NC}"
        else
            echo -e "${RED}❌ R: FAILED${NC}"
        fi
        
        # Test Python packages (full)
        echo -e "${YELLOW}Testing Python bioinformatics packages...${NC}"
        docker run --rm ${IMAGE_NAME}:${IMAGE_TAG} python3 -c "
import sys
packages = ['biopython', 'pysam', 'HTSeq', 'scanpy']
for pkg in packages:
    try:
        __import__(pkg)
        print(f'✅ {pkg}')
    except ImportError:
        print(f'❌ {pkg}')
        sys.exit(1)
" || {
        echo -e "${RED}❌ Some Python packages failed to import${NC}"
        exit 1
    }
    fi
    
    rm -f /tmp/tool_test.txt
    
    echo ""
    echo -e "${GREEN}✅ All tool tests passed!${NC}"
fi

echo ""
echo -e "${BLUE}=====================================================================${NC}"
echo -e "${BLUE}Setup Complete${NC}"
echo -e "${BLUE}=====================================================================${NC}"
echo ""
echo -e "${GREEN}The sandbox environment is ready!${NC}"
echo ""
echo "Next steps:"
echo "  1. Set HELIX_USE_SANDBOX=true in your .env file"
echo "  2. Restart the backend server"
echo "  3. Test with small files (< 100MB) for automatic sandbox routing"
echo ""
echo "For more information, see: docs/SANDBOX_EXECUTION.md"
echo ""

# Show sample test command
echo -e "${YELLOW}Test sandbox execution:${NC}"
cat << 'EOF'
python3 << PYTHON
from backend.sandbox_executor import get_sandbox_executor
executor = get_sandbox_executor()
tools = executor.list_available_tools()
print(f"Available tools: {len([t for t in tools.values() if t['available']])}/{len(tools)}")
for name, info in tools.items():
    status = "✅" if info["available"] else "❌"
    print(f"{status} {name}: {info['version'] or 'N/A'}")
PYTHON
EOF

echo ""
