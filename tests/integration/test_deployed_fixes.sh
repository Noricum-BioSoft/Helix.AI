#!/bin/bash
# Test script to verify the deployed fixes are working

set -e

BACKEND_URL="${BACKEND_URL:-http://HelixA-ALBAE-D7BksiQIynZb-1051248867.us-west-1.elb.amazonaws.com}"
# Alternative: use CloudFront URL
# BACKEND_URL="https://d2a8mt5n89vos4.cloudfront.net"

echo "🧪 Testing Deployed Fixes"
echo "=========================="
echo "Backend URL: $BACKEND_URL"
echo ""

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test 1: Health check
echo "1️⃣ Testing health endpoint..."
HEALTH_RESPONSE=$(curl -s -w "\n%{http_code}" --max-time 10 "${BACKEND_URL}/health" || echo "FAILED")
HTTP_CODE=$(echo "$HEALTH_RESPONSE" | tail -n1)
BODY=$(echo "$HEALTH_RESPONSE" | sed '$d')

if [ "$HTTP_CODE" = "200" ]; then
    echo -e "${GREEN}✅ Health check passed (HTTP $HTTP_CODE)${NC}"
    echo "   Response: $BODY"
else
    echo -e "${RED}❌ Health check failed (HTTP $HTTP_CODE)${NC}"
    exit 1
fi
echo ""

# Test 2: Test sequence_alignment tool (Bio.Align.Applications fix)
echo "2️⃣ Testing sequence_alignment tool (Bio.Align.Applications fix)..."
ALIGNMENT_PAYLOAD=$(cat <<EOF
{
  "command": "Align these sequences: >seq1 ATGCGATCG >seq2 ATGCGATCA",
  "session_id": "test-session-$(date +%s)"
}
EOF
)

ALIGNMENT_RESPONSE=$(curl -s -w "\n%{http_code}" -X POST \
    -H "Content-Type: application/json" \
    --max-time 30 \
    -d "$ALIGNMENT_PAYLOAD" \
    "${BACKEND_URL}/execute" || echo "FAILED")

ALIGNMENT_HTTP_CODE=$(echo "$ALIGNMENT_RESPONSE" | tail -n1)
ALIGNMENT_BODY=$(echo "$ALIGNMENT_RESPONSE" | sed '$d')

if [ "$ALIGNMENT_HTTP_CODE" = "200" ]; then
    # Check if response contains error about Bio.Align.Applications
    if echo "$ALIGNMENT_BODY" | grep -q "Bio.Align.Applications\|No module named"; then
        echo -e "${RED}❌ sequence_alignment failed - still has import error${NC}"
        echo "   Response: $ALIGNMENT_BODY"
    else
        echo -e "${GREEN}✅ sequence_alignment tool works (HTTP $ALIGNMENT_HTTP_CODE)${NC}"
        echo "   Response preview: $(echo "$ALIGNMENT_BODY" | head -c 200)..."
    fi
else
    echo -e "${YELLOW}⚠️  sequence_alignment returned HTTP $ALIGNMENT_HTTP_CODE${NC}"
    echo "   Response: $ALIGNMENT_BODY"
fi
echo ""

# Test 3: Test mutate_sequence tool (langchain_core.tools fix)
echo "3️⃣ Testing mutate_sequence tool (langchain_core.tools fix)..."
MUTATION_PAYLOAD=$(cat <<EOF
{
  "command": "Mutate sequence ATGCGATCG to create 3 variants",
  "session_id": "test-session-$(date +%s)"
}
EOF
)

MUTATION_RESPONSE=$(curl -s -w "\n%{http_code}" -X POST \
    -H "Content-Type: application/json" \
    --max-time 30 \
    -d "$MUTATION_PAYLOAD" \
    "${BACKEND_URL}/execute" || echo "FAILED")

MUTATION_HTTP_CODE=$(echo "$MUTATION_RESPONSE" | tail -n1)
MUTATION_BODY=$(echo "$MUTATION_RESPONSE" | sed '$d')

if [ "$MUTATION_HTTP_CODE" = "200" ]; then
    # Check if response contains error about langchain.agents
    if echo "$MUTATION_BODY" | grep -q "langchain.agents\|cannot import name 'tool'"; then
        echo -e "${RED}❌ mutate_sequence failed - still has import error${NC}"
        echo "   Response: $MUTATION_BODY"
    else
        echo -e "${GREEN}✅ mutate_sequence tool works (HTTP $MUTATION_HTTP_CODE)${NC}"
        echo "   Response preview: $(echo "$MUTATION_BODY" | head -c 200)..."
    fi
else
    echo -e "${YELLOW}⚠️  mutate_sequence returned HTTP $MUTATION_HTTP_CODE${NC}"
    echo "   Response: $MUTATION_BODY"
fi
echo ""

# Test 4: List available tools via MCP endpoint
echo "4️⃣ Testing tool listing..."
TOOLS_RESPONSE=$(curl -s -w "\n%{http_code}" --max-time 10 "${BACKEND_URL}/mcp/tools" || echo "FAILED")
TOOLS_HTTP_CODE=$(echo "$TOOLS_RESPONSE" | tail -n1)
TOOLS_BODY=$(echo "$TOOLS_RESPONSE" | sed '$d')

if [ "$TOOLS_HTTP_CODE" = "200" ]; then
    echo -e "${GREEN}✅ Tools endpoint works (HTTP $TOOLS_HTTP_CODE)${NC}"
    # Check if both tools are listed
    if echo "$TOOLS_BODY" | grep -q "sequence_alignment" && echo "$TOOLS_BODY" | grep -q "mutate_sequence"; then
        echo -e "${GREEN}✅ Both sequence_alignment and mutate_sequence are available${NC}"
    else
        echo -e "${YELLOW}⚠️  Some tools may be missing from the list${NC}"
    fi
else
    echo -e "${YELLOW}⚠️  Tools endpoint returned HTTP $TOOLS_HTTP_CODE${NC}"
fi
echo ""

echo "=========================="
echo "✅ Testing complete!"
echo ""
echo "If all tests passed, your fixes are working correctly! 🎉"

