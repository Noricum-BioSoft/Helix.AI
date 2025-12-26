#!/usr/bin/env bash
set -euo pipefail

# Script to help set up DNS records for Helix.AI custom domain
# This script outputs the DNS records you need to add to your DNS provider

# Usage:
#   ./scripts/aws/setup-dns-records.sh <STACK_NAME> [REGION] [DOMAIN]
#
# Example:
#   ./scripts/aws/setup-dns-records.sh HelixAIStack us-west-1 helix-ai.noricum-biosoft.com

STACK_NAME="${1:-HelixAIStack}"
REGION="${2:-us-west-1}"
DOMAIN="${3:-helix-ai.noricum-biosoft.com}"

echo "=========================================="
echo "Helix.AI DNS Records Setup"
echo "=========================================="
echo ""
echo "Stack: ${STACK_NAME}"
echo "Region: ${REGION}"
echo "Domain: ${DOMAIN}"
echo ""

# Get ALB DNS name
echo "Fetching ALB DNS name..."
ALB_DNS=$(aws cloudformation describe-stacks \
  --stack-name "${STACK_NAME}" \
  --region "${REGION}" \
  --query 'Stacks[0].Outputs[?OutputKey==`ALBDNSName`].OutputValue' \
  --output text 2>/dev/null || echo "")

if [ -z "${ALB_DNS}" ]; then
  echo "❌ Error: Could not find ALB DNS name. Is the stack deployed?"
  exit 1
fi

# Get CloudFront domain
echo "Fetching CloudFront domain..."
CF_DOMAIN=$(aws cloudformation describe-stacks \
  --stack-name "${STACK_NAME}" \
  --region "${REGION}" \
  --query 'Stacks[0].Outputs[?OutputKey==`CloudFrontDomainName`].OutputValue' \
  --output text 2>/dev/null || echo "")

if [ -z "${CF_DOMAIN}" ]; then
  echo "❌ Error: Could not find CloudFront domain. Is the stack deployed?"
  exit 1
fi

# Extract subdomain (e.g., "helix-ai" from "helix-ai.noricum-biosoft.com")
SUBDOMAIN=$(echo "${DOMAIN}" | cut -d'.' -f1)
PARENT_DOMAIN=$(echo "${DOMAIN}" | cut -d'.' -f2-)

echo ""
echo "=========================================="
echo "DNS Records to Add"
echo "=========================================="
echo ""
echo "Add these records to your DNS provider for: ${PARENT_DOMAIN}"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "1. Frontend (CloudFront) - CNAME Record"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "Type:    CNAME"
echo "Name:    ${SUBDOMAIN}"
echo "Value:   ${CF_DOMAIN}"
echo "TTL:     300 (or your provider's default)"
echo ""
echo "This will make ${DOMAIN} point to your CloudFront distribution."
echo ""

# Check if custom API domain is configured
API_DOMAIN=$(aws cloudformation describe-stacks \
  --stack-name "${STACK_NAME}" \
  --region "${REGION}" \
  --query 'Stacks[0].Outputs[?OutputKey==`BackendCustomDomain`].OutputValue' \
  --output text 2>/dev/null || echo "")

if [ -n "${API_DOMAIN}" ] && [ "${API_DOMAIN}" != "None" ]; then
  API_SUBDOMAIN=$(echo "${API_DOMAIN}" | cut -d'.' -f1)
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "2. Backend API (ALB) - CNAME or ALIAS Record"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo ""
  echo "Type:    CNAME (or ALIAS if your provider supports it)"
  echo "Name:    ${API_SUBDOMAIN}"
  echo "Value:   ${ALB_DNS}"
  echo "TTL:     300 (or your provider's default)"
  echo ""
  echo "This will make ${API_DOMAIN} point to your Application Load Balancer."
  echo ""
  echo "Note: If your DNS provider supports ALIAS records (Route 53, Cloudflare, etc.),"
  echo "      use ALIAS instead of CNAME for better performance."
else
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "2. Backend API (ALB) - Optional CNAME Record"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo ""
  echo "If you want a custom domain for the API (e.g., api.${DOMAIN}), add:"
  echo ""
  echo "Type:    CNAME (or ALIAS if your provider supports it)"
  echo "Name:    api (or your preferred subdomain)"
  echo "Value:   ${ALB_DNS}"
  echo "TTL:     300"
  echo ""
  echo "Then update HELIX_API_DOMAIN in your CDK deployment."
  echo ""
fi

echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "Next Steps"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "1. Add the DNS records above to your DNS provider"
echo "2. Wait for DNS propagation (usually 5-30 minutes, can take up to 48 hours)"
echo "3. Test DNS resolution:"
echo "   dig ${DOMAIN}"
echo "   nslookup ${DOMAIN}"
echo "4. Test the endpoints:"
echo "   curl -I https://${DOMAIN}"
if [ -n "${API_DOMAIN}" ] && [ "${API_DOMAIN}" != "None" ]; then
  echo "   curl https://${API_DOMAIN}/health"
fi
echo ""
echo "For more information, see: docs/CUSTOM_DOMAIN_SETUP.md"
echo ""


