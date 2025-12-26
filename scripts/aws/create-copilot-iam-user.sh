#!/usr/bin/env bash
set -euo pipefail

# Script to create IAM user for AWS Copilot
# This fixes the "Roles may not be assumed by root accounts" error

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

print_info() {
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

USER_NAME="${1:-copilot-admin}"
POLICY_ARN="arn:aws:iam::aws:policy/AdministratorAccess"

print_info "Creating IAM user for AWS Copilot..."
print_warning "This will create a user with AdministratorAccess policy"
echo ""

# Check if user already exists
if aws iam get-user --user-name "$USER_NAME" &>/dev/null; then
    print_warning "IAM user '$USER_NAME' already exists"
    read -p "Continue anyway? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
else
    print_info "Creating IAM user: $USER_NAME"
    aws iam create-user --user-name "$USER_NAME"
    print_success "User created"
fi

# Attach administrator policy
print_info "Attaching AdministratorAccess policy..."
if aws iam list-attached-user-policies --user-name "$USER_NAME" --query "AttachedPolicies[?PolicyArn=='$POLICY_ARN']" --output text | grep -q "$POLICY_ARN"; then
    print_warning "Policy already attached"
else
    aws iam attach-user-policy \
        --user-name "$USER_NAME" \
        --policy-arn "$POLICY_ARN"
    print_success "Policy attached"
fi

# Create access key
print_info "Creating access key..."
CREDS_FILE="/tmp/copilot-${USER_NAME}-credentials.json"
aws iam create-access-key --user-name "$USER_NAME" \
    --output json > "$CREDS_FILE"

ACCESS_KEY=$(jq -r '.AccessKey.AccessKeyId' "$CREDS_FILE")
SECRET_KEY=$(jq -r '.AccessKey.SecretAccessKey' "$CREDS_FILE")

print_success "Access key created"
echo ""
print_info "Credentials saved to: $CREDS_FILE"
print_warning "⚠️  Keep this file secure! Delete it after configuring AWS CLI."
echo ""

# Show next steps
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
print_success "IAM User Setup Complete!"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
print_info "Next steps:"
echo ""
echo "1. Configure AWS CLI with new credentials:"
echo "   ${BLUE}aws configure --profile copilot-admin${NC}"
echo ""
echo "   When prompted, enter:"
echo "   - AWS Access Key ID: ${ACCESS_KEY}"
echo "   - AWS Secret Access Key: ${SECRET_KEY}"
echo "   - Default region: us-west-1"
echo "   - Default output: json"
echo ""
echo "2. Verify credentials:"
echo "   ${BLUE}aws sts get-caller-identity --profile copilot-admin${NC}"
echo ""
echo "3. Use Copilot with new profile:"
echo "   ${BLUE}export AWS_PROFILE=copilot-admin${NC}"
echo "   ${BLUE}copilot env deploy --name production${NC}"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
print_warning "Security reminder:"
echo "  - Delete $CREDS_FILE after configuring AWS CLI"
echo "  - Rotate access keys every 90 days"
echo "  - Never commit credentials to git"
echo ""


