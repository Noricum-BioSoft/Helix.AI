#!/usr/bin/env bash
set -euo pipefail

# Script to create IAM user and access keys for GitHub Actions
# This is a quick alternative to OIDC setup

# Colors
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
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

USER_NAME="github-actions-deploy"
POLICY_NAME="GitHubActionsDeployPolicy"

print_info "Setting up IAM user and access keys for GitHub Actions..."
print_warning "This is less secure than OIDC. For production, use OIDC instead."
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

# Create policy JSON
print_info "Creating IAM policy..."
POLICY_DOC=$(cat <<'EOF'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "ecr:GetAuthorizationToken",
        "ecr:BatchCheckLayerAvailability",
        "ecr:GetDownloadUrlForLayer",
        "ecr:BatchGetImage",
        "ecr:PutImage",
        "ecr:InitiateLayerUpload",
        "ecr:UploadLayerPart",
        "ecr:CompleteLayerUpload",
        "ecr:DescribeRepositories",
        "ecr:DescribeImages"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "s3:PutObject",
        "s3:GetObject",
        "s3:DeleteObject",
        "s3:ListBucket"
      ],
      "Resource": [
        "arn:aws:s3:::*",
        "arn:aws:s3:::*/*"
      ]
    },
    {
      "Effect": "Allow",
      "Action": [
        "cloudfront:CreateInvalidation",
        "cloudfront:GetInvalidation",
        "cloudfront:ListInvalidations"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "ecs:DescribeTaskDefinition",
        "ecs:RegisterTaskDefinition",
        "ecs:UpdateService",
        "ecs:DescribeServices",
        "ecs:DescribeTasks",
        "ecs:ListTasks"
      ],
      "Resource": "*"
    },
    {
      "Effect": "Allow",
      "Action": [
        "logs:CreateLogGroup",
        "logs:CreateLogStream",
        "logs:PutLogEvents",
        "logs:DescribeLogStreams"
      ],
      "Resource": "*"
    }
  ]
}
EOF
)

# Check if policy exists
POLICY_ARN=$(aws iam list-policies --scope Local --query "Policies[?PolicyName=='$POLICY_NAME'].Arn" --output text 2>/dev/null || echo "")

if [ -z "$POLICY_ARN" ]; then
    print_info "Creating policy: $POLICY_NAME"
    POLICY_ARN=$(aws iam create-policy \
        --policy-name "$POLICY_NAME" \
        --policy-document "$POLICY_DOC" \
        --query 'Policy.Arn' \
        --output text)
    print_success "Policy created: $POLICY_ARN"
else
    print_info "Policy already exists: $POLICY_ARN"
    print_info "Updating policy..."
    aws iam create-policy-version \
        --policy-arn "$POLICY_ARN" \
        --policy-document "$POLICY_DOC" \
        --set-as-default > /dev/null
    print_success "Policy updated"
fi

# Attach policy to user
print_info "Attaching policy to user..."
aws iam attach-user-policy \
    --user-name "$USER_NAME" \
    --policy-arn "$POLICY_ARN" \
    > /dev/null 2>&1 || true
print_success "Policy attached"

# Create access keys
print_info "Creating access keys..."
KEYS=$(aws iam create-access-key --user-name "$USER_NAME" --output json)
ACCESS_KEY_ID=$(echo "$KEYS" | jq -r '.AccessKey.AccessKeyId')
SECRET_ACCESS_KEY=$(echo "$KEYS" | jq -r '.AccessKey.SecretAccessKey')

print_success "Access keys created!"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Add these to your GitHub repository secrets:"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
echo "  AWS_ACCESS_KEY_ID = $ACCESS_KEY_ID"
echo "  AWS_SECRET_ACCESS_KEY = $SECRET_ACCESS_KEY"
echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo ""
print_warning "Important: Save these keys now - you won't be able to see the secret key again!"
echo ""
echo "To update the workflow to use access keys instead of OIDC:"
echo "  1. Copy .github/workflows/deploy.yml.access-keys to .github/workflows/deploy.yml"
echo "  2. Add the secrets above to GitHub"
echo ""
echo "Or see: docs/GITHUB_OIDC_SETUP.md for OIDC setup instructions"

