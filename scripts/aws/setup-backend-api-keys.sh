#!/bin/bash
# Script to set up API keys for CDK-deployed backend in AWS Secrets Manager

set -euo pipefail

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

export AWS_REGION="${AWS_REGION:-us-west-1}"

print_info "Setting up API keys for CDK-deployed backend..."
echo ""

# Secret names (matching CDK stack)
OPENAI_SECRET_NAME="helix-ai-backend-OPENAI_API_KEY"
DEEPSEEK_SECRET_NAME="helix-ai-backend-DEEPSEEK_API_KEY"

# OpenAI Key
if aws secretsmanager describe-secret --secret-id "$OPENAI_SECRET_NAME" --region "$AWS_REGION" &>/dev/null; then
    print_warning "OpenAI secret already exists: $OPENAI_SECRET_NAME"
    read -p "Update it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        read -sp "Enter your OpenAI API key: " OPENAI_KEY
        echo
        aws secretsmanager update-secret \
            --secret-id "$OPENAI_SECRET_NAME" \
            --secret-string "$OPENAI_KEY" \
            --region "$AWS_REGION" > /dev/null
        print_success "OpenAI key updated"
    fi
else
    read -sp "Enter your OpenAI API key (or press Enter to skip): " OPENAI_KEY
    echo
    if [[ -n "$OPENAI_KEY" ]]; then
        aws secretsmanager create-secret \
            --name "$OPENAI_SECRET_NAME" \
            --secret-string "$OPENAI_KEY" \
            --region "$AWS_REGION" > /dev/null
        print_success "OpenAI key stored in Secrets Manager"
    fi
fi

# DeepSeek Key
if aws secretsmanager describe-secret --secret-id "$DEEPSEEK_SECRET_NAME" --region "$AWS_REGION" &>/dev/null; then
    print_warning "DeepSeek secret already exists: $DEEPSEEK_SECRET_NAME"
    read -p "Update it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        read -sp "Enter your DeepSeek API key: " DEEPSEEK_KEY
        echo
        aws secretsmanager update-secret \
            --secret-id "$DEEPSEEK_SECRET_NAME" \
            --secret-string "$DEEPSEEK_KEY" \
            --region "$AWS_REGION" > /dev/null
        print_success "DeepSeek key updated"
    fi
else
    read -sp "Enter your DeepSeek API key (or press Enter to skip): " DEEPSEEK_KEY
    echo
    if [[ -n "$DEEPSEEK_KEY" ]]; then
        aws secretsmanager create-secret \
            --name "$DEEPSEEK_SECRET_NAME" \
            --secret-string "$DEEPSEEK_KEY" \
            --region "$AWS_REGION" > /dev/null
        print_success "DeepSeek key stored in Secrets Manager"
    fi
fi

echo ""
print_success "API keys configured!"
echo ""
print_info "Next steps:"
echo "  1. Update CDK stack to reference these secrets (already done in helix_stack.py)"
echo "  2. Deploy the updated stack: cd infrastructure && cdk deploy"
echo "  3. Force ECS service to redeploy: aws ecs update-service --cluster helix-ai-cluster --service HelixAIStack-ServiceD69D759B-G23PygXszFhF --force-new-deployment --region us-west-1"


