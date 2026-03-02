#!/bin/bash
# User data script for EC2 instance running Helix.AI backend
set -e
exec > >(tee /var/log/user-data.log|logger -t user-data -s 2>/dev/console) 2>&1

# Update system
dnf update -y

# Install Docker and AWS CLI
dnf install -y docker awscli curl jq

# Start and enable Docker
systemctl start docker
systemctl enable docker

# Add ec2-user to docker group
usermod -aG docker ec2-user

# Install Docker Compose
curl -L "https://github.com/docker/compose/releases/latest/download/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
chmod +x /usr/local/bin/docker-compose

# Configure AWS CLI region (will be replaced by CDK)
AWS_REGION_PLACEHOLDER

# Login to ECR (will be replaced by CDK)
ECR_LOGIN_PLACEHOLDER

# Get secrets from Secrets Manager
OPENAI_KEY=$(aws secretsmanager get-secret-value --secret-id helix-ai-production-OPENAI_API_KEY --region ${AWS_REGION} --query SecretString --output text 2>/dev/null || echo '')
DEEPSEEK_KEY=$(aws secretsmanager get-secret-value --secret-id helix-ai-production-DEEPSEEK_API_KEY --region ${AWS_REGION} --query SecretString --output text 2>/dev/null || echo '')

# Create environment file
mkdir -p /opt/helix-ai
cat > /opt/helix-ai/.env <<EOF
OPENAI_API_KEY=${OPENAI_KEY}
DEEPSEEK_API_KEY=${DEEPSEEK_KEY}
PYTHONUNBUFFERED=1
EOF

# Pull and run the Docker image (will be replaced by CDK)
DOCKER_PULL_PLACEHOLDER

# Create docker-compose.yml for easy management (will be replaced by CDK)
DOCKER_COMPOSE_PLACEHOLDER

# Run the container
cd /opt/helix-ai
docker-compose up -d

# Set up log rotation for Docker
cat > /etc/logrotate.d/docker-containers <<EOF
/var/lib/docker/containers/*/*.log {
    rotate 7
    daily
    compress
    size=1M
    missingok
    delaycompress
    copytruncate
}
EOF

# Wait for container to be healthy
echo "Waiting for backend to be healthy..."
for i in {1..30}; do
    if docker exec helix-ai-backend curl -f http://localhost:8001/health >/dev/null 2>&1; then
        echo "Backend is healthy!"
        break
    fi
    echo "Waiting... ($i/30)"
    sleep 5
done
