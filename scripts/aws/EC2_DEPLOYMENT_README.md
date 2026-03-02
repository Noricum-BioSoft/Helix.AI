# EC2 Deployment Quick Start

This directory contains scripts for deploying Helix.AI to AWS using EC2 instances for the backend.

## Quick Start

1. **Set up configuration**:
   ```bash
   ./setup-deployment.sh
   # Or copy deploy.config.example to deploy.config and edit it
   ```

2. **Deploy everything**:
   ```bash
   ./deploy-all.sh
   ```

This will deploy:
- Infrastructure (CDK stack)
- Backend Docker image to ECR
- EC2 instance with backend container
- Frontend to S3/CloudFront
- Configure ALB to use EC2

## Scripts Overview

### `deploy-all.sh`
Complete deployment script that does everything in one go.

**Usage:**
```bash
./deploy-all.sh [deploy.config]
```

### `deploy-to-ec2.sh`
Deploys the backend to an EC2 instance.

**Usage:**
```bash
./deploy-to-ec2.sh [deploy.config]
```

This script:
- Builds and pushes Docker image to ECR
- Waits for EC2 instance to be ready
- Provides instructions for updating the container

### `switch-alb-to-ec2.sh`
Switches the ALB listener to use the EC2 target group instead of ECS.

**Usage:**
```bash
./switch-alb-to-ec2.sh [deploy.config]
```

### `deploy.sh`
Original deployment script for ECS/Fargate deployments. Still works for frontend deployment.

**Usage:**
```bash
./deploy.sh [deploy.config]
```

## Configuration

Create a `deploy.config` file (copy from `deploy.config.example`) with your AWS details:

```bash
# Required
AWS_REGION=us-east-1
AWS_ACCOUNT_ID=123456789012
STACK_NAME=HelixAIStack
ECR_REPOSITORY=helix-ai-backend
S3_BUCKET=helix-ai-frontend-123456789012-us-east-1

# Optional
CLOUDFRONT_DISTRIBUTION_ID=
VITE_API_BASE_URL=http://your-alb-url.elb.amazonaws.com
```

## Prerequisites

Before deploying:

1. **Set up API keys in Secrets Manager**:
   ```bash
   aws secretsmanager create-secret \
       --name helix-ai-production-OPENAI_API_KEY \
       --secret-string "your-key-here" \
       --region us-east-1
   ```

2. **Ensure AWS credentials are configured**:
   ```bash
   aws configure
   ```

3. **Bootstrap CDK** (first time only):
   ```bash
   cd ../../infrastructure
   cdk bootstrap
   ```

## Step-by-Step Deployment

If you prefer manual steps:

### 1. Deploy Infrastructure
```bash
cd ../../infrastructure
cdk deploy
```

### 2. Deploy Backend to EC2
```bash
cd ../../scripts/aws
./deploy-to-ec2.sh
```

### 3. Switch ALB to EC2
```bash
./switch-alb-to-ec2.sh
```

### 4. Deploy Frontend
```bash
./deploy.sh
```

## Viewing Stack Outputs

After deployment, view stack outputs:

```bash
aws cloudformation describe-stacks \
    --stack-name HelixAIStack \
    --query 'Stacks[0].Outputs'
```

Important outputs:
- `BackendAPIURL`: Backend API URL
- `FrontendURL`: Frontend CloudFront URL  
- `EC2InstanceId`: EC2 instance ID
- `EC2TargetGroupARN`: ALB target group ARN

## Troubleshooting

### EC2 Instance Not Starting
Check user data logs:
```bash
aws ssm send-command \
    --instance-ids <INSTANCE_ID> \
    --document-name "AWS-RunShellScript" \
    --parameters 'commands=["cat /var/log/user-data.log"]'
```

### Backend Not Responding
Check container logs:
```bash
aws ssm send-command \
    --instance-ids <INSTANCE_ID> \
    --document-name "AWS-RunShellScript" \
    --parameters 'commands=["docker logs helix-ai-backend"]'
```

### ALB Health Checks Failing
Check target group health:
```bash
aws elbv2 describe-target-health \
    --target-group-arn <TARGET_GROUP_ARN>
```

## More Information

See the [complete deployment guide](../../docs/deployment/AWS_EC2_DEPLOYMENT_GUIDE.md) for detailed information.
