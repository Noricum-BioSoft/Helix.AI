# AWS Deployment Guide for Helix.AI

This guide provides step-by-step instructions for deploying Helix.AI to AWS using automated deployment scripts and infrastructure as code.

## Table of Contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Quick Start](#quick-start)
4. [Infrastructure Setup](#infrastructure-setup)
5. [Configuration](#configuration)
6. [Deployment Methods](#deployment-methods)
7. [Post-Deployment](#post-deployment)
8. [Troubleshooting](#troubleshooting)
9. [CI/CD Setup](#cicd-setup)

## Overview

Helix.AI can be deployed to AWS using two main approaches:

1. **Infrastructure as Code (IaC)**: AWS CDK for provisioning all AWS resources
2. **Automated Deployment**: Shell scripts for building and deploying application code

The deployment architecture includes:
- **Backend**: FastAPI application running on ECS Fargate
- **Frontend**: React application hosted on S3 with CloudFront CDN
- **Load Balancer**: Application Load Balancer (ALB) for backend traffic
- **Container Registry**: ECR for storing Docker images

## Prerequisites

### AWS Account Setup

1. **AWS Account**: Active AWS account with appropriate permissions
2. **AWS CLI**: Installed and configured with credentials
   ```bash
   aws configure
   ```
3. **AWS CDK**: Installed globally
   ```bash
   npm install -g aws-cdk
   cdk --version
   ```
4. **Docker**: Installed and running (for building images)
5. **Node.js**: Version 18+ (for frontend builds and CDK)

### Local Setup

1. **Python 3.9+**: For backend and CDK infrastructure
2. **Git**: For cloning and version control
3. **jq**: For JSON processing in scripts (optional but helpful)
   ```bash
   # macOS
   brew install jq
   
   # Linux
   sudo apt-get install jq
   ```

### IAM Permissions

Your AWS credentials need permissions for:
- EC2 (VPC, Security Groups, Load Balancers)
- ECS (Cluster, Services, Task Definitions)
- ECR (Repository operations)
- S3 (Bucket creation and management)
- CloudFront (Distribution management)
- CloudFormation (Stack operations)
- IAM (Role creation)
- CloudWatch (Logs and monitoring)

**Recommended**: Use an IAM user with AdministratorAccess for initial setup, then restrict permissions for CI/CD.

## Quick Start

### Option 1: Automated Full Deployment

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd Helix.AI
   ```

2. **Set up infrastructure** (first time only):
   ```bash
   cd infrastructure
   pip install -r requirements.txt
   cdk bootstrap
   cdk deploy
   ```

3. **Configure deployment settings**:
   ```bash
   cd scripts/aws
   cp deploy.config.example deploy.config
   # Edit deploy.config with your AWS details
   ```

4. **Deploy application**:
   ```bash
   ./deploy.sh
   ```

### Option 2: Manual Step-by-Step

See [Deployment Methods](#deployment-methods) section below.

## Infrastructure Setup

### Step 1: Deploy AWS Infrastructure

The CDK stack creates all necessary AWS resources:

```bash
cd infrastructure

# Install dependencies
pip install -r requirements.txt

# Bootstrap CDK (first time only in a region)
cdk bootstrap

# Preview what will be created
cdk synth

# Deploy infrastructure
cdk deploy
```

**Note**: The first deployment takes 10-15 minutes. Subsequent deployments are faster.

### Step 2: Capture Output Values

After deployment, note these values from the stack outputs:

- `ECRRepositoryURI`: ECR repository URI
- `ECSClusterName`: ECS cluster name
- `ECSServiceName`: ECS service name
- `ECSTaskDefinitionFamily`: Task definition family name
- `ALBDNSName`: Load balancer DNS name
- `FrontendBucketName`: S3 bucket name
- `CloudFrontDistributionId`: CloudFront distribution ID
- `CloudFrontDomainName`: CloudFront domain name

You can view outputs with:
```bash
aws cloudformation describe-stacks \
  --stack-name HelixAIStack \
  --query 'Stacks[0].Outputs' \
  --output table
```

## Configuration

### Deployment Configuration File

Create `scripts/aws/deploy.config`:

```bash
cd scripts/aws
cp deploy.config.example deploy.config
```

Edit `deploy.config` with your values:

```bash
# AWS Configuration
AWS_REGION=us-east-1
AWS_ACCOUNT_ID=123456789012

# ECR Configuration
ECR_REPOSITORY=helix-backend
IMAGE_TAG=latest

# S3 Configuration
S3_BUCKET=helix-frontend-bucket-123456789012-us-east-1

# CloudFront Configuration (optional)
CLOUDFRONT_DISTRIBUTION_ID=E1234567890ABC

# ECS Configuration (for automatic service updates)
ECS_CLUSTER_NAME=helix-ai-cluster
ECS_SERVICE_NAME=HelixStack-Service-XXXXXXXXX
ECS_TASK_DEFINITION_FAMILY=HelixStack-TaskDefinition-XXXXXXXXX

# Frontend Configuration
VITE_API_BASE_URL=http://HelixStack-ALB-XXXXXXXXX.us-east-1.elb.amazonaws.com

# Deployment Options
SKIP_BACKEND=false
SKIP_FRONTEND=false
SKIP_CLOUDFRONT_INVALIDATION=false
SKIP_ECS_UPDATE=false
```

### Frontend Environment Variables

The frontend needs the backend API URL at build time. Set `VITE_API_BASE_URL` in `deploy.config`.

**Easiest way to get the URL**: Use the helper script:
```bash
cd scripts/aws
./get-alb-url.sh
```

This will:
- Fetch the ALB DNS name from your CloudFormation stack
- Test if the backend is accessible
- Show you the exact URL to use

**Manual method**: Get it from stack outputs:
```bash
aws cloudformation describe-stacks \
  --stack-name HelixAIStack \
  --query 'Stacks[0].Outputs[?OutputKey==`ALBDNSName`].OutputValue' \
  --output text
```

Then add to `deploy.config`:
```bash
VITE_API_BASE_URL=http://your-alb-dns-name.region.elb.amazonaws.com
```

For production with HTTPS:
```bash
VITE_API_BASE_URL=https://api.yourdomain.com
```

**See also**: [Backend API URL Setup Guide](BACKEND_API_URL_SETUP.md) for detailed instructions.

## Deployment Methods

### Method 1: Master Deployment Script (Recommended)

The master script automates everything:

```bash
cd scripts/aws
./deploy.sh
```

Or with a custom config file:
```bash
./deploy.sh my-custom-config
```

The script will:
1. Validate prerequisites
2. Build and push backend Docker image to ECR
3. Update ECS service (if configured)
4. Build and deploy frontend to S3
5. Invalidate CloudFront cache (if configured)

### Method 2: Individual Scripts

Deploy components separately:

**Backend only**:
```bash
./scripts/aws/ecr_push_backend.sh \
  123456789012 \
  us-east-1 \
  helix-backend \
  latest
```

**Frontend only**:
```bash
./scripts/aws/s3_sync_frontend.sh \
  us-east-1 \
  helix-frontend-bucket \
  E1234567890ABC
```

### Method 3: Manual AWS CLI Commands

See individual script files for manual commands.

## Post-Deployment

### Verify Backend Deployment

1. **Check ECS service status**:
   ```bash
   aws ecs describe-services \
     --cluster helix-ai-cluster \
     --services <service-name> \
     --query 'services[0].{Status:status,Running:runningCount,Desired:desiredCount}'
   ```

2. **Test health endpoint**:
   ```bash
   curl http://<alb-dns-name>/health
   ```

3. **Check container logs**:
   ```bash
   # Get log stream name
   aws logs describe-log-streams \
     --log-group-name /ecs/helix-ai \
     --order-by LastEventTime \
     --descending \
     --max-items 1
   
   # View logs
   aws logs get-log-events \
     --log-group-name /ecs/helix-ai \
     --log-stream-name <stream-name>
   ```

### Verify Frontend Deployment

1. **Check S3 bucket**:
   ```bash
   aws s3 ls s3://<bucket-name>/
   ```

2. **Test CloudFront distribution**:
   ```bash
   curl -I https://<cloudfront-domain>/index.html
   ```

3. **Check CloudFront cache invalidation status**:
   ```bash
   aws cloudfront list-invalidations \
     --distribution-id <distribution-id>
   ```

### Update Environment Variables

If you need to update environment variables for the ECS service:

1. **Update task definition** (via AWS Console or CLI)
2. **Or update via CDK** and redeploy:
   ```bash
   cd infrastructure
   # Edit helix_infrastructure/helix_stack.py
   cdk deploy
   ```

### Scaling

**Manual scaling**:
```bash
aws ecs update-service \
  --cluster helix-ai-cluster \
  --service <service-name> \
  --desired-count 2 \
  --region us-east-1
```

**Auto-scaling**: Configure via CDK or AWS Console (ECS Service → Auto Scaling tab)

## Troubleshooting

### Common Issues

#### 1. ECS Service Won't Start

**Symptoms**: Tasks keep stopping or won't start

**Solutions**:
- Check CloudWatch Logs for errors
- Verify ECR image exists and is accessible
- Check security group rules (ALB → ECS)
- Verify health check endpoint is working
- Check task definition resource limits (CPU/memory)

#### 2. Backend Not Accessible

**Symptoms**: Can't reach backend API

**Solutions**:
- Verify ALB security group allows inbound traffic
- Check ECS tasks are running and healthy
- Verify target group health checks are passing
- Check ALB listener rules

#### 3. Frontend Build Fails

**Symptoms**: Build errors during deployment

**Solutions**:
- Verify Node.js version (18+)
- Check for dependency issues: `npm ci` in frontend/
- Verify `VITE_API_BASE_URL` is set correctly
- Check frontend build logs

#### 4. CloudFront Not Updating

**Symptoms**: Changes not visible after deployment

**Solutions**:
- Verify CloudFront invalidation was created
- Check invalidation status (can take 5-15 minutes)
- Verify S3 bucket contents were updated
- Check CloudFront cache policies

#### 5. CDK Deployment Fails

**Symptoms**: Stack deployment errors

**Solutions**:
- Verify AWS credentials are configured
- Check IAM permissions are sufficient
- Review CloudFormation events in AWS Console
- Ensure CDK is bootstrapped: `cdk bootstrap`
- Check for resource naming conflicts

### Getting Help

1. **Check CloudWatch Logs**: Most errors are logged here
2. **Review CloudFormation Events**: For infrastructure issues
3. **Check ECS Service Events**: For container/service issues
4. **View ALB Access Logs**: For traffic/routing issues

## CI/CD Setup

### GitHub Actions

The project includes a GitHub Actions workflow (`.github/workflows/deploy.yml`) that automatically deploys on push to `main` branch.

### Required GitHub Secrets

Add these secrets in your GitHub repository settings:

1. `AWS_REGION`: AWS region (e.g., `us-east-1`)
2. `AWS_ACCOUNT_ID`: Your AWS account ID
3. `AWS_OIDC_ROLE_ARN`: IAM role ARN for GitHub OIDC authentication
4. `ECR_REPOSITORY`: ECR repository name (e.g., `helix-backend`)
5. `S3_BUCKET`: S3 bucket name for frontend
6. `VITE_API_BASE_URL`: Backend API URL for frontend build
7. `CLOUDFRONT_DISTRIBUTION_ID`: CloudFront distribution ID (optional)
8. `ECS_CLUSTER_NAME`: ECS cluster name (optional - for auto ECS updates)
9. `ECS_SERVICE_NAME`: ECS service name (optional)
10. `ECS_TASK_DEFINITION_FAMILY`: Task definition family (optional)

### Setting Up GitHub OIDC

1. **Create IAM OIDC Identity Provider** (if not exists):
   ```bash
   aws iam create-open-id-connect-provider \
     --url https://token.actions.githubusercontent.com \
     --client-id-list sts.amazonaws.com \
     --thumbprint-list 6938fd4d98bab03faadb97b34396831e3780aea1
   ```

2. **Create IAM Role for GitHub Actions**:
   - Trust policy: Allow `token.actions.githubusercontent.com`
   - Conditions: Match repository and branch
   - Permissions: ECR, ECS, S3, CloudFront access

3. **Get Role ARN** and add to GitHub secrets as `AWS_OIDC_ROLE_ARN`

### Workflow Trigger

The workflow triggers on:
- Push to `main` branch
- Manual trigger via GitHub Actions UI

## Next Steps

### Production Hardening

1. **HTTPS Setup**:
   - Request ACM certificate
   - Add HTTPS listener to ALB
   - Configure CloudFront with certificate

2. **Custom Domain**:
   - Set up Route53 hosted zone
   - Create DNS records pointing to ALB/CloudFront
   - Update CDK to use custom domains

3. **Secrets Management**:
   - Store API keys in AWS Secrets Manager
   - Update ECS task definition to load secrets
   - Remove hardcoded credentials

4. **Monitoring & Alarms**:
   - Set up CloudWatch alarms
   - Configure SNS notifications
   - Set up X-Ray tracing (optional)

5. **Backup & Disaster Recovery**:
   - Enable ECR image scanning
   - Set up backup for session data (if using S3/EFS)
   - Document recovery procedures

6. **Security**:
   - Enable WAF on CloudFront
   - Restrict S3 bucket policies
   - Review security group rules
   - Enable AWS GuardDuty

### Cost Optimization

1. **Right-size resources**: Adjust ECS task CPU/memory based on actual usage
2. **Auto-scaling**: Scale down during off-hours
3. **CloudFront caching**: Optimize cache policies
4. **NAT Gateway**: Use NAT instances for lower cost (less HA)
5. **Reserved Capacity**: Consider EC2 Reserved Instances if moving from Fargate

## Additional Resources

- [AWS ECS Documentation](https://docs.aws.amazon.com/ecs/)
- [AWS CDK Documentation](https://docs.aws.amazon.com/cdk/)
- [AWS CloudFront Documentation](https://docs.aws.amazon.com/cloudfront/)
- [Helix.AI Infrastructure README](../infrastructure/README.md)

---

**Questions or Issues?** Open an issue on GitHub or contact the development team.

