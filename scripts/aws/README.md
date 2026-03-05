# AWS Deployment Scripts

This directory contains scripts for automating AWS deployments of Helix.AI.

## Quick Start

1. **Set up deployment configuration**:
   ```bash
   ./setup-deployment.sh
   ```
   This interactive script will guide you through creating your `deploy.config` file.

2. **Deploy infrastructure** (first time only):
   ```bash
   cd ../../infrastructure
   pip install -r requirements.txt
   cdk bootstrap
   cdk deploy
   ```

3. **Deploy application**:
   ```bash
   ./deploy.sh
   ```

## Scripts Overview

### `deploy.sh`
Master deployment script that orchestrates the complete deployment process:
- Builds and pushes backend Docker image to ECR
- Updates ECS service (if configured)
- Builds and deploys frontend to S3
- Invalidates CloudFront cache (if configured)

**Usage**:
```bash
./deploy.sh [config-file]
```

### `setup-deployment.sh`
Interactive setup script that helps you create your deployment configuration file.

**Usage**:
```bash
./setup-deployment.sh
```

### `ecr_push_backend.sh`
Builds and pushes the backend Docker image to AWS ECR.

**Usage**:
```bash
./ecr_push_backend.sh <AWS_ACCOUNT_ID> <AWS_REGION> <ECR_REPO_NAME> [<IMAGE_TAG>]
```

**Example**:
```bash
./ecr_push_backend.sh 123456789012 us-east-1 helix-backend latest
```

### `s3_sync_frontend.sh`
Builds the frontend and syncs it to S3, optionally invalidating CloudFront.

**Usage**:
```bash
./s3_sync_frontend.sh <AWS_REGION> <S3_BUCKET_NAME> [<CLOUDFRONT_DIST_ID>]
```

**Example**:
```bash
./s3_sync_frontend.sh us-east-1 helix-frontend-bucket E123ABC456DEF
```

### `get-alb-url.sh`
Helper script to fetch the ALB DNS name from your CloudFormation stack and test if the backend is accessible.

**Usage**:
```bash
./get-alb-url.sh [STACK_NAME]
```

**Example**:
```bash
./get-alb-url.sh HelixAIStack
```

This script is useful for finding the backend API URL to use in `deploy.config` for `VITE_API_BASE_URL`.

## Configuration File

The `deploy.config` file (created from `deploy.config.example`) contains all AWS deployment settings:

```bash
# Copy example config
cp deploy.config.example deploy.config

# Edit with your values
nano deploy.config
```

### Required Settings

- `AWS_REGION`: AWS region (e.g., `us-east-1`)
- `AWS_ACCOUNT_ID`: Your AWS account ID
- `ECR_REPOSITORY`: ECR repository name
- `S3_BUCKET`: S3 bucket name for frontend

### Optional Settings

- `CLOUDFRONT_DISTRIBUTION_ID`: For CloudFront cache invalidation
- `ECS_CLUSTER_NAME`: For automatic ECS service updates
- `ECS_SERVICE_NAME`: For automatic ECS service updates
- `ECS_TASK_DEFINITION_FAMILY`: For automatic ECS service updates
- `VITE_API_BASE_URL`: Backend URL for frontend build

## Deployment Workflow

### Initial Setup

1. **Deploy Infrastructure**:
   ```bash
   cd ../../infrastructure
   cdk deploy
   ```

2. **Capture Output Values**:
   Note the outputs from the CDK deployment:
   - ECR Repository URI
   - ECS Cluster/Service names
   - ALB DNS name
   - S3 Bucket name
   - CloudFront Distribution ID

3. **Configure Deployment**:
   ```bash
   cd ../scripts/aws
   ./setup-deployment.sh
   # Or manually edit deploy.config
   ```

### Regular Deployments

After initial setup, deployments are simple:

```bash
./deploy.sh
```

This will:
1. Build and push backend Docker image
2. Update ECS service (if configured)
3. Build and deploy frontend
4. Invalidate CloudFront (if configured)

## Advanced Usage

### Deploy Only Backend

Set in `deploy.config`:
```bash
SKIP_FRONTEND=true
```

Or use the individual script:
```bash
./ecr_push_backend.sh <args>
```

### Deploy Only Frontend

Set in `deploy.config`:
```bash
SKIP_BACKEND=true
```

Or use the individual script:
```bash
./s3_sync_frontend.sh <args>
```

### Custom Configuration File

Use a different config file:
```bash
./deploy.sh production.config
```

### Skip ECS Update

If you want to manually update ECS later:
```bash
SKIP_ECS_UPDATE=true ./deploy.sh
```

## Troubleshooting

### AWS Credentials Not Configured

Ensure AWS CLI is configured:
```bash
aws configure
aws sts get-caller-identity  # Test credentials
```

### ECR Repository Not Found

The scripts will create the repository if it doesn't exist. If you get permission errors, ensure your AWS user has ECR permissions.

### Docker Build Fails

Check that:
- Docker is running
- Backend Dockerfile is valid
- Sufficient disk space available

### Frontend Build Fails

Check that:
- Node.js 18+ is installed
- Frontend dependencies are installed: `cd frontend && npm ci`
- `VITE_API_BASE_URL` is set correctly

### ECS Service Update Fails

Ensure:
- ECS cluster/service names are correct
- Task definition family name is correct
- IAM permissions allow ECS updates

## See Also

- [AWS Deployment Guide](../../docs/AWS_DEPLOYMENT_GUIDE.md) - Comprehensive deployment documentation
- [Infrastructure README](../../infrastructure/README.md) - AWS CDK infrastructure setup
- [GitHub Actions Workflow](../../.github/workflows/deploy.yml) - CI/CD deployment

