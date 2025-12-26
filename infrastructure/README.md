# Helix.AI AWS Infrastructure

This directory contains AWS CDK code to deploy the complete infrastructure for Helix.AI on AWS.

## Prerequisites

1. **AWS CLI** configured with appropriate credentials
2. **AWS CDK CLI** installed:
   ```bash
   npm install -g aws-cdk
   ```
3. **Python 3.9+** with pip
4. **Docker** (for building and testing locally)

## Setup

1. **Install Python dependencies**:
   ```bash
   cd infrastructure
   pip install -r requirements.txt
   ```

2. **Bootstrap CDK** (first time only):
   ```bash
   cdk bootstrap
   ```

3. **Set environment variables** (optional):
   ```bash
   export CDK_DEFAULT_ACCOUNT=$(aws sts get-caller-identity --query Account --output text)
   export CDK_DEFAULT_REGION=us-east-1
   export STACK_NAME=HelixAIStack
   ```

## Deploy Infrastructure

1. **Synthesize CloudFormation template** (preview changes):
   ```bash
   cdk synth
   ```

2. **Deploy the stack**:
   ```bash
   cdk deploy
   ```

3. **Deploy to a specific account/region**:
   ```bash
   cdk deploy --profile your-profile --region us-west-2
   ```

## What Gets Deployed

The stack creates:

1. **ECR Repository**: Docker image repository for backend
2. **VPC**: Virtual Private Cloud with public and private subnets
3. **ECS Cluster**: Fargate cluster for running containers
4. **Application Load Balancer**: HTTP/HTTPS load balancer for backend
5. **ECS Service**: Fargate service running the backend container
6. **S3 Bucket**: Static website hosting for frontend
7. **CloudFront Distribution**: CDN for frontend with HTTPS
8. **Security Groups**: Network security configurations
9. **IAM Roles**: Permissions for ECS tasks
10. **CloudWatch Logs**: Centralized logging

## Configuration

### Adjusting Resources

Edit `helix_infrastructure/helix_stack.py` to customize:

- **ECS Task CPU/Memory**: Modify `memory_limit_mib` and `cpu` in `FargateTaskDefinition`
- **Desired Task Count**: Change `desired_count` in `FargateService`
- **VPC Configuration**: Modify `max_azs` and `nat_gateways` in `Vpc`
- **CloudFront Cache Policy**: Change `cache_policy` in `Distribution`

### Environment Variables

You can pass environment variables via CDK context:

```bash
cdk deploy -c region=us-west-2 -c project_name=helix-ai
```

## Post-Deployment

After deploying the infrastructure:

1. **Note the outputs** from the stack deployment (ECR URI, ALB DNS, etc.)
2. **Update `scripts/aws/deploy.config`** with:
   - ECR repository name
   - S3 bucket name
   - CloudFront distribution ID
   - ECS cluster/service names
   - ALB URL for `VITE_API_BASE_URL`
3. **Build and push backend image**:
   ```bash
   ./scripts/aws/deploy.sh
   ```

## Useful Commands

```bash
# List all stacks
cdk list

# View stack diff
cdk diff

# Deploy specific stack
cdk deploy HelixAIStack

# Destroy stack (removes all resources)
cdk destroy

# View stack outputs
aws cloudformation describe-stacks --stack-name HelixAIStack --query 'Stacks[0].Outputs'
```

## Cleanup

To remove all infrastructure:

```bash
cdk destroy
```

**Note**: S3 buckets and ECR repositories are set to `RETAIN` by default, so they won't be deleted. You'll need to manually empty and delete them if desired.

## Cost Optimization

- **NAT Gateways**: Currently set to 1. For production, consider using 2+ for high availability, or use NAT Instances for lower cost.
- **ECS Tasks**: Start with 1 task. Scale based on usage.
- **CloudWatch Logs**: Retention is set to 7 days. Adjust based on requirements.
- **CloudFront**: Using `PRICE_CLASS_100` (US/Europe only). Use `PRICE_CLASS_ALL` for global distribution.

## Troubleshooting

### CDK Bootstrap Issues

If you get bootstrap errors:
```bash
cdk bootstrap aws://ACCOUNT-ID/REGION
```

### Permission Issues

Ensure your AWS credentials have permissions for:
- EC2 (VPC, Security Groups, Load Balancers)
- ECS (Cluster, Services, Task Definitions)
- ECR (Repository operations)
- S3 (Bucket creation)
- CloudFront (Distribution creation)
- CloudFormation (Stack operations)
- IAM (Role creation)

### ECS Service Won't Start

1. Check ECS service events in AWS Console
2. Check CloudWatch Logs for container logs
3. Verify ECR image exists and is accessible
4. Check security group rules
5. Verify health check endpoint is responding

## Next Steps

1. Set up HTTPS certificate with ACM for ALB
2. Configure custom domain names
3. Set up auto-scaling for ECS service
4. Configure CloudWatch alarms
5. Set up secrets in AWS Secrets Manager
6. Configure WAF rules for CloudFront
7. Set up backup and disaster recovery




