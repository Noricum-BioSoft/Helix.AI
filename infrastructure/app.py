#!/usr/bin/env python3
"""
Helix.AI AWS CDK Application
Deploys complete infrastructure for Helix.AI including:
- ECR repository for backend Docker images
- ECS Fargate cluster with ALB
- S3 bucket for frontend static hosting
- CloudFront distribution for CDN
- VPC, security groups, and IAM roles
"""

import os
import aws_cdk as cdk
from helix_infrastructure.helix_stack import HelixStack

app = cdk.App()

# Get configuration from environment or use defaults
env = cdk.Environment(
    account=os.environ.get("CDK_DEFAULT_ACCOUNT"),
    region=os.environ.get("CDK_DEFAULT_REGION", "us-east-1")
)

stack_name = os.environ.get("STACK_NAME", "HelixAIStack")

HelixStack(
    app,
    stack_name,
    env=env,
    description="Helix.AI Bioinformatics Platform Infrastructure"
)

app.synth()




