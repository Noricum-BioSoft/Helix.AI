"""
Helix.AI CDK Stack
Creates all AWS infrastructure components for the Helix.AI platform
"""

from aws_cdk import (
    Stack,
    aws_ecr as ecr,
    aws_ecs as ecs,
    aws_ec2 as ec2,
    aws_elasticloadbalancingv2 as elbv2,
    aws_iam as iam,
    aws_s3 as s3,
    aws_s3_deployment as s3_deployment,
    aws_cloudfront as cloudfront,
    aws_cloudfront_origins as origins,
    aws_logs as logs,
    aws_certificatemanager as acm,
    aws_secretsmanager as secretsmanager,
    CfnOutput,
    Duration,
    RemovalPolicy,
)
from constructs import Construct
import os


class HelixStack(Stack):
    """Main stack for Helix.AI infrastructure"""

    def __init__(self, scope: Construct, construct_id: str, **kwargs) -> None:
        super().__init__(scope, construct_id, **kwargs)

        # Configuration
        project_name = "helix-ai"
        backend_port = 8001
        
        # Custom domain configuration (optional)
        # Can be set via environment variables or CDK context
        # Environment variables take precedence
        helix_domain = os.environ.get("HELIX_DOMAIN") or self.node.try_get_context("helixDomain") or None
        helix_api_domain = os.environ.get("HELIX_API_DOMAIN") or self.node.try_get_context("helixApiDomain") or helix_domain
        acm_cert_arn_cloudfront = os.environ.get("ACM_CERTIFICATE_ARN_CLOUDFRONT") or self.node.try_get_context("acmCertificateArnCloudFront") or None
        acm_cert_arn_alb = os.environ.get("ACM_CERTIFICATE_ARN_ALB") or self.node.try_get_context("acmCertificateArnAlb") or None
        
        # Import certificates if provided (explicitly check for non-empty strings)
        alb_certificate = None
        cloudfront_certificate = None
        
        if acm_cert_arn_alb and acm_cert_arn_alb.strip():
            # Certificate must be in the same region as the ALB
            try:
                alb_certificate = acm.Certificate.from_certificate_arn(
                    self, "ALBCertificate", acm_cert_arn_alb.strip()
                )
            except Exception as e:
                # Log but don't fail - certificate might not exist yet
                print(f"Warning: Could not import ALB certificate: {e}")
                alb_certificate = None
        
        if acm_cert_arn_cloudfront and acm_cert_arn_cloudfront.strip():
            # CloudFront certificates must be in us-east-1
            # We need to import from us-east-1 region
            try:
                cloudfront_certificate = acm.Certificate.from_certificate_arn(
                    self, "CloudFrontCertificate", acm_cert_arn_cloudfront.strip()
                )
            except Exception as e:
                # Log but don't fail - certificate might not exist yet
                print(f"Warning: Could not import CloudFront certificate: {e}")
                cloudfront_certificate = None
        
        # ==========================================
        # 1. ECR Repository for Backend
        # ==========================================
        # Reference existing repository instead of creating new one
        # (to avoid ResourceExistenceCheck validation error)
        ecr_repo = ecr.Repository.from_repository_name(
            self,
            "BackendRepository",
            repository_name=f"{project_name}-backend",
        )
        
        CfnOutput(
            self,
            "ECRRepositoryURI",
            value=ecr_repo.repository_uri,
            description="ECR Repository URI for backend Docker images",
        )

        # ==========================================
        # 2. VPC and Networking
        # ==========================================
        # Use available AZs for us-west-1 (us-west-1b and us-west-1c are available in this account)
        vpc = ec2.Vpc(
            self,
            "VPC",
            availability_zones=["us-west-1b", "us-west-1c"],  # Match available AZs in the stack
            nat_gateways=1,  # Single NAT gateway for cost savings (increase for production)
            subnet_configuration=[
                ec2.SubnetConfiguration(
                    subnet_type=ec2.SubnetType.PUBLIC,
                    name="Public",
                    cidr_mask=24,
                ),
                ec2.SubnetConfiguration(
                    subnet_type=ec2.SubnetType.PRIVATE_WITH_EGRESS,
                    name="Private",
                    cidr_mask=24,
                ),
            ],
        )
        
        CfnOutput(
            self,
            "VPCId",
            value=vpc.vpc_id,
            description="VPC ID",
        )

        # ==========================================
        # 3. Security Groups
        # ==========================================
        # ALB Security Group
        alb_sg = ec2.SecurityGroup(
            self,
            "ALBSecurityGroup",
            vpc=vpc,
            description="Security group for Application Load Balancer",
            allow_all_outbound=True,
        )
        alb_sg.add_ingress_rule(
            peer=ec2.Peer.any_ipv4(),
            connection=ec2.Port.tcp(80),
            description="Allow HTTP traffic",
        )
        alb_sg.add_ingress_rule(
            peer=ec2.Peer.any_ipv4(),
            connection=ec2.Port.tcp(443),
            description="Allow HTTPS traffic",
        )

        # ECS Service Security Group
        ecs_sg = ec2.SecurityGroup(
            self,
            "ECSSecurityGroup",
            vpc=vpc,
            description="Security group for ECS tasks",
            allow_all_outbound=True,
        )
        ecs_sg.add_ingress_rule(
            peer=alb_sg,
            connection=ec2.Port.tcp(backend_port),
            description="Allow traffic from ALB",
        )

        # ==========================================
        # 4. ECS Cluster and Fargate Service
        # ==========================================
        cluster = ecs.Cluster(
            self,
            "Cluster",
            vpc=vpc,
            cluster_name=f"{project_name}-cluster",
            # Note: container_insights is deprecated in CDK
            # Container Insights can be enabled via console or using enable_container_insights_v2 property
            # For now, we'll omit it to avoid the deprecation warning
        )
        
        CfnOutput(
            self,
            "ECSClusterName",
            value=cluster.cluster_name,
            description="ECS Cluster Name",
        )

        # Task Execution Role
        task_execution_role = iam.Role(
            self,
            "TaskExecutionRole",
            assumed_by=iam.ServicePrincipal("ecs-tasks.amazonaws.com"),
            managed_policies=[
                iam.ManagedPolicy.from_aws_managed_policy_name(
                    "service-role/AmazonECSTaskExecutionRolePolicy"
                ),
            ],
        )
        ecr_repo.grant_pull(task_execution_role)

        # Task Role
        task_role = iam.Role(
            self,
            "TaskRole",
            assumed_by=iam.ServicePrincipal("ecs-tasks.amazonaws.com"),
            description="Role for ECS tasks to access AWS services",
        )
        
        # Grant S3 permissions for file operations
        # The backend needs to:
        # - Read: Download uploaded files, access job results from EMR, access dataset files
        # - Write: Copy job results to session paths, upload HTML visualizations
        # Note: For production, consider restricting to specific buckets
        task_role.add_to_policy(
            iam.PolicyStatement(
                effect=iam.Effect.ALLOW,
                actions=[
                    "s3:GetObject",
                    "s3:HeadObject",
                    "s3:PutObject",
                    "s3:CopyObject",
                    "s3:ListBucket",
                ],
                resources=["*"],  # Allow access to all buckets - restrict in production if needed
            )
        )
        
        # Add permissions for accessing other AWS services as needed
        # e.g., Secrets Manager, Parameter Store, EMR, etc.
        
        # Grant task execution role permission to read secrets from Secrets Manager
        # This allows the container to access API keys stored in Secrets Manager
        # Use the existing production secrets
        task_execution_role.add_to_policy(
            iam.PolicyStatement(
                effect=iam.Effect.ALLOW,
                actions=["secretsmanager:GetSecretValue"],
                resources=[
                    f"arn:aws:secretsmanager:{self.region}:{self.account}:secret:helix-ai-production-*"
                ]
            )
        )

        # Log Group
        log_group = logs.LogGroup(
            self,
            "LogGroup",
            log_group_name=f"/ecs/{project_name}",
            retention=logs.RetentionDays.ONE_WEEK,  # Adjust retention as needed
            removal_policy=RemovalPolicy.DESTROY,
        )

        # Task Definition
        task_definition = ecs.FargateTaskDefinition(
            self,
            "TaskDefinition",
            execution_role=task_execution_role,
            task_role=task_role,
            memory_limit_mib=2048,  # 2GB - adjust based on requirements
            cpu=1024,  # 1 vCPU - adjust based on requirements
        )

        # Reference API key secrets from Secrets Manager
        # Use the existing secrets that were created for Copilot deployment
        # These secrets already exist: helix-ai-production-OPENAI_API_KEY and helix-ai-production-DEEPSEEK_API_KEY
        openai_secret_name = "helix-ai-production-OPENAI_API_KEY"
        deepseek_secret_name = "helix-ai-production-DEEPSEEK_API_KEY"
        
        # Reference the existing secrets
        openai_secret = secretsmanager.Secret.from_secret_name_v2(
            self,
            "OpenAISecret",
            secret_name=openai_secret_name
        )
        
        deepseek_secret = secretsmanager.Secret.from_secret_name_v2(
            self,
            "DeepSeekSecret",
            secret_name=deepseek_secret_name
        )
        
        # Build secrets dict for container
        container_secrets = {
            "OPENAI_API_KEY": ecs.Secret.from_secrets_manager(openai_secret),
            "DEEPSEEK_API_KEY": ecs.Secret.from_secrets_manager(deepseek_secret),
        }

        # Container Definition
        container = task_definition.add_container(
            "BackendContainer",
            image=ecs.ContainerImage.from_ecr_repository(
                repository=ecr_repo,
                tag="latest",
            ),
            logging=ecs.LogDrivers.aws_logs(
                stream_prefix="helix-backend",
                log_group=log_group,
            ),
            environment={
                "PYTHONUNBUFFERED": "1",
            },
            secrets=container_secrets if container_secrets else None,  # Only add if secrets exist
            health_check=ecs.HealthCheck(
                command=["CMD-SHELL", f"curl -f http://localhost:{backend_port}/health || exit 1"],
                interval=Duration.seconds(30),
                timeout=Duration.seconds(5),
                retries=3,
                start_period=Duration.seconds(60),
            ),
        )
        container.add_port_mappings(
            ecs.PortMapping(
                container_port=backend_port,
                protocol=ecs.Protocol.TCP,
            )
        )

        # ECS Service
        service = ecs.FargateService(
            self,
            "Service",
            cluster=cluster,
            task_definition=task_definition,
            desired_count=0,  # Set to 0 temporarily to allow stack deployment, then debug ECS service separately
            security_groups=[ecs_sg],
            vpc_subnets=ec2.SubnetSelection(subnet_type=ec2.SubnetType.PRIVATE_WITH_EGRESS),
            health_check_grace_period=Duration.seconds(60),
            # Deployment configuration
            max_healthy_percent=200,
            min_healthy_percent=100,
        )
        
        CfnOutput(
            self,
            "ECSServiceName",
            value=service.service_name,
            description="ECS Service Name",
        )
        
        CfnOutput(
            self,
            "ECSTaskDefinitionFamily",
            value=task_definition.family,
            description="ECS Task Definition Family",
        )

        # ==========================================
        # 5. Application Load Balancer
        # ==========================================
        alb = elbv2.ApplicationLoadBalancer(
            self,
            "ALB",
            vpc=vpc,
            internet_facing=True,
            security_group=alb_sg,
        )
        
        # Health Check Target Group
        target_group = elbv2.ApplicationTargetGroup(
            self,
            "TargetGroup",
            vpc=vpc,
            port=backend_port,
            protocol=elbv2.ApplicationProtocol.HTTP,
            target_type=elbv2.TargetType.IP,
            health_check=elbv2.HealthCheck(
                path="/health",
                interval=Duration.seconds(30),
                timeout=Duration.seconds(5),
                healthy_threshold_count=2,
                unhealthy_threshold_count=3,
            ),
            deregistration_delay=Duration.seconds(30),
        )
        
        # Register ECS service with target group
        service.attach_to_application_target_group(target_group)
        
        # HTTPS Listener (if certificate is provided)
        https_listener = None
        if alb_certificate:
            https_listener = alb.add_listener(
                "HTTPSListener",
                port=443,
                protocol=elbv2.ApplicationProtocol.HTTPS,
                certificates=[alb_certificate],
                default_target_groups=[target_group],
            )
        
        # HTTP Listener
        # If HTTPS is configured, redirect HTTP to HTTPS
        # Otherwise, forward to target group
        if https_listener:
            # Redirect HTTP to HTTPS
            http_listener = alb.add_listener(
                "HTTPListener",
                port=80,
                default_action=elbv2.ListenerAction.redirect(
                    protocol="HTTPS",
                    port="443",
                    permanent=True,
                ),
            )
        else:
            # Forward HTTP to target group (no certificate configured)
            http_listener = alb.add_listener(
                "HTTPListener",
                port=80,
                default_target_groups=[target_group],
            )
        
        CfnOutput(
            self,
            "ALBDNSName",
            value=alb.load_balancer_dns_name,
            description="Application Load Balancer DNS Name",
        )
        
        # Output custom domain if configured
        if helix_api_domain:
            CfnOutput(
                self,
                "BackendCustomDomain",
                value=helix_api_domain,
                description="Backend API Custom Domain",
            )

        # ==========================================
        # 6. S3 Bucket for Frontend
        # ==========================================
        # Reference existing bucket instead of creating new one
        # (to avoid ResourceExistenceCheck validation error)
        frontend_bucket = s3.Bucket.from_bucket_name(
            self,
            "FrontendBucket",
            bucket_name=f"{project_name}-frontend-{self.account}-{self.region}",
        )
        
        CfnOutput(
            self,
            "FrontendBucketName",
            value=frontend_bucket.bucket_name,
            description="S3 Bucket Name for Frontend",
        )
        
        CfnOutput(
            self,
            "FrontendBucketWebsiteURL",
            value=frontend_bucket.bucket_website_url,
            description="Frontend Website URL (S3)",
        )

        # ==========================================
        # 7. CloudFront Distribution
        # ==========================================
        # CloudFront Origin Access Identity (OAI) - for private S3 access
        oai = cloudfront.OriginAccessIdentity(
            self,
            "CloudFrontOAI",
            comment=f"OAI for {project_name} frontend",
        )
        
        # Grant CloudFront access to S3 bucket
        frontend_bucket.grant_read(oai)

        # CloudFront Distribution
        # Note: Default behavior with '*' pattern handles all routes for React Router
        distribution_props = {
            "default_behavior": cloudfront.BehaviorOptions(
                origin=origins.S3Origin(
                    bucket=frontend_bucket,
                    origin_access_identity=oai,
                ),
                viewer_protocol_policy=cloudfront.ViewerProtocolPolicy.REDIRECT_TO_HTTPS,
                allowed_methods=cloudfront.AllowedMethods.ALLOW_GET_HEAD_OPTIONS,
                cached_methods=cloudfront.CachedMethods.CACHE_GET_HEAD,
                compress=True,
                cache_policy=cloudfront.CachePolicy.CACHING_OPTIMIZED,
            ),
            "error_responses": [
                cloudfront.ErrorResponse(
                    http_status=404,
                    response_http_status=200,
                    response_page_path="/index.html",
                    ttl=Duration.minutes(5),
                ),
                cloudfront.ErrorResponse(
                    http_status=403,
                    response_http_status=200,
                    response_page_path="/index.html",
                    ttl=Duration.minutes(5),
                ),
            ],
            "price_class": cloudfront.PriceClass.PRICE_CLASS_100,  # Use only North America and Europe
            "comment": f"CloudFront distribution for {project_name} frontend",
            "enable_logging": False,  # Enable if you want access logs
        }
        
        # Add custom domain and certificate if provided
        if helix_domain and cloudfront_certificate:
            distribution_props["domain_names"] = [helix_domain]
            distribution_props["certificate"] = cloudfront_certificate
        
        distribution = cloudfront.Distribution(
            self,
            "CloudFrontDistribution",
            **distribution_props
        )
        
        CfnOutput(
            self,
            "CloudFrontDistributionId",
            value=distribution.distribution_id,
            description="CloudFront Distribution ID",
        )
        
        CfnOutput(
            self,
            "CloudFrontDomainName",
            value=distribution.distribution_domain_name,
            description="CloudFront Distribution Domain Name",
        )
        
        # Output custom domain if configured
        if helix_domain:
            CfnOutput(
                self,
                "FrontendCustomDomain",
                value=helix_domain,
                description="Frontend Custom Domain",
            )

        # ==========================================
        # 8. Outputs
        # ==========================================
        # Backend URL - use HTTPS if certificate is configured, otherwise HTTP
        backend_protocol = "https" if alb_certificate else "http"
        CfnOutput(
            self,
            "BackendAPIURL",
            value=f"{backend_protocol}://{alb.load_balancer_dns_name}",
            description="Backend API URL",
        )
        
        # If custom domain is configured, also output that
        if helix_api_domain and alb_certificate:
            CfnOutput(
                self,
                "BackendAPIURLCustom",
                value=f"https://{helix_api_domain}",
                description="Backend API URL (Custom Domain)",
            )
        
        # Frontend URL - always HTTPS (CloudFront)
        CfnOutput(
            self,
            "FrontendURL",
            value=f"https://{distribution.distribution_domain_name}",
            description="Frontend URL (CloudFront)",
        )
        
        # If custom domain is configured, also output that
        if helix_domain:
            CfnOutput(
                self,
                "FrontendURLCustom",
                value=f"https://{helix_domain}",
                description="Frontend URL (Custom Domain)",
            )

