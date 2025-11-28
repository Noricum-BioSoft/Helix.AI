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
    CfnOutput,
    Duration,
    RemovalPolicy,
)
from constructs import Construct


class HelixStack(Stack):
    """Main stack for Helix.AI infrastructure"""

    def __init__(self, scope: Construct, construct_id: str, **kwargs) -> None:
        super().__init__(scope, construct_id, **kwargs)

        # Configuration
        project_name = "helix-ai"
        backend_port = 8001
        
        # ==========================================
        # 1. ECR Repository for Backend
        # ==========================================
        ecr_repo = ecr.Repository(
            self,
            "BackendRepository",
            repository_name=f"{project_name}-backend",
            image_scan_on_push=True,
            removal_policy=RemovalPolicy.RETAIN,  # Keep images when stack is deleted
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
        vpc = ec2.Vpc(
            self,
            "VPC",
            max_azs=2,  # Use 2 availability zones for high availability
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
        # Add permissions for accessing other AWS services as needed
        # e.g., S3, Secrets Manager, Parameter Store, etc.

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
            desired_count=1,  # Start with 1, scale as needed
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
        
        # HTTP Listener (redirect to HTTPS in production)
        listener = alb.add_listener(
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

        # ==========================================
        # 6. S3 Bucket for Frontend
        # ==========================================
        frontend_bucket = s3.Bucket(
            self,
            "FrontendBucket",
            bucket_name=f"{project_name}-frontend-{self.account}-{self.region}",
            versioned=False,
            public_read_access=True,
            block_public_access=s3.BlockPublicAccess(
                block_public_acls=False,
                block_public_policy=False,
                ignore_public_acls=False,
                restrict_public_buckets=False,
            ),
            website_index_document="index.html",
            website_error_document="index.html",  # For React Router
            removal_policy=RemovalPolicy.RETAIN,  # Keep bucket when stack is deleted
            auto_delete_objects=False,  # Set to True if you want auto-cleanup
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
        distribution = cloudfront.Distribution(
            self,
            "CloudFrontDistribution",
            default_behavior=cloudfront.BehaviorOptions(
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
            # Custom error responses for React Router (SPA routing)
            error_responses=[
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
            price_class=cloudfront.PriceClass.PRICE_CLASS_100,  # Use only North America and Europe
            comment=f"CloudFront distribution for {project_name} frontend",
            enable_logging=False,  # Enable if you want access logs
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

        # ==========================================
        # 8. Outputs
        # ==========================================
        CfnOutput(
            self,
            "BackendAPIURL",
            value=f"http://{alb.load_balancer_dns_name}",
            description="Backend API URL",
        )
        
        CfnOutput(
            self,
            "FrontendURL",
            value=f"https://{distribution.distribution_domain_name}",
            description="Frontend URL (CloudFront)",
        )

