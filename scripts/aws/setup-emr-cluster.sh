#!/bin/bash
# Setup EMR cluster for PySpark FASTQ processing
# This script creates necessary IAM roles and launches an EMR cluster

set -e

# Configuration
CLUSTER_NAME="${EMR_CLUSTER_NAME:-helix-ai-fastqc}"
REGION="${AWS_REGION:-us-east-1}"
S3_BUCKET="${S3_DATASET_BUCKET:-noricum-ngs-data}"
RELEASE_LABEL="emr-6.15.0"  # Latest stable EMR with Spark 3.5
INSTANCE_TYPE="${EMR_INSTANCE_TYPE:-m5.xlarge}"
INSTANCE_COUNT="${EMR_INSTANCE_COUNT:-3}"
EC2_KEY_NAME="${EC2_KEY_NAME:-}"  # Optional: your EC2 key pair name

echo "=========================================="
echo "EMR Cluster Setup for Helix.AI FASTQ Processing"
echo "=========================================="
echo "Cluster Name: $CLUSTER_NAME"
echo "Region: $REGION"
echo "S3 Bucket: $S3_BUCKET"
echo "Instance Type: $INSTANCE_TYPE"
echo "Instance Count: $INSTANCE_COUNT"
echo ""

# Check AWS CLI
if ! command -v aws &> /dev/null; then
    echo "Error: AWS CLI not found. Please install it first."
    exit 1
fi

# Check AWS credentials
if ! aws sts get-caller-identity &> /dev/null; then
    echo "Error: AWS credentials not configured. Run 'aws configure' first."
    exit 1
fi

# Step 1: Create EMR Service Role (if it doesn't exist)
echo "Step 1: Setting up EMR Service Role..."
SERVICE_ROLE_NAME="EMR_DefaultRole"
SERVICE_ROLE_ARN=""

if aws iam get-role --role-name "$SERVICE_ROLE_NAME" &>/dev/null; then
    echo "✅ EMR Service Role already exists"
    SERVICE_ROLE_ARN=$(aws iam get-role --role-name "$SERVICE_ROLE_NAME" --query 'Role.Arn' --output text)
else
    echo "Creating EMR Service Role..."
    
    # Create trust policy
    cat > /tmp/emr-service-trust-policy.json << 'EOF'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "elasticmapreduce.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

    # Create role
    aws iam create-role \
        --role-name "$SERVICE_ROLE_NAME" \
        --assume-role-policy-document file:///tmp/emr-service-trust-policy.json

    # Attach managed policy
    aws iam attach-role-policy \
        --role-name "$SERVICE_ROLE_NAME" \
        --policy-arn arn:aws:iam::aws:policy/service-role/AmazonElasticMapReduceRole

    SERVICE_ROLE_ARN=$(aws iam get-role --role-name "$SERVICE_ROLE_NAME" --query 'Role.Arn' --output text)
    echo "✅ Created EMR Service Role: $SERVICE_ROLE_ARN"
fi

# Step 2: Create EMR EC2 Instance Profile (if it doesn't exist)
echo ""
echo "Step 2: Setting up EMR EC2 Instance Profile..."
INSTANCE_PROFILE_NAME="EMR_EC2_DefaultRole"
INSTANCE_PROFILE_ARN=""

if aws iam get-instance-profile --instance-profile-name "$INSTANCE_PROFILE_NAME" &>/dev/null; then
    echo "✅ EMR EC2 Instance Profile already exists"
    INSTANCE_PROFILE_ARN=$(aws iam get-instance-profile --instance-profile-name "$INSTANCE_PROFILE_NAME" --query 'InstanceProfile.Arn' --output text)
else
    echo "Creating EMR EC2 Instance Profile..."
    
    # Create EC2 role
    EC2_ROLE_NAME="EMR_EC2_DefaultRole"
    
    if ! aws iam get-role --role-name "$EC2_ROLE_NAME" &>/dev/null; then
        # Create trust policy for EC2
        cat > /tmp/emr-ec2-trust-policy.json << 'EOF'
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Principal": {
        "Service": "ec2.amazonaws.com"
      },
      "Action": "sts:AssumeRole"
    }
  ]
}
EOF

        aws iam create-role \
            --role-name "$EC2_ROLE_NAME" \
            --assume-role-policy-document file:///tmp/emr-ec2-trust-policy.json

        # Attach managed policies
        aws iam attach-role-policy \
            --role-name "$EC2_ROLE_NAME" \
            --policy-arn arn:aws:iam::aws:policy/service-role/AmazonElasticMapReduceforEC2Role

        # Add S3 access policy
        cat > /tmp/emr-s3-policy.json << EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Effect": "Allow",
      "Action": [
        "s3:GetObject",
        "s3:PutObject",
        "s3:ListBucket"
      ],
      "Resource": [
        "arn:aws:s3:::${S3_BUCKET}/*",
        "arn:aws:s3:::${S3_BUCKET}",
        "arn:aws:s3:::helix-ai-frontend-*/*",
        "arn:aws:s3:::helix-ai-frontend-*"
      ]
    }
  ]
}
EOF

        aws iam put-role-policy \
            --role-name "$EC2_ROLE_NAME" \
            --policy-name S3Access \
            --policy-document file:///tmp/emr-s3-policy.json
    fi

    # Create instance profile
    aws iam create-instance-profile --instance-profile-name "$INSTANCE_PROFILE_NAME"
    aws iam add-role-to-instance-profile \
        --instance-profile-name "$INSTANCE_PROFILE_NAME" \
        --role-name "$EC2_ROLE_NAME"

    # Wait for instance profile to be ready
    echo "Waiting for instance profile to be ready..."
    sleep 10

    INSTANCE_PROFILE_ARN=$(aws iam get-instance-profile --instance-profile-name "$INSTANCE_PROFILE_NAME" --query 'InstanceProfile.Arn' --output text)
    echo "✅ Created EMR EC2 Instance Profile: $INSTANCE_PROFILE_ARN"
fi

# Step 3: Create bootstrap script for Python dependencies
echo ""
echo "Step 3: Creating bootstrap script..."
BOOTSTRAP_SCRIPT="s3://${S3_BUCKET}/emr-bootstrap/fastqc-bootstrap.sh"

# Create bootstrap script locally
mkdir -p /tmp/emr-bootstrap
cat > /tmp/emr-bootstrap/fastqc-bootstrap.sh << 'BOOTSTRAP_EOF'
#!/bin/bash
# Bootstrap script for EMR cluster
# Installs Python dependencies for FASTQ processing

# Update system
sudo yum update -y

# Install additional Python packages
sudo pip3 install biopython pandas numpy matplotlib seaborn plotly

# Install FASTQ parsing libraries
sudo pip3 install pyspark[fastq] 2>/dev/null || echo "pyspark already installed"

# Create directories for logs
sudo mkdir -p /mnt/var/log/fastqc
sudo chmod 777 /mnt/var/log/fastqc

echo "Bootstrap script completed successfully"
BOOTSTRAP_EOF

chmod +x /tmp/emr-bootstrap/fastqc-bootstrap.sh

# Upload bootstrap script to S3
echo "Uploading bootstrap script to S3..."
aws s3 cp /tmp/emr-bootstrap/fastqc-bootstrap.sh "$BOOTSTRAP_SCRIPT" || {
    # Create bucket if it doesn't exist
    aws s3 mb "s3://${S3_BUCKET}" --region "$REGION" 2>/dev/null || true
    aws s3 cp /tmp/emr-bootstrap/fastqc-bootstrap.sh "$BOOTSTRAP_SCRIPT"
}
echo "✅ Bootstrap script uploaded"

# Step 4: Launch EMR cluster
echo ""
echo "Step 4: Launching EMR cluster..."
echo "This may take 5-10 minutes..."


# Create configurations file
cat > /tmp/emr-configurations.json << 'CONFIG_EOF'
[
  {
    "Classification": "spark-defaults",
    "Properties": {
      "spark.sql.adaptive.enabled": "true",
      "spark.sql.adaptive.coalescePartitions.enabled": "true",
      "spark.serializer": "org.apache.spark.serializer.KryoSerializer",
      "spark.sql.execution.arrow.pyspark.enabled": "true"
    }
  },
  {
    "Classification": "spark-env",
    "Configurations": [
      {
        "Classification": "export",
        "Properties": {
          "PYSPARK_PYTHON": "/usr/bin/python3",
          "PYSPARK_DRIVER_PYTHON": "/usr/bin/python3"
        }
      }
    ]
  }
]
CONFIG_EOF

# Create instance groups JSON file
# Note: BidPrice indicates spot instance; no BidPrice = on-demand
cat > /tmp/emr-instance-groups.json << EOF
[
  {
    "Name": "Master",
    "InstanceGroupType": "MASTER",
    "InstanceType": "$INSTANCE_TYPE",
    "InstanceCount": 1
  },
  {
    "Name": "Core",
    "InstanceGroupType": "CORE",
    "InstanceType": "$INSTANCE_TYPE",
    "InstanceCount": $((INSTANCE_COUNT - 1)),
    "BidPrice": "0.10"
  }
]
EOF

# Launch cluster
# Note: --ec2-attributes InstanceProfile specifies the EC2 instance profile
CLUSTER_ID=$(aws emr create-cluster \
    --name "$CLUSTER_NAME" \
    --release-label "$RELEASE_LABEL" \
    --applications Name=Spark Name=Hadoop Name=Zeppelin \
    --instance-groups file:///tmp/emr-instance-groups.json \
    --bootstrap-actions Path="$BOOTSTRAP_SCRIPT",Name="Install Python Dependencies" \
    --configurations file:///tmp/emr-configurations.json \
    --service-role "$SERVICE_ROLE_ARN" \
    --ec2-attributes InstanceProfile="$INSTANCE_PROFILE_ARN" \
    --log-uri "s3://${S3_BUCKET}/emr-logs/" \
    --region "$REGION" \
    --tags Project=Helix.AI Purpose="FASTQ Quality Analysis" \
    --query 'ClusterId' \
    --output text)

echo "✅ EMR Cluster launched!"
echo ""
echo "Cluster ID: $CLUSTER_ID"
echo ""
echo "To check cluster status:"
echo "  aws emr describe-cluster --cluster-id $CLUSTER_ID --region $REGION --query 'Cluster.Status.State' --output text"
echo ""
echo "To get master node public DNS (for SSH):"
echo "  aws emr describe-cluster --cluster-id $CLUSTER_ID --region $REGION --query 'Cluster.MasterPublicDnsName' --output text"
echo ""
echo "To terminate cluster:"
echo "  aws emr terminate-clusters --cluster-ids $CLUSTER_ID --region $REGION"
echo ""
echo "Cluster will be ready in approximately 5-10 minutes."
echo "You can monitor progress in the AWS Console or using the describe-cluster command above."

