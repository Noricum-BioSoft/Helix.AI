#!/bin/bash
# Create and configure an EC2 instance for Helix.AI bioinformatics tools
# This script launches an EC2 instance and installs all necessary tools

set -e

# Configuration
REGION="${AWS_REGION:-us-east-1}"
INSTANCE_TYPE="${INSTANCE_TYPE:-t3.medium}"
KEY_NAME="${EC2_KEY_NAME:-}"  # Your EC2 key pair name
AMI_ID="${AMI_ID:-}"  # Optional: custom AMI ID
SECURITY_GROUP_NAME="helix-bioinformatics-tools"

echo "=========================================="
echo "Helix.AI Bioinformatics EC2 Instance Setup"
echo "=========================================="
echo "Region: $REGION"
echo "Instance Type: $INSTANCE_TYPE"
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

# Get default VPC
echo "Step 1: Finding default VPC..."
VPC_ID=$(aws ec2 describe-vpcs \
    --filters "Name=isDefault,Values=true" \
    --query 'Vpcs[0].VpcId' \
    --output text \
    --region $REGION)

if [ "$VPC_ID" == "None" ] || [ -z "$VPC_ID" ]; then
    echo "Error: No default VPC found. Please specify a VPC ID."
    exit 1
fi

echo "✅ Found VPC: $VPC_ID"

# Get or create security group
echo "Step 2: Setting up security group..."
SG_ID=$(aws ec2 describe-security-groups \
    --filters "Name=group-name,Values=$SECURITY_GROUP_NAME" "Name=vpc-id,Values=$VPC_ID" \
    --query 'SecurityGroups[0].GroupId' \
    --output text \
    --region $REGION)

if [ "$SG_ID" == "None" ] || [ -z "$SG_ID" ]; then
    echo "Creating security group..."
    SG_ID=$(aws ec2 create-security-group \
        --group-name $SECURITY_GROUP_NAME \
        --description "Security group for Helix.AI bioinformatics tools" \
        --vpc-id $VPC_ID \
        --query 'GroupId' \
        --output text \
        --region $REGION)
    
    # Add SSH access (restrict to your IP in production!)
    MY_IP=$(curl -s https://checkip.amazonaws.com)
    aws ec2 authorize-security-group-ingress \
        --group-id $SG_ID \
        --protocol tcp \
        --port 22 \
        --cidr "$MY_IP/32" \
        --region $REGION
    
    echo "✅ Created security group: $SG_ID"
    echo "⚠️  WARNING: SSH access is currently open to your IP only. Adjust as needed for production."
else
    echo "✅ Using existing security group: $SG_ID"
fi

# Get default AMI (Amazon Linux 2)
if [ -z "$AMI_ID" ]; then
    echo "Step 3: Finding Amazon Linux 2 AMI..."
    AMI_ID=$(aws ec2 describe-images \
        --owners amazon \
        --filters "Name=name,Values=amzn2-ami-hvm-*-x86_64-gp2" "Name=state,Values=available" \
        --query 'Images | sort_by(@, &CreationDate) | [-1].ImageId' \
        --output text \
        --region $REGION)
fi

echo "✅ Using AMI: $AMI_ID"

# Create user data script
USER_DATA=$(cat <<'EOF'
#!/bin/bash
# Install bioinformatics tools on instance startup

# Update system
yum update -y

# Install basic dependencies
yum install -y wget curl git gcc gcc-c++ make

# Install Miniconda
cd /tmp
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
rm -f Miniconda3-latest-Linux-x86_64.sh

# Add conda to PATH
export PATH="/opt/conda/bin:$PATH"

# Initialize conda
/opt/conda/bin/conda init bash

# Add bioconda channels
/opt/conda/bin/conda config --add channels bioconda
/opt/conda/bin/conda config --add channels conda-forge
/opt/conda/bin/conda config --set channel_priority strict

# Install bioinformatics tools (this takes 10-15 minutes)
/opt/conda/bin/conda install -y bbtools samtools bcftools bwa bowtie2 fastqc biopython pandas numpy

# Create symlinks
ln -sf /opt/conda/bin/bbmerge.sh /usr/local/bin/bbmerge.sh
ln -sf /opt/conda/bin/samtools /usr/local/bin/samtools
ln -sf /opt/conda/bin/bcftools /usr/local/bin/bcftools
ln -sf /opt/conda/bin/bwa /usr/local/bin/bwa
ln -sf /opt/conda/bin/bowtie2 /usr/local/bin/bowtie2
ln -sf /opt/conda/bin/fastqc /usr/local/bin/fastqc

# Create working directory
mkdir -p /opt/helix-tools
chmod 777 /opt/helix-tools

# Signal completion
echo "Bioinformatics tools installation completed" > /opt/helix-tools/install.log
EOF
)

# Base64 encode user data
USER_DATA_B64=$(echo "$USER_DATA" | base64)

# Launch instance
echo "Step 4: Launching EC2 instance..."
LAUNCH_ARGS=(
    --image-id "$AMI_ID"
    --instance-type "$INSTANCE_TYPE"
    --security-group-ids "$SG_ID"
    --user-data "$USER_DATA_B64"
    --tag-specifications "ResourceType=instance,Tags=[{Key=Project,Value=Helix.AI},{Key=Purpose,Value=BioinformaticsTools},{Key=Name,Value=Helix.AI-Bioinformatics-Tools}]"
    --query 'Instances[0].InstanceId'
    --output text
    --region $REGION
)

if [ -n "$KEY_NAME" ]; then
    LAUNCH_ARGS+=(--key-name "$KEY_NAME")
fi

INSTANCE_ID=$(aws ec2 run-instances "${LAUNCH_ARGS[@]}")

if [ -z "$INSTANCE_ID" ] || [ "$INSTANCE_ID" == "None" ]; then
    echo "Error: Failed to launch instance"
    exit 1
fi

echo "✅ Instance launched: $INSTANCE_ID"

# Wait for instance to be running
echo "Step 5: Waiting for instance to be running..."
aws ec2 wait instance-running --instance-ids "$INSTANCE_ID" --region $REGION
echo "✅ Instance is running"

# Get public IP
echo "Step 6: Getting instance details..."
PUBLIC_IP=$(aws ec2 describe-instances \
    --instance-ids "$INSTANCE_ID" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' \
    --output text \
    --region $REGION)

echo "✅ Public IP: $PUBLIC_IP"

echo ""
echo "=========================================="
echo "✅ EC2 Instance Setup Complete!"
echo "=========================================="
echo ""
echo "Instance ID: $INSTANCE_ID"
echo "Public IP: $PUBLIC_IP"
echo ""
echo "⚠️  IMPORTANT: Tools are being installed via user data script."
echo "   This takes 10-15 minutes. Wait before using the instance."
echo ""
echo "To check installation status, SSH to the instance and run:"
echo "  cat /opt/helix-tools/install.log"
echo ""
echo "To use this instance with Helix.AI, set these environment variables:"
echo "  export HELIX_EC2_INSTANCE_ID=$INSTANCE_ID"
echo "  export HELIX_EC2_KEY_NAME=$KEY_NAME"
echo "  export HELIX_EC2_KEY_FILE=~/.ssh/${KEY_NAME}.pem"
echo "  export HELIX_USE_EC2=true"
echo ""
echo "To SSH to the instance:"
echo "  ssh -i ~/.ssh/${KEY_NAME}.pem ec2-user@$PUBLIC_IP"
echo ""

