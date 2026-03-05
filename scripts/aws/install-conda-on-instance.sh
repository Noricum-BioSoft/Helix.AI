#!/bin/bash
# Install conda and bioinformatics tools on an existing EC2 instance

set -e

INSTANCE_ID="${1:-${HELIX_EC2_INSTANCE_ID}}"
KEY_FILE="${2:-${HELIX_EC2_KEY_FILE}}"
REGION="${AWS_REGION:-us-west-1}"

if [ -z "$INSTANCE_ID" ] || [ -z "$KEY_FILE" ]; then
    echo "Usage: $0 <instance-id> <key-file>"
    echo "   or set HELIX_EC2_INSTANCE_ID and HELIX_EC2_KEY_FILE"
    exit 1
fi

# Get instance IP
PUBLIC_IP=$(aws ec2 describe-instances \
    --instance-ids "$INSTANCE_ID" \
    --region "$REGION" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' \
    --output text)

if [ "$PUBLIC_IP" == "None" ] || [ -z "$PUBLIC_IP" ]; then
    echo "❌ Could not get public IP for instance $INSTANCE_ID"
    exit 1
fi

echo "Installing conda and bioinformatics tools on instance..."
echo "Instance ID: $INSTANCE_ID"
echo "Public IP: $PUBLIC_IP"
echo ""

# Create installation script
cat > /tmp/install_conda.sh << 'INSTALL_SCRIPT'
#!/bin/bash
set -e

echo "Installing Miniconda and bioinformatics tools..."

# Update system
sudo yum update -y

# Install basic dependencies
sudo yum install -y wget curl git gcc gcc-c++ make

# Install Miniconda (using version compatible with GLIBC 2.26)
cd /tmp
# Use Miniconda3-py38_4.9.2 which works with GLIBC 2.26
wget -q https://repo.anaconda.com/miniconda/Miniconda3-py38_4.9.2-Linux-x86_64.sh
chmod +x Miniconda3-py38_4.9.2-Linux-x86_64.sh
sudo bash Miniconda3-py38_4.9.2-Linux-x86_64.sh -b -p /opt/conda
rm -f Miniconda3-py38_4.9.2-Linux-x86_64.sh

# Add conda to PATH
export PATH="/opt/conda/bin:$PATH"

# Initialize conda
/opt/conda/bin/conda init bash

# Add bioconda channels
/opt/conda/bin/conda config --add channels bioconda
/opt/conda/bin/conda config --add channels conda-forge
/opt/conda/bin/conda config --set channel_priority strict

# Install common bioinformatics tools and Python packages
/opt/conda/bin/conda install -y bbtools samtools bcftools bwa bowtie2 fastqc biopython pandas numpy boto3 botocore

# Create symlinks for easy access
sudo ln -sf /opt/conda/bin/bbmerge.sh /usr/local/bin/bbmerge.sh || true
sudo ln -sf /opt/conda/bin/samtools /usr/local/bin/samtools || true
sudo ln -sf /opt/conda/bin/bcftools /usr/local/bin/bcftools || true
sudo ln -sf /opt/conda/bin/bwa /usr/local/bin/bwa || true
sudo ln -sf /opt/conda/bin/bowtie2 /usr/local/bin/bowtie2 || true
sudo ln -sf /opt/conda/bin/fastqc /usr/local/bin/fastqc || true

# Create working directory
sudo mkdir -p /opt/helix-tools
sudo chmod 777 /opt/helix-tools

# Signal completion
echo "Bioinformatics tools installation completed" | sudo tee /opt/helix-tools/install.log

echo ""
echo "✅ Installation complete!"
echo "Python version: $(/opt/conda/bin/python3 --version)"
INSTALL_SCRIPT

# Copy and execute installation script
echo "Copying installation script to instance..."
scp -i "$KEY_FILE" -o StrictHostKeyChecking=no /tmp/install_conda.sh ec2-user@$PUBLIC_IP:/tmp/

echo "Executing installation script (this will take 10-15 minutes)..."
ssh -i "$KEY_FILE" -o StrictHostKeyChecking=no ec2-user@$PUBLIC_IP "bash /tmp/install_conda.sh"

echo ""
echo "✅ Installation complete!"
echo ""
echo "Verifying installation..."
ssh -i "$KEY_FILE" -o StrictHostKeyChecking=no ec2-user@$PUBLIC_IP "/opt/conda/bin/python3 --version && /opt/conda/bin/python3 -c 'import boto3; print(\"boto3 version:\", boto3.__version__)'"

# Cleanup
rm -f /tmp/install_conda.sh

