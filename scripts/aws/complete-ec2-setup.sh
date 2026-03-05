#!/bin/bash
# Complete EC2 setup: Install Python packages and bioinformatics tools, create symlinks

set -e

INSTANCE_ID="${1:-${HELIX_EC2_INSTANCE_ID}}"
KEY_FILE="${2:-${HELIX_EC2_KEY_FILE}}"
REGION="${AWS_REGION:-us-west-1}"

if [ -z "$INSTANCE_ID" ] || [ -z "$KEY_FILE" ]; then
    echo "Usage: $0 <instance-id> <key-file>"
    echo "   or set HELIX_EC2_INSTANCE_ID and HELIX_EC2_KEY_FILE"
    exit 1
fi

PUBLIC_IP=$(aws ec2 describe-instances \
    --instance-ids "$INSTANCE_ID" \
    --region "$REGION" \
    --query 'Reservations[0].Instances[0].PublicIpAddress' \
    --output text)

if [ "$PUBLIC_IP" == "None" ] || [ -z "$PUBLIC_IP" ]; then
    echo "❌ Could not get public IP for instance $INSTANCE_ID"
    exit 1
fi

echo "Completing EC2 setup..."
echo "Instance ID: $INSTANCE_ID"
echo "Public IP: $PUBLIC_IP"
echo ""

# Create setup script
cat > /tmp/complete_setup.sh << 'SETUP_SCRIPT'
#!/bin/bash
set -e

export PATH="/opt/conda/bin:$PATH"

# Increase resource limits to avoid thread creation errors
ulimit -u 4096 2>/dev/null || true

echo "Step 1: Installing Python packages via pip (more reliable)..."
/opt/conda/bin/pip install --no-cache-dir boto3 botocore biopython pandas numpy

echo ""
echo "Step 2: Installing bioinformatics tools via conda (one at a time to avoid resource issues)..."
# Install tools one at a time to avoid threading issues
for tool in bbtools samtools bcftools bwa bowtie2 fastqc; do
    echo "Installing $tool..."
    /opt/conda/bin/conda install -y -c bioconda $tool || echo "Warning: Failed to install $tool via conda"
done

echo ""
echo "Step 3: Creating symlinks..."
sudo ln -sf /opt/conda/bin/bbmerge.sh /usr/local/bin/bbmerge.sh || true
sudo ln -sf /opt/conda/bin/samtools /usr/local/bin/samtools || true
sudo ln -sf /opt/conda/bin/bcftools /usr/local/bin/bcftools || true
sudo ln -sf /opt/conda/bin/bwa /usr/local/bin/bwa || true
sudo ln -sf /opt/conda/bin/bowtie2 /usr/local/bin/bowtie2 || true
sudo ln -sf /opt/conda/bin/fastqc /usr/local/bin/fastqc || true

echo ""
echo "Step 4: Verifying installation..."
echo "Python: $(/opt/conda/bin/python3 --version)"
echo "boto3: $(/opt/conda/bin/python3 -c 'import boto3; print(boto3.__version__)' 2>/dev/null || echo 'not installed')"
echo "bbmerge.sh: $(which bbmerge.sh || echo 'not found')"
echo "samtools: $(which samtools || echo 'not found')"

echo ""
echo "✅ Setup complete!"
SETUP_SCRIPT

# Copy and execute
echo "Copying setup script to instance..."
scp -i "$KEY_FILE" -o StrictHostKeyChecking=no /tmp/complete_setup.sh ec2-user@$PUBLIC_IP:/tmp/

echo "Executing setup (this will take 10-15 minutes)..."
ssh -i "$KEY_FILE" -o StrictHostKeyChecking=no ec2-user@$PUBLIC_IP "bash /tmp/complete_setup.sh"

# Cleanup
rm -f /tmp/complete_setup.sh

echo ""
echo "✅ Setup complete! Run ./scripts/aws/run-ec2-test.sh to verify."

