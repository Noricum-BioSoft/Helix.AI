# EC2 Quick Start Guide

## Prerequisites

1. **AWS Account** with EC2 permissions
2. **EC2 Key Pair** - Create one if you don't have it:
   ```bash
   aws ec2 create-key-pair --key-name helix-bioinformatics --query 'KeyMaterial' --output text > ~/.ssh/helix-bioinformatics.pem
   chmod 600 ~/.ssh/helix-bioinformatics.pem
   ```
3. **AWS CLI configured** with credentials

## Option 1: Automated Setup (Recommended)

Run the test script to check your setup and optionally create an instance:

```bash
# Set required environment variables
export HELIX_EC2_KEY_NAME=helix-bioinformatics  # Your EC2 key pair name
export HELIX_EC2_KEY_FILE=~/.ssh/helix-bioinformatics.pem  # Path to your private key

# Optional: Set region
export AWS_REGION=us-east-1

# Run the test script
python scripts/aws/test-ec2-setup.py
```

The script will:
- ✅ Check all required configuration
- ✅ Verify AWS credentials
- ✅ Check for existing instances
- ✅ Create a new instance if needed
- ✅ Provide next steps

## Option 2: Setup Existing Instance (If you already have an instance)

If you have an existing EC2 instance that needs conda and IAM role setup:

```bash
# Set your instance details
export HELIX_EC2_INSTANCE_ID=i-0cc452521209457d4
export HELIX_EC2_KEY_FILE=/Users/eoberortner/.ssh/helix-ai-ec2.pem
export AWS_REGION=us-west-1

# Run the complete setup script
./scripts/aws/setup-existing-instance.sh
```

This will:
- ✅ Create and attach an IAM role with S3 permissions
- ✅ Install conda with Python 3.8+
- ✅ Install bioinformatics tools (BBTools, samtools, etc.)
- ✅ Install Python packages (boto3, biopython, etc.)
- ⏱️ **Takes 10-15 minutes** - the script shows progress

## Option 3: Manual Setup with Bash Script

```bash
# Set environment variables
export EC2_KEY_NAME=helix-bioinformatics
export AWS_REGION=us-east-1
export INSTANCE_TYPE=t3.medium

# Run the setup script
cd scripts/aws
./create-bioinformatics-ec2.sh
```

This will:
- Create a security group
- Launch an EC2 instance
- Install bioinformatics tools (takes 10-15 minutes)
- Return the instance ID

## Option 4: Let the System Auto-Create

If you just want the system to automatically create an instance when needed:

```bash
# Set required variables
export HELIX_EC2_KEY_NAME=helix-bioinformatics
export HELIX_EC2_KEY_FILE=~/.ssh/helix-bioinformatics.pem
export HELIX_EC2_AUTO_CREATE=true  # or '1', 'yes', 'on', 'auto-create'
export HELIX_USE_EC2=true

# Don't set HELIX_EC2_INSTANCE_ID - let it auto-create
```

## After Setup

1. **Wait for tools to install** (10-15 minutes after instance launch)
   ```bash
   # SSH to check status
   ssh -i ~/.ssh/helix-bioinformatics.pem ec2-user@<instance-ip>
   cat /opt/helix-tools/install.log
   ```

2. **Set environment variables** for your backend:
   ```bash
   export HELIX_EC2_INSTANCE_ID=i-xxxxxxxxx  # From the setup output
   export HELIX_EC2_KEY_NAME=helix-bioinformatics
   export HELIX_EC2_KEY_FILE=~/.ssh/helix-bioinformatics.pem
   export HELIX_USE_EC2=true
   ```

3. **Test your setup** (recommended):
   ```bash
   # Run the comprehensive test suite
   ./scripts/aws/run-ec2-test.sh
   ```
   This will test:
   - Python version and packages
   - AWS credentials and S3 access
   - Bioinformatics tools availability
   - Tool execution
   - S3 upload/download operations

4. **Restart your backend server**

## Testing Your Setup

After setting up your EC2 instance, verify everything works:

```bash
# Set your instance details
export HELIX_EC2_INSTANCE_ID=i-0cc452521209457d4
export HELIX_EC2_KEY_FILE=/Users/eoberortner/.ssh/helix-ai-ec2.pem
export AWS_REGION=us-west-1

# Run the test suite
./scripts/aws/run-ec2-test.sh
```

The test suite verifies:
- ✅ Python version (should be 3.8+)
- ✅ Python packages (boto3, biopython, pandas, numpy)
- ✅ AWS credentials and S3 access
- ✅ Bioinformatics tools (samtools, bbmerge, bwa, etc.)
- ✅ Tool execution
- ✅ S3 upload/download/delete operations

## Troubleshooting

### "InvalidInstanceID.NotFound"
- The instance ID doesn't exist
- **Solution**: Set `HELIX_EC2_AUTO_CREATE=true` to auto-create, or use a valid instance ID

### "Key pair not found"
- The key pair name doesn't exist in your AWS region
- **Solution**: Create a key pair or use an existing one:
  ```bash
  aws ec2 describe-key-pairs --region us-east-1
  ```

### "Cannot create EC2 instance: HELIX_EC2_KEY_NAME is not set"
- Missing required configuration
- **Solution**: Set `HELIX_EC2_KEY_NAME` environment variable

### Tools not available on instance
- Tools are still installing (takes 10-15 minutes)
- **Solution**: Wait and check `/opt/helix-tools/install.log` on the instance

## Cost Estimate

- **t3.medium**: ~$0.04/hour (~$30/month if running 24/7)
- **t3.large**: ~$0.08/hour (~$60/month if running 24/7)
- Consider stopping the instance when not in use to save costs

