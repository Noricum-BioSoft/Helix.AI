# EC2 Execution Setup for Tool Generator Agent

This guide explains how to set up an EC2 instance with pre-installed bioinformatics tools for executing generated code.

## Overview

The tool-generator-agent can execute generated code on a dedicated EC2 instance that has bioinformatics tools pre-installed. This ensures that tools like BBMerge, samtools, bcftools, etc. are always available without needing installation.

## Quick Start

### Option 1: Use Existing EC2 Instance

If you already have an EC2 instance with bioinformatics tools:

```bash
export HELIX_EC2_INSTANCE_ID=i-1234567890abcdef0
export HELIX_EC2_KEY_NAME=my-key-pair
export HELIX_EC2_KEY_FILE=~/.ssh/my-key-pair.pem
export HELIX_USE_EC2=true
```

### Option 2: Create New EC2 Instance

1. **Create the instance with tools pre-installed**:
   ```bash
   cd scripts/aws
   ./create-bioinformatics-ec2.sh
   ```

2. **Wait for tools to install** (10-15 minutes):
   ```bash
   # SSH to the instance
   ssh -i ~/.ssh/your-key.pem ec2-user@<public-ip>
   
   # Check installation status
   cat /opt/helix-tools/install.log
   ```

3. **Set environment variables**:
   ```bash
   export HELIX_EC2_INSTANCE_ID=<instance-id-from-script>
   export HELIX_EC2_KEY_NAME=<your-key-name>
   export HELIX_EC2_KEY_FILE=~/.ssh/<your-key-name>.pem
   export HELIX_USE_EC2=true
   ```

### Option 3: Install Tools on Existing Instance

If you have an existing EC2 instance:

1. **SSH to the instance**
2. **Run the setup script**:
   ```bash
   curl -s https://raw.githubusercontent.com/your-repo/Helix.AI/main/scripts/aws/setup-bioinformatics-ec2.sh | bash
   ```
   
   Or copy and run locally:
   ```bash
   scp scripts/aws/setup-bioinformatics-ec2.sh ec2-user@<instance-ip>:/tmp/
   ssh ec2-user@<instance-ip> "bash /tmp/setup-bioinformatics-ec2.sh"
   ```

3. **Set environment variables** (as in Option 1)

## Environment Variables

| Variable | Description | Required |
|----------|-------------|----------|
| `HELIX_USE_EC2` | Enable EC2 execution (`true`/`false`) | Yes |
| `HELIX_EC2_INSTANCE_ID` | EC2 instance ID to use | No (auto-created if not set) |
| `HELIX_EC2_KEY_NAME` | EC2 key pair name | Yes |
| `HELIX_EC2_KEY_FILE` | Path to SSH private key file | Yes |
| `HELIX_EC2_AUTO_CREATE` | Auto-create instance if configured ID doesn't exist (`true`, `1`, `yes`, `on`, `auto-create`, or any non-empty value except `false`) | No |
| `HELIX_EC2_INSTANCE_TYPE` | Instance type for new instances (default: `t3.medium`) | No |
| `HELIX_EC2_AMI_ID` | Custom AMI ID (default: Amazon Linux 2) | No |
| `HELIX_EC2_SECURITY_GROUP_ID` | Security group ID (auto-created if not provided) | No |
| `HELIX_EC2_SUBNET_ID` | Subnet ID (uses default VPC if not provided) | No |

## Installed Tools

The EC2 instance includes:

- **Read Processing**: BBTools (bbmerge.sh), FLASH, PEAR
- **Sequence Analysis**: samtools, bcftools, BWA, Bowtie2
- **Quality Control**: FastQC
- **Python Libraries**: BioPython, pandas, numpy

All tools are installed via conda/bioconda and available in `/opt/conda/bin/` with symlinks in `/usr/local/bin/`.

## How It Works

1. **Code Generation**: The tool-generator-agent generates Python code
2. **Execution Decision**: If `HELIX_USE_EC2=true`, code executes on EC2
3. **Instance Management**: 
   - If `HELIX_EC2_INSTANCE_ID` is set:
     - Uses that instance if it exists and is running
     - If instance doesn't exist and `HELIX_EC2_AUTO_CREATE=true`, creates a new instance
     - If instance is stopped, attempts to start it
     - Otherwise, falls back to local execution with a warning
   - If `HELIX_EC2_INSTANCE_ID` is not set:
     - Finds existing tagged instance (with `Project=Helix.AI` and `Purpose=BioinformaticsTools`)
     - If none found, creates a new instance
4. **Code Execution**: 
   - Code is copied to EC2 via SCP
   - Executed via SSH with Python from `/opt/conda/bin/python`
   - Results are returned to the backend

## Security Considerations

1. **SSH Key Security**: Keep your private key file secure (chmod 600)
2. **Security Groups**: Restrict SSH access to your IP in production
3. **IAM Permissions**: 
   - **EC2 executor** needs permissions to:
     - `ec2:DescribeInstances`
     - `ec2:RunInstances` (if creating new instances)
     - `ec2:StartInstances`
     - `ec2:CreateSecurityGroup` (if creating security groups)
     - `ec2:AuthorizeSecurityGroupIngress`
   - **EC2 Instance** needs an IAM role with S3 permissions:
     - `s3:GetObject` (to read from S3)
     - `s3:PutObject` (to write to S3)
     - `s3:ListBucket` (to list bucket contents)
     - The instance will automatically use IAM role credentials - no manual credential configuration needed

## Troubleshooting

### Instance Not Found

If you see `InvalidInstanceID.NotFound` error:

1. **Verify the instance ID**:
   ```bash
   aws ec2 describe-instances --instance-ids i-xxxxxxxxx
   ```

2. **Auto-create option**: Set `HELIX_EC2_AUTO_CREATE=true` to automatically create a new instance if the configured one doesn't exist:
   ```bash
   export HELIX_EC2_AUTO_CREATE=true
   ```

3. **Manual fix**: Update `HELIX_EC2_INSTANCE_ID` with a valid instance ID, or remove it to let the system find/create one automatically

4. **Region check**: Verify instance is in the same region as configured (`AWS_REGION` env var or default `us-east-1`)

5. **Permissions**: Check AWS credentials have permission to describe instances

### SSH Connection Failed
- Verify security group allows SSH from your IP
- Check that key file path is correct and has correct permissions (chmod 600)
- Ensure instance has a public IP address

### Tools Not Available
- Wait for user data script to complete (10-15 minutes after instance launch)
- SSH to instance and check: `cat /opt/helix-tools/install.log`
- Manually run setup script if needed: `bash scripts/aws/setup-bioinformatics-ec2.sh`

### Execution Timeout
- Increase timeout in `ec2_executor.py` if needed
- Check instance has sufficient resources (CPU/memory)
- Verify network connectivity between backend and EC2

## Cost Considerations

- **t3.medium**: ~$0.04/hour (~$30/month if running 24/7)
- **t3.large**: ~$0.08/hour (~$60/month if running 24/7)
- Consider stopping instance when not in use to save costs
- Use spot instances for development to reduce costs

## Best Practices

1. **Use a dedicated instance** for production with tools pre-installed
2. **Create a custom AMI** after tools are installed to speed up new instance launches
3. **Tag instances** properly for easy management
4. **Monitor costs** and stop instances when not needed
5. **Use security groups** to restrict access appropriately
6. **Set up CloudWatch alarms** for instance health monitoring

