# EMR Cluster Setup Guide for FASTQ Processing

## Quick Start

### 1. Run Setup Script

```bash
# Basic setup (uses defaults)
./scripts/aws/setup-emr-cluster.sh

# With custom configuration
export EMR_CLUSTER_NAME="my-fastqc-cluster"
export EMR_INSTANCE_TYPE="m5.2xlarge"
export EMR_INSTANCE_COUNT=5
export EC2_KEY_NAME="my-key-pair"  # Optional, for SSH access
./scripts/aws/setup-emr-cluster.sh
```

### 2. Monitor Cluster Creation

```bash
# Get cluster ID from output, then check status
CLUSTER_ID="j-XXXXXXXXXXXXX"
aws emr describe-cluster --cluster-id $CLUSTER_ID --query 'Cluster.Status.State' --output text

# Or use helper function
source scripts/aws/emr-commands.sh
emr_status $CLUSTER_ID
```

### 3. Wait for Cluster to be Ready

Cluster status will progress: `STARTING` → `BOOTSTRAPPING` → `RUNNING`

Typically takes 5-10 minutes.

---

## Manual Setup (Step by Step)

### Prerequisites

1. **AWS CLI configured** with appropriate credentials
2. **S3 bucket** for storing bootstrap scripts and logs
3. **IAM permissions** for EMR, EC2, S3, IAM

### Step 1: Create IAM Roles

The setup script does this automatically, but here's the manual process:

#### EMR Service Role

```bash
# Create trust policy
cat > emr-service-trust.json << 'EOF'
{
  "Version": "2012-10-17",
  "Statement": [{
    "Effect": "Allow",
    "Principal": {"Service": "elasticmapreduce.amazonaws.com"},
    "Action": "sts:AssumeRole"
  }]
}
EOF

# Create role
aws iam create-role \
  --role-name EMR_DefaultRole \
  --assume-role-policy-document file://emr-service-trust.json

# Attach managed policy
aws iam attach-role-policy \
  --role-name EMR_DefaultRole \
  --policy-arn arn:aws:iam::aws:policy/service-role/AmazonElasticMapReduceRole
```

#### EMR EC2 Instance Profile

```bash
# Create EC2 role
aws iam create-role \
  --role-name EMR_EC2_DefaultRole \
  --assume-role-policy-document file://emr-ec2-trust.json

# Attach managed policy
aws iam attach-role-policy \
  --role-name EMR_EC2_DefaultRole \
  --policy-arn arn:aws:iam::aws:policy/service-role/AmazonElasticMapReduceforEC2Role

# Create instance profile
aws iam create-instance-profile --instance-profile-name EMR_EC2_DefaultRole
aws iam add-role-to-instance-profile \
  --instance-profile-name EMR_EC2_DefaultRole \
  --role-name EMR_EC2_DefaultRole
```

### Step 2: Create Bootstrap Script

```bash
# Create bootstrap script
cat > fastqc-bootstrap.sh << 'EOF'
#!/bin/bash
sudo pip3 install biopython pandas numpy matplotlib seaborn plotly
sudo mkdir -p /mnt/var/log/fastqc
sudo chmod 777 /mnt/var/log/fastqc
EOF

# Upload to S3
aws s3 cp fastqc-bootstrap.sh s3://noricum-ngs-data/emr-bootstrap/fastqc-bootstrap.sh
```

### Step 3: Launch EMR Cluster

```bash
aws emr create-cluster \
  --name "helix-ai-fastqc" \
  --release-label emr-6.15.0 \
  --applications Name=Spark Name=Hadoop Name=Zeppelin \
  --instance-groups \
    InstanceGroupType=MASTER,InstanceCount=1,InstanceType=m5.xlarge \
    InstanceGroupType=CORE,InstanceCount=2,InstanceType=m5.xlarge,MarketType=SPOT,BidPrice=0.10 \
  --bootstrap-actions Path=s3://noricum-ngs-data/emr-bootstrap/fastqc-bootstrap.sh,Name="Install Dependencies" \
  --service-role EMR_DefaultRole \
  --job-flow-role EMR_EC2_DefaultRole \
  --log-uri s3://noricum-ngs-data/emr-logs/ \
  --region us-east-1 \
  --tags Project=Helix.AI Purpose="FASTQ Quality Analysis"
```

---

## Cluster Configuration Options

### Instance Types

| Type | vCPU | RAM | Best For |
|------|------|-----|----------|
| `m5.xlarge` | 4 | 16GB | Small-medium datasets |
| `m5.2xlarge` | 8 | 32GB | Medium-large datasets |
| `m5.4xlarge` | 16 | 64GB | Large datasets |
| `r5.2xlarge` | 8 | 64GB | Memory-intensive |
| `c5.2xlarge` | 8 | 16GB | CPU-intensive |

### Cost Optimization

1. **Use Spot Instances** for core nodes (save 50-90%)
   ```bash
   MarketType=SPOT,BidPrice=0.10
   ```

2. **Auto-terminate** when idle
   ```bash
   --auto-terminate
   ```

3. **Use smaller instances** for development
   ```bash
   InstanceType=m5.large
   ```

4. **Terminate immediately** after job completes
   ```bash
   --auto-terminate
   ```

### Spark Configuration

Add to cluster creation:

```bash
--configurations '[
  {
    "Classification": "spark-defaults",
    "Properties": {
      "spark.sql.adaptive.enabled": "true",
      "spark.sql.adaptive.coalescePartitions.enabled": "true",
      "spark.serializer": "org.apache.spark.serializer.KryoSerializer",
      "spark.executor.memory": "8g",
      "spark.executor.cores": "4",
      "spark.driver.memory": "4g"
    }
  }
]'
```

---

## Common Operations

### List Clusters

```bash
aws emr list-clusters --active --output table
```

### Get Cluster Status

```bash
CLUSTER_ID="j-XXXXXXXXXXXXX"
aws emr describe-cluster --cluster-id $CLUSTER_ID \
  --query 'Cluster.[Status.State,Status.StateChangeReason.Message]' \
  --output table
```

### Get Master Node DNS (for SSH)

```bash
aws emr describe-cluster --cluster-id $CLUSTER_ID \
  --query 'Cluster.MasterPublicDnsName' \
  --output text
```

### SSH into Cluster

```bash
MASTER_DNS=$(aws emr describe-cluster --cluster-id $CLUSTER_ID \
  --query 'Cluster.MasterPublicDnsName' --output text)

ssh -i ~/.ssh/your-key.pem hadoop@$MASTER_DNS
```

### Submit Spark Job

```bash
aws emr add-steps \
  --cluster-id $CLUSTER_ID \
  --steps Type=SPARK,Name="FASTQ Analysis",ActionOnFailure=CONTINUE,\
Args="[--deploy-mode,cluster,--py-files,s3://bucket/scripts/fastqc.py,s3://bucket/input/*.fq]" \
  --query 'StepIds[0]' \
  --output text
```

### Check Step Status

```bash
STEP_ID="s-XXXXXXXXXXXXX"
aws emr describe-step \
  --cluster-id $CLUSTER_ID \
  --step-id $STEP_ID \
  --query 'Step.Status.State' \
  --output text
```

### Terminate Cluster

```bash
aws emr terminate-clusters --cluster-ids $CLUSTER_ID
```

---

## Using Helper Scripts

Source the helper functions:

```bash
source scripts/aws/emr-commands.sh

# List clusters
emr_list

# Check status
emr_status j-XXXXXXXXXXXXX

# Get master DNS
emr_master_dns j-XXXXXXXXXXXXX

# Terminate cluster
emr_terminate j-XXXXXXXXXXXXX
```

---

## Integration with Helix.AI

### Backend API Endpoint

Create endpoint to submit EMR jobs:

```python
@app.post("/emr/submit-fastqc")
async def submit_fastqc_job(
    session_id: str,
    dataset_id: str,
    r1_path: str,
    r2_path: str
):
    """Submit FASTQ quality analysis job to EMR."""
    # 1. Get or create EMR cluster
    # 2. Submit Spark job
    # 3. Return job ID
    # 4. Poll for completion
    pass
```

### Job Status Endpoint

```python
@app.get("/emr/job/{job_id}/status")
async def get_job_status(job_id: str):
    """Get EMR job status."""
    # Query EMR step status
    # Return progress, results location
    pass
```

---

## Troubleshooting

### Cluster Fails to Start

1. **Check IAM roles**: Ensure service role and instance profile exist
2. **Check S3 permissions**: EC2 role needs S3 access
3. **Check subnet**: May need to specify VPC subnet
4. **Check limits**: AWS account may have instance limits

### Bootstrap Script Fails

1. **Check S3 path**: Bootstrap script must be accessible
2. **Check script permissions**: Must be executable
3. **View logs**: Check `/mnt/var/log/bootstrap-actions/` on master node

### Job Fails

1. **Check Spark logs**: `yarn logs -applicationId <app-id>`
2. **Check S3 access**: Job needs read/write access to S3
3. **Check Python dependencies**: May need to install in bootstrap

### High Costs

1. **Use Spot instances**: Can save 50-90%
2. **Terminate when done**: Use `--auto-terminate`
3. **Use smaller instances**: Start with m5.large for testing
4. **Monitor usage**: Set up CloudWatch alarms

---

## Cost Estimates

### Small Cluster (3 m5.xlarge)

- **On-Demand**: ~$0.50/hour
- **With Spot (core nodes)**: ~$0.25/hour
- **Typical job (1 hour)**: $0.25-0.50

### Medium Cluster (5 m5.2xlarge)

- **On-Demand**: ~$2.00/hour
- **With Spot**: ~$1.00/hour
- **Typical job (1 hour)**: $1.00-2.00

### Cost Optimization Tips

1. Use Spot instances for 50-90% savings
2. Auto-terminate clusters after jobs complete
3. Use smaller instances for development
4. Monitor with CloudWatch to avoid runaway costs

---

## Next Steps

1. **Run setup script**: `./scripts/aws/setup-emr-cluster.sh`
2. **Wait for cluster**: Monitor until status is "RUNNING"
3. **Submit test job**: Use helper functions or AWS CLI
4. **Integrate with backend**: Create API endpoints for job submission
5. **Monitor costs**: Set up CloudWatch alarms










