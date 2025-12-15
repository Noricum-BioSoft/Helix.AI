# EMR FASTQ Analysis Scripts

## Quick Start

### 1. Set Cluster ID

```bash
export EMR_CLUSTER_ID="j-12QYDE51Q9LDP"
```

### 2. Submit Job

```bash
# Use default paths (GRCh38.p12.MafHi dataset)
./scripts/emr/submit-fastqc-job.sh

# Or specify custom paths
./scripts/emr/submit-fastqc-job.sh \
  s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/mate_R1.fq \
  s3://noricum-ngs-data/datasets/GRCh38.p12.MafHi/mate_R2.fq \
  s3://noricum-ngs-data/fastqc-results/GRCh38.p12.MafHi/
```

### 3. Monitor Job

```bash
# Check step status
STEP_ID="s-XXXXXXXXXXXXX"
aws emr describe-step \
  --cluster-id $EMR_CLUSTER_ID \
  --step-id $STEP_ID \
  --region us-east-1 \
  --query 'Step.Status.State' \
  --output text

# Or use helper
source scripts/aws/emr-commands.sh
emr_step_status $EMR_CLUSTER_ID $STEP_ID
```

### 4. View Results

```bash
# List output files
aws s3 ls s3://noricum-ngs-data/fastqc-results/GRCh38.p12.MafHi/ --recursive

# Download results
aws s3 cp s3://noricum-ngs-data/fastqc-results/GRCh38.p12.MafHi/part-00000 ./results.json
```

## Job Status States

- `PENDING` - Job queued, waiting to start
- `RUNNING` - Job currently executing
- `COMPLETED` - Job finished successfully
- `CANCELLED` - Job was cancelled
- `FAILED` - Job failed (check logs)
- `INTERRUPTED` - Job was interrupted

## Viewing Logs

Logs are stored in S3:

```bash
# List log files
aws s3 ls s3://noricum-ngs-data/emr-logs/$EMR_CLUSTER_ID/steps/$STEP_ID/ --recursive

# Download stderr (error logs)
aws s3 cp s3://noricum-ngs-data/emr-logs/$EMR_CLUSTER_ID/steps/$STEP_ID/stderr.gz ./stderr.gz
gunzip stderr.gz
cat stderr

# Download stdout (output logs)
aws s3 cp s3://noricum-ngs-data/emr-logs/$EMR_CLUSTER_ID/steps/$STEP_ID/stdout.gz ./stdout.gz
gunzip stdout.gz
cat stdout
```

## Troubleshooting

### Job Fails Immediately

1. Check cluster state: `aws emr describe-cluster --cluster-id $CLUSTER_ID`
2. Check S3 permissions for input/output paths
3. Verify script uploaded correctly

### Job Takes Too Long

1. Check Spark UI (if accessible via SSH tunnel)
2. Monitor CloudWatch metrics
3. Consider increasing cluster size

### Results Not Appearing

1. Check job status (should be COMPLETED)
2. Verify output path permissions
3. Check for errors in stderr logs







