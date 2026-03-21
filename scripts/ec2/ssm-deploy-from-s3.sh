#!/bin/bash
# Run on EC2 as root (e.g. via SSM). Fetches tarball from S3 and runs bootstrap-ec2.sh
set -euxo pipefail
BUCKET="${1:-helix-ai-frontend-794270057041-us-west-1}"
KEY="${2:-bootstrap/helix-ai-source.tar.gz}"
cd /opt
mkdir -p helix
cd helix
aws s3 cp "s3://${BUCKET}/${KEY}" ./helix-ai-source.tar.gz
# git archive extracts at cwd (no top-level Helix.AI/ dir). Clear old tree except the tarball.
find . -mindepth 1 -maxdepth 1 ! -name 'helix-ai-source.tar.gz' -exec rm -rf {} +
tar xzf helix-ai-source.tar.gz
exec ./scripts/ec2/bootstrap-ec2.sh
