# Deployment Status

## ✅ Deployment Complete — Feb 22, 2026

All components are live and healthy.

---

## URLs

| Service | URL |
|---------|-----|
| **Frontend (CloudFront)** | https://d2a8mt5n89vos4.cloudfront.net |
| **Backend API (ALB)** | http://HelixA-ALBAE-D7BksiQIynZb-1051248867.us-west-1.elb.amazonaws.com |
| **Backend Health** | http://HelixA-ALBAE-D7BksiQIynZb-1051248867.us-west-1.elb.amazonaws.com/health |

---

## Infrastructure

| Component | Value | Status |
|-----------|-------|--------|
| CDK Stack | `HelixAIStack` | ✅ UPDATE_COMPLETE |
| EC2 Instance | `i-08f3780e06b932be4` (t3.medium, us-west-1) | ✅ Running |
| Docker | 25.0.14 | ✅ Installed |
| Backend Container | `helix-backend` on port 8001 | ✅ Healthy |
| ECR Image | `helix-ai-backend:latest` (Dec 26, 2025) | ✅ Deployed |
| ALB Target Group | EC2 target group | ✅ Healthy |
| S3 Frontend Bucket | `helix-ai-frontend-794270057041-us-west-1` | ✅ Deployed |
| CloudFront | `E1F1SKENKU7JKN` → `d2a8mt5n89vos4.cloudfront.net` | ✅ Deployed |
| Secrets Manager | `OPENAI_API_KEY`, `DEEPSEEK_API_KEY` | ✅ Configured |

---

## Re-deploying

### Backend only (new Docker image)
```bash
cd scripts/aws
./ecr_push_backend.sh 794270057041 us-west-1 helix-ai-backend latest
# Then pull on EC2:
aws ssm send-command \
  --instance-id i-08f3780e06b932be4 \
  --region us-west-1 \
  --document-name AWS-RunShellScript \
  --parameters '{"commands":[
    "aws ecr get-login-password --region us-west-1 | docker login --username AWS --password-stdin 794270057041.dkr.ecr.us-west-1.amazonaws.com",
    "docker pull 794270057041.dkr.ecr.us-west-1.amazonaws.com/helix-ai-backend:latest",
    "docker stop helix-backend && docker rm helix-backend",
    "OPENAI_KEY=$(aws secretsmanager get-secret-value --secret-id helix-ai-production-OPENAI_API_KEY --region us-west-1 --query SecretString --output text)",
    "DEEPSEEK_KEY=$(aws secretsmanager get-secret-value --secret-id helix-ai-production-DEEPSEEK_API_KEY --region us-west-1 --query SecretString --output text)",
    "docker run -d --name helix-backend --restart unless-stopped -p 8001:8001 -e OPENAI_API_KEY=$OPENAI_KEY -e DEEPSEEK_API_KEY=$DEEPSEEK_KEY -e AWS_REGION=us-west-1 -e AWS_DEFAULT_REGION=us-west-1 -e CORS_ORIGINS=* 794270057041.dkr.ecr.us-west-1.amazonaws.com/helix-ai-backend:latest"
  ]}'
```

### Frontend only
```bash
cd frontend
VITE_API_BASE_URL=http://HelixA-ALBAE-D7BksiQIynZb-1051248867.us-west-1.elb.amazonaws.com npm run build
aws s3 sync dist/ s3://helix-ai-frontend-794270057041-us-west-1/ --region us-west-1 --delete \
  --cache-control "public,max-age=31536000,immutable" --exclude "index.html"
aws s3 cp dist/index.html s3://helix-ai-frontend-794270057041-us-west-1/index.html \
  --region us-west-1 --cache-control "no-cache,no-store,must-revalidate" --content-type "text/html"
aws cloudfront create-invalidation --distribution-id E1F1SKENKU7JKN --paths "/*"
```

### Full deploy (backend + frontend)
```bash
cd scripts/aws
./deploy.sh deploy.config
```
