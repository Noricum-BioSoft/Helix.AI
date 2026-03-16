# Helix.AI Architecture: Frontend vs Backend Deployment

## Quick Answer

**In Production:**
- ✅ **Backend**: Runs on AWS (ECS Fargate) - **You have this deployed**
- ✅ **Frontend**: Should also run on AWS (S3 + CloudFront) - **Needs to be deployed**

**In Development:**
- Both can run locally (frontend on localhost:5173, backend on localhost:8001)

---

## Architecture Overview

### Production Architecture (Recommended)

```
┌─────────────────────────────────────────────────────────────┐
│                        Users/Internet                        │
└────────────────────┬────────────────────────────────────────┘
                     │
         ┌───────────┴───────────┐
         │                       │
         ▼                       ▼
┌─────────────────┐    ┌──────────────────────┐
│   CloudFront     │    │  Application Load   │
│   (CDN)          │    │  Balancer (ALB)      │
│                  │    │                      │
│  Frontend (S3)   │    │  Backend (ECS)       │
└─────────────────┘    └──────────────────────┘
         │                       │
         │                       │
         ▼                       ▼
┌─────────────────┐    ┌──────────────────────┐
│   S3 Bucket      │    │  ECS Fargate         │
│   (Static Files) │    │  (Container)         │
│                  │    │                      │
│  - index.html    │    │  - FastAPI App       │
│  - JS bundles    │    │  - Port 8001         │
│  - CSS, assets   │    │  - Health: /health   │
└─────────────────┘    └──────────────────────┘
```

### Current Status

| Component | Status | Location |
|-----------|--------|----------|
| **Backend** | ✅ Deployed | AWS ECS Fargate (via Copilot) |
| **Frontend** | ⚠️ Not deployed | Needs deployment to S3/CloudFront |

---

## How It Works

### Frontend (React App)

**What it is:**
- Static files (HTML, CSS, JavaScript)
- Built with `npm run build` → creates `frontend/dist/` folder
- No server needed - just files

**Where it runs:**
- **Development**: `npm run dev` → runs locally on `localhost:5173`
- **Production**: Deployed to S3 → served via CloudFront CDN

**How it connects to backend:**
- Frontend makes API calls to backend URL
- Set via `VITE_API_BASE_URL` at build time
- All API calls go to: `${API_BASE_URL}/health`, `${API_BASE_URL}/execute`, etc.

### Backend (FastAPI App)

**What it is:**
- Python application (FastAPI)
- Runs in Docker container
- Needs a server to run

**Where it runs:**
- **Development**: `uvicorn main_with_mcp:app` → runs locally on `localhost:8001`
- **Production**: ✅ Deployed to ECS Fargate → accessible via ALB

**Current URL:**
```
http://helix--Publi-7CO3kaz1lrRR-443122446.us-west-1.elb.amazonaws.com
```

---

## Development vs Production

### Development Mode

```bash
# Terminal 1: Backend
cd backend
uvicorn main_with_mcp:app --reload --port 8001

# Terminal 2: Frontend  
cd frontend
npm run dev
```

- Frontend: `http://localhost:5173`
- Backend: `http://localhost:8001`
- Frontend connects to `http://localhost:8001` (default)

### Production Mode (What You Should Have)

```
Frontend: https://your-cloudfront-domain.cloudfront.net
         ↓ (makes API calls to)
Backend:  http://helix--Publi-7CO3kaz1lrRR-443122446.us-west-1.elb.amazonaws.com
```

- Frontend: Served from S3/CloudFront (static files)
- Backend: Running on ECS Fargate (container)
- Frontend connects to ALB URL (set via `VITE_API_BASE_URL`)

---

## Why This Architecture?

### Frontend on S3/CloudFront
- ✅ **Cost-effective**: S3 is cheap for static files
- ✅ **Fast**: CloudFront CDN delivers files globally
- ✅ **Scalable**: Handles unlimited traffic
- ✅ **No server needed**: Just static files
- ✅ **HTTPS included**: CloudFront provides SSL

### Backend on ECS Fargate
- ✅ **Dynamic**: Runs Python application
- ✅ **Scalable**: Can scale containers up/down
- ✅ **Managed**: AWS handles infrastructure
- ✅ **Isolated**: Runs in containers

---

## What You Need to Do

### Option 1: Deploy Frontend to S3/CloudFront (Recommended)

1. **Build frontend with API URL:**
   ```bash
   cd frontend
   export VITE_API_BASE_URL=http://helix--Publi-7CO3kaz1lrRR-443122446.us-west-1.elb.amazonaws.com
   npm run build
   ```

2. **Deploy to S3:**
   ```bash
   # Use existing bucket or create new one
   aws s3 sync frontend/dist s3://helix-ai-frontend-794270057041-us-west-1 \
     --region us-west-1 \
     --delete \
     --acl public-read
   ```

3. **Set up CloudFront** (optional, for HTTPS and CDN):
   - Create CloudFront distribution pointing to S3 bucket
   - Get CloudFront URL
   - Users access frontend via CloudFront URL

### Option 2: Run Frontend Locally (Development Only)

If you're just testing:
```bash
cd frontend
export VITE_API_BASE_URL=http://helix--Publi-7CO3kaz1lrRR-443122446.us-west-1.elb.amazonaws.com
npm run dev
```

Then open `http://localhost:5173` - it will connect to your AWS backend.

---

## Summary

| Question | Answer |
|----------|--------|
| **Does frontend run locally?** | Only in development. In production, it should be on S3/CloudFront |
| **Does backend run on AWS?** | ✅ Yes - Already deployed to ECS Fargate |
| **Can both run locally?** | ✅ Yes - For development |
| **Should both be on AWS in production?** | ✅ Yes - Frontend on S3/CloudFront, Backend on ECS |

---

## Next Steps

1. **Deploy frontend to S3** (see `docs/FRONTEND_API_URL_SETUP.md`)
2. **Set up CloudFront** (optional, for HTTPS)
3. **Update DNS** (optional, for custom domain)

The backend is already working on AWS! You just need to deploy the frontend.



