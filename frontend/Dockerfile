# syntax=docker/dockerfile:1.7

FROM node:18-alpine AS builder
WORKDIR /app

COPY frontend/package.json frontend/package-lock.json* ./ 
RUN npm ci

COPY frontend /app
RUN npm run build

# NGINX stage to serve static files
FROM nginx:1.25-alpine AS runtime
WORKDIR /usr/share/nginx/html
COPY --from=builder /app/dist /usr/share/nginx/html

# Basic security headers (optional simple config)
COPY <<'EOF' /etc/nginx/conf.d/default.conf
server {
  listen 80;
  server_name _;

  root /usr/share/nginx/html;
  index index.html;

  location / {
    try_files $uri $uri/ /index.html;
  }

  # Security headers
  add_header X-Content-Type-Options nosniff;
  add_header X-Frame-Options DENY;
  add_header X-XSS-Protection "1; mode=block";
  add_header Referrer-Policy no-referrer-when-downgrade;
  add_header Permissions-Policy "geolocation=(), microphone=(), camera=()";
}
EOF

EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]


