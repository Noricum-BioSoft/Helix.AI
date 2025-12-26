# syntax=docker/dockerfile:1.7

FROM python:3.11-slim AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# System deps (build essentials for some libs)
# Install alignment tools (muscle5, clustalw) for sequence alignment
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential curl ca-certificates wget unzip \
    clustalw && \
    # Download and install muscle5.3 (modern version with -align flag support)
    # Determine architecture and download appropriate binary
    ARCH=$(dpkg --print-architecture) && \
    if [ "$ARCH" = "amd64" ]; then \
        MUSCLE_URL="https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-linux-x86.v5.3"; \
    elif [ "$ARCH" = "arm64" ] || [ "$ARCH" = "aarch64" ]; then \
        MUSCLE_URL="https://github.com/rcedgar/muscle/releases/download/v5.3/muscle-aarch64.v5.3"; \
    else \
        echo "Unsupported architecture: $ARCH" && exit 1; \
    fi && \
    echo "Downloading Muscle from $MUSCLE_URL for architecture $ARCH" && \
    wget "$MUSCLE_URL" -O /usr/local/bin/muscle && \
    chmod +x /usr/local/bin/muscle && \
    rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY backend/requirements.txt /app/backend/requirements.txt
RUN pip install --upgrade pip && pip install -r /app/backend/requirements.txt

# Copy backend code to /app/backend to maintain import structure
COPY backend /app/backend

# Copy agents directory for tool-generator-agent.md
COPY agents /app/agents

# Copy tools directory for tool imports (read_merging, etc.)
COPY tools /app/tools

# Set PYTHONPATH so imports work correctly
# Include /app/tools so tools can be imported directly (e.g., "import read_merging")
ENV PYTHONPATH=/app:/app/tools

# Default port
EXPOSE 8001

# Healthcheck (optional)
HEALTHCHECK --interval=30s --timeout=3s --start-period=20s --retries=3 \
  CMD curl -f http://localhost:8001/health || exit 1

# Run the app (from backend directory)
# Enable access logs and application logs
# --log-level info: Show INFO level logs (includes access logs and application logs)
# --access-log: Enable access log (HTTP requests)
CMD ["uvicorn", "backend.main_with_mcp:app", "--host", "0.0.0.0", "--port", "8001", "--log-level", "info", "--access-log"]


