# HBAT Web Application Dockerfile
# Multi-stage build for optimized production deployment

# Stage 1: Builder
FROM python:3.12-slim as builder

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies required for building Python packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    cmake \
    git \
    libgraphviz-dev \
    graphviz \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# Create and activate virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy requirements and install Python dependencies
WORKDIR /app
COPY requirements.txt .

# Install Python dependencies
RUN pip install --upgrade pip setuptools wheel && \
    pip install -r requirements.txt && \
    pip install nicegui>=1.4.0

# Copy project files for installation
COPY pyproject.toml .
COPY hbat/ ./hbat/
COPY README.md .
COPY .git .git

# Install the package
# Set SETUPTOOLS_SCM_PRETEND_VERSION as fallback if git is not available
RUN pip install -e .

# Stage 2: Runtime
FROM python:3.12-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:$PATH" \
    HBAT_ENV=production

# Install runtime system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    graphviz \
    libgraphviz-dev \
    && rm -rf /var/lib/apt/lists/*

# Create non-root user for security
RUN useradd -m -u 1000 hbat && \
    mkdir -p /app/uploads && \
    chown -R hbat:hbat /app

# Copy virtual environment from builder
COPY --from=builder /opt/venv /opt/venv

# Copy application code
WORKDIR /app
COPY --from=builder /app/hbat ./hbat
COPY --from=builder /app/pyproject.toml .
COPY --from=builder /app/README.md .

# Set ownership
RUN chown -R hbat:hbat /app

# Switch to non-root user
USER hbat

# Expose port 8080 (NiceGUI default)
EXPOSE 8080

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=5s --retries=3 \
    CMD python -c "import urllib.request; urllib.request.urlopen('http://localhost:8080').read()" || exit 1

# Set the entrypoint to run the web application
CMD ["python", "-m", "hbat.server.app"]
