# HBAT Web Application Dockerfile
# Multi-stage build for optimized production deployment

# Stage 1: Builder
FROM python:3.13-slim as builder

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies required for building Python packages and Graphviz
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    g++ \
    cmake \
    git \
    pkg-config \
    wget \
    ca-certificates \
    libexpat1-dev \
    libgd-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libltdl-dev \
    flex \
    bison \
    && rm -rf /var/lib/apt/lists/*

# Build and install latest Graphviz from source using CMake
ARG GRAPHVIZ_VERSION=13.1.1
RUN wget -O graphviz.tar.gz "https://gitlab.com/graphviz/graphviz/-/archive/${GRAPHVIZ_VERSION}/graphviz-${GRAPHVIZ_VERSION}.tar.gz" && \
    tar -xzf graphviz.tar.gz && \
    cd graphviz-${GRAPHVIZ_VERSION} && \
    mkdir build && cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local/ && \
    make -j$(nproc) && \
    make install && \
    cd ../.. && \
    rm -rf graphviz-${GRAPHVIZ_VERSION} graphviz.tar.gz && \
    ldconfig

# Create and activate virtual environment
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"

# Copy requirements and install Python dependencies
WORKDIR /app
COPY requirements.txt .

# Install Python dependencies
RUN pip install --upgrade pip setuptools wheel && \
    pip install -r requirements.txt

# Copy project files for installation
COPY pyproject.toml .
COPY hbat/ ./hbat/
COPY README.md .
COPY .git .git

# Install the package
# Set SETUPTOOLS_SCM_PRETEND_VERSION as fallback if git is not available
RUN pip install -e .

# Stage 2: Runtime
FROM python:3.13-slim

# Set environment variables
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PATH="/opt/venv/bin:$PATH" \
    HBAT_ENV=production

# Install runtime dependencies for Graphviz
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libexpat1 \
        libgd3 \
        libcairo2 \
        libpango-1.0-0 \
        libpangocairo-1.0-0 \
        libltdl7 \
    && rm -rf /var/lib/apt/lists/*

# Copy Graphviz installation from builder
COPY --from=builder /usr/local/bin/dot /usr/local/bin/
COPY --from=builder /usr/local/bin/neato /usr/local/bin/
COPY --from=builder /usr/local/bin/fdp /usr/local/bin/
COPY --from=builder /usr/local/bin/sfdp /usr/local/bin/
COPY --from=builder /usr/local/bin/circo /usr/local/bin/
COPY --from=builder /usr/local/bin/twopi /usr/local/bin/
COPY --from=builder /usr/local/lib/libcdt* /usr/local/lib/
COPY --from=builder /usr/local/lib/libcgraph* /usr/local/lib/
COPY --from=builder /usr/local/lib/libgvc* /usr/local/lib/
COPY --from=builder /usr/local/lib/libpathplan* /usr/local/lib/
COPY --from=builder /usr/local/lib/libxdot* /usr/local/lib/
COPY --from=builder /usr/local/lib/graphviz /usr/local/lib/graphviz
RUN ldconfig

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
