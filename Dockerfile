FROM python:3.9-slim

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    cmake \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    samtools \
    tabix \
    bgzip \
    awscli \
    && rm -rf /var/lib/apt/lists/*

# Install Python packages directly
RUN pip install --no-cache-dir \
    pysam>=0.19.0 \
    pandas>=1.3.0 \
    numpy>=1.21.0 \
    whatshap>=1.4 \
    click>=8.0.0 \
    pyyaml>=6.0 \
    jsonschema>=4.0.0

# Note: Clair3 and HapCUT2 should be installed on the host system
# This container is for the pipeline only

# Set working directory
WORKDIR /app

# Copy application code
COPY src/ ./src/
COPY setup.py pyproject.toml ./
COPY README.md ./

# Install the application
RUN pip install -e .

# Create directories for data and results
RUN mkdir -p /data /results /config

# Set environment variables
ENV ONT_REFERENCE_GENOME="/data/reference/genome.fna"
ENV ONT_CLAIR3_PATH="/opt/Clair3"
ENV ONT_HAPCUT2_PATH="/opt/HapCUT2/build"
CMD ["ont-amplicon-phase", "--help"]

# Labels
LABEL maintainer="Clinical Genomics Team <support@example.com>"
LABEL version="1.0.0"
LABEL description="Professional pipeline for ONT amplicon sequencing analysis and variant phasing"