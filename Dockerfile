FROM continuumio/miniconda3:23.5.2-0

# Prevent interactive prompts
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
    libhts-dev \
    samtools \
    tabix \
    awscli \
    docker.io \
    && rm -rf /var/lib/apt/lists/*

# Create conda environment for ONT tools
RUN conda create -n ONT python=3.9 -y
ENV PATH /opt/conda/envs/ONT/bin:$PATH

# Install Python packages in ONT environment
RUN /opt/conda/envs/ONT/bin/pip install --no-cache-dir \
    pysam>=0.19.0 \
    pandas>=1.3.0 \
    numpy>=1.21.0 \
    whatshap>=1.4 \
    click>=8.0.0 \
    pyyaml>=6.0 \
    jsonschema>=4.0.0

# Note: Clair3 Docker image will be pulled on first use

# Install HapCUT2 (pinned version)
RUN cd /opt && \
    wget https://github.com/vibansal/HapCUT2/archive/v1.3.4.tar.gz && \
    tar -zxf v1.3.4.tar.gz && \
    cd HapCUT2-1.3.4 && \
    make

# Download reference genome
RUN mkdir -p /opt/reference && \
    cd /opt/reference && \
    wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz && \
    gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz && \
    samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Create Clair3 wrapper script that pulls image if needed
RUN mkdir -p /opt/bin && \
    echo '#!/bin/bash' > /opt/bin/run_clair3.sh && \
    echo 'if ! docker image inspect hkubal/clair3:latest >/dev/null 2>&1; then' >> /opt/bin/run_clair3.sh && \
    echo '  echo "Pulling Clair3 Docker image..."' >> /opt/bin/run_clair3.sh && \
    echo '  docker pull hkubal/clair3:latest' >> /opt/bin/run_clair3.sh && \
    echo 'fi' >> /opt/bin/run_clair3.sh && \
    echo 'docker run --rm -v "$PWD":"$PWD" hkubal/clair3:latest /opt/bin/run_clair3.sh "$@"' >> /opt/bin/run_clair3.sh && \
    chmod +x /opt/bin/run_clair3.sh

# Set working directory
WORKDIR /app

# Copy application code
COPY src/ ./src/
COPY setup.py pyproject.toml ./
COPY README.md ./

# Install the application in ONT environment
RUN /opt/conda/envs/ONT/bin/pip install -e .

# Setup conda initialization and default environment
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate ONT" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

CMD ["ont-amplicon-phase", "--help"]

# Labels
LABEL maintainer="Javad Jamshidi <javad.jamshidi@neura.edu.au>"
LABEL version="1.0.0"
LABEL description="ONT amplicon sequencing analysis pipeline with all dependencies included"