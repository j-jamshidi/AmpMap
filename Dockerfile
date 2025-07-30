FROM ubuntu:20.04

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
    python3 \
    python3-pip \
    python3-dev \
    samtools \
    tabix \
    bgzip \
    awscli \
    && rm -rf /var/lib/apt/lists/*

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

# Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# Create conda environment
RUN conda create -n ont-amplicon python=3.9 -y && \
    conda clean -a

# Activate environment and install Python packages
SHELL ["conda", "run", "-n", "ont-amplicon", "/bin/bash", "-c"]

RUN pip install --no-cache-dir \
    pysam>=0.19.0 \
    pandas>=1.3.0 \
    numpy>=1.21.0 \
    whatshap>=1.4 \
    click>=8.0.0 \
    pyyaml>=6.0 \
    jsonschema>=4.0.0

# Install Clair3
WORKDIR /opt
RUN git clone https://github.com/HKU-BAL/Clair3.git && \
    cd Clair3 && \
    conda env create -f environment.yml && \
    echo "source activate clair3" >> ~/.bashrc

# Install HapCUT2
RUN git clone https://github.com/vibansal/HapCUT2.git && \
    cd HapCUT2 && \
    make && \
    make install

# Set working directory
WORKDIR /app

# Copy application code
COPY src/ ./src/
COPY setup.py pyproject.toml ./
COPY README.md ./

# Install the application
RUN conda run -n ont-amplicon pip install -e .

# Create directories for data and results
RUN mkdir -p /data /results /config

# Set environment variables
ENV ONT_REFERENCE_GENOME="/data/reference/genome.fna"
ENV ONT_CLAIR3_PATH="/opt/Clair3"
ENV ONT_HAPCUT2_PATH="/opt/HapCUT2/build"
ENV CONDA_DEFAULT_ENV="ont-amplicon"

# Create entrypoint script
RUN echo '#!/bin/bash\n\
source /opt/conda/etc/profile.d/conda.sh\n\
conda activate ont-amplicon\n\
exec "$@"' > /entrypoint.sh && \
    chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["ont-amplicon-phase", "--help"]

# Labels
LABEL maintainer="Clinical Genomics Team <support@example.com>"
LABEL version="1.0.0"
LABEL description="Professional pipeline for ONT amplicon sequencing analysis and variant phasing"