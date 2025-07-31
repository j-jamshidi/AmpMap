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
    samtools \
    tabix \
    awscli \
    && rm -rf /var/lib/apt/lists/*

# Install htslib from source for HapCUT2
RUN cd /opt && \
    wget https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 && \
    tar -xjf htslib-1.19.1.tar.bz2 && \
    cd htslib-1.19.1 && \
    ./configure --prefix=/usr/local && \
    make && make install && \
    cd .. && rm -rf htslib-1.19.1*

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

# Install Clair3 directly with conda
RUN conda create -c conda-forge -c bioconda -n clair3 python=3.9 tensorflow=2.15.0 whatshap samtools parallel xz zlib bzip2 automake curl pigz cffi make gcc -y

# Clone and build Clair3
RUN cd /opt && \
    git clone --branch v1.0.8 --depth 1 https://github.com/HKU-BAL/Clair3.git && \
    cd Clair3 && \
    /opt/conda/envs/clair3/bin/make PREFIX=/opt/conda/envs/clair3 && \
    mkdir -p models && \
    cd models && \
    wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r1041_e82_400bps_sup_v500.tar.gz && \
    tar -zxf r1041_e82_400bps_sup_v500.tar.gz && \
    rm r1041_e82_400bps_sup_v500.tar.gz

# Install pypy in clair3 environment
RUN cd /opt/conda/envs/clair3/bin && \
    wget https://downloads.python.org/pypy/pypy3.10-v7.3.19-linux64.tar.bz2 && \
    tar -jxf pypy3.10-v7.3.19-linux64.tar.bz2 && \
    ln -sf pypy3.10-v7.3.19-linux64/bin/pypy3 pypy3 && \
    rm pypy3.10-v7.3.19-linux64.tar.bz2

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

# Set working directory
WORKDIR /app

# Copy application code
COPY src/ ./src/
COPY setup.py pyproject.toml ./
COPY README.md ./

# Install the application in ONT environment
RUN /opt/conda/envs/ONT/bin/pip install -e .

# Set environment variables to use built-in tools and reference
ENV ONT_REFERENCE_GENOME=/opt/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
ENV ONT_CLAIR3_PATH=/opt/Clair3
ENV ONT_HAPCUT2_PATH=/opt/HapCUT2-1.3.4/build

# Setup conda initialization and default environment
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate ONT" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

CMD ["ont-amplicon-phase", "--help"]

# Labels
LABEL maintainer="Javad Jamshidi <javad.jamshidi@neura.edu.au>"
LABEL version="1.0.0"
LABEL description="ONT amplicon sequencing analysis pipeline with all dependencies included"