# Installation Guide

This guide provides detailed installation instructions for the ONT Amplicon Phase Pipeline on Linux and macOS systems.

## System Requirements

### Minimum Requirements
- **OS**: Linux (Ubuntu 18.04+, CentOS 7+) or macOS 10.15+
- **CPU**: 4 cores
- **RAM**: 8 GB
- **Storage**: 50 GB free space
- **Python**: 3.8 or higher

### Recommended Requirements
- **CPU**: 8+ cores
- **RAM**: 16+ GB
- **Storage**: 100+ GB SSD
- **Network**: High-speed internet for initial setup

## Installation Methods

### Method 1: Docker Installation (Recommended)

Docker provides the most reliable and reproducible installation method.

#### Prerequisites
```bash
# Install Docker (Ubuntu/Debian)
sudo apt-get update
sudo apt-get install -y docker.io docker-compose

# Install Docker (CentOS/RHEL)
sudo yum install -y docker docker-compose

# Install Docker (macOS)
# Download Docker Desktop from https://www.docker.com/products/docker-desktop
```

#### Build and Run
```bash
# Clone repository
git clone https://github.com/your-org/ont-amplicon-phase.git
cd ont-amplicon-phase

# Build Docker image
docker build -t ont-amplicon-phase:latest .

# Verify installation
docker run --rm ont-amplicon-phase:latest ont-amplicon-phase --version
```

### Method 2: Conda Installation

#### Install Miniconda
```bash
# Download and install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# For macOS
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
bash Miniconda3-latest-MacOSX-x86_64.sh
```

#### Create Environment and Install
```bash
# Create conda environment
conda create -n ont-amplicon python=3.9 -y
conda activate ont-amplicon

# Install system dependencies
conda install -c bioconda samtools=1.15 tabix=1.15 -y

# Install Python package
pip install ont-amplicon-phase
```

### Method 3: System Installation

#### Install System Dependencies

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y \
    python3 python3-pip python3-dev \
    samtools tabix bgzip \
    build-essential cmake \
    zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-openssl-dev libssl-dev
```

**CentOS/RHEL:**
```bash
sudo yum install -y \
    python3 python3-pip python3-devel \
    samtools tabix \
    gcc gcc-c++ make cmake \
    zlib-devel bzip2-devel xz-devel \
    openssl-devel libcurl-devel
```

**macOS:**
```bash
# Install Homebrew if not already installed
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Install dependencies
brew install python@3.9 samtools htslib cmake
```

#### Install Python Package
```bash
pip3 install ont-amplicon-phase
```

## External Tool Installation

The pipeline requires several external bioinformatics tools:

### Clair3 Installation

```bash
# Clone Clair3
git clone https://github.com/HKU-BAL/Clair3.git
cd Clair3

# Install using conda
conda env create -f environment.yml
conda activate clair3

# Download models
mkdir models && cd models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/r1041_e82_400bps_sup_v500.tar.gz
tar -zxvf r1041_e82_400bps_sup_v500.tar.gz

# Set environment variable
export ONT_CLAIR3_PATH=/path/to/Clair3
```

### HapCUT2 Installation

```bash
# Clone HapCUT2
git clone https://github.com/vibansal/HapCUT2.git
cd HapCUT2

# Build
make
make install

# Set environment variable
export ONT_HAPCUT2_PATH=/path/to/HapCUT2/build
```

### WhatsHap Installation

```bash
# Install via pip (included in main installation)
pip install whatshap

# Or via conda
conda install -c bioconda whatshap
```

## Reference Genome Setup

### Download Reference Genome

```bash
# Create reference directory
mkdir -p /data/reference
cd /data/reference

# Download GRCh38 reference
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Decompress
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# Index with samtools
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Set environment variable
export ONT_REFERENCE_GENOME=/data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

## Configuration

### Environment Variables

Add to your shell profile (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
# ONT Amplicon Phase Pipeline Configuration
export ONT_REFERENCE_GENOME="/data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
export ONT_CLAIR3_PATH="/opt/Clair3"
export ONT_HAPCUT2_PATH="/opt/HapCUT2/build"
export ONT_THREADS=8
export ONT_LOG_LEVEL="INFO"
```

### Configuration File

Create a configuration file at `~/.ont_amplicon_phase/config.yaml`:

```yaml
paths:
  reference_genome: "/data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  clair3_path: "/opt/Clair3"
  hapcut2_path: "/opt/HapCUT2/build"

qc:
  min_read_length: 1000
  min_mapping_quality: 20
  min_depth: 50

variant_calling:
  clair3:
    threads: 8
    min_coverage: 20

logging:
  level: "INFO"
```

## Verification

### Test Installation

```bash
# Check version
ont-amplicon-phase --version

# Validate configuration
ont-amplicon-phase config

# Run help
ont-amplicon-phase --help
```

### Test with Sample Data

```bash
# Create test directory
mkdir -p test_run/bam_pass/barcode01

# Create sample sheet
cat > test_run/sample_sheet.csv << EOF
Batch,Barcode,Episode,Coordinate,Variant1,Variant2,EpisodeWES
1,01,TEST_SAMPLE,chr1:1000000-1002000,chr1:1001000 A>T,,TEST_WES
EOF

# Run pipeline (dry run)
ont-amplicon-phase run TEST_RUN \
    --input-dir test_run \
    --output-dir test_results \
    --dry-run
```

## Troubleshooting

### Common Issues

#### Permission Errors
```bash
# Fix permissions for Docker
sudo usermod -aG docker $USER
# Log out and back in

# Fix file permissions
chmod +x /path/to/ont-amplicon-phase
```

#### Missing Dependencies
```bash
# Check Python version
python3 --version

# Check samtools
samtools --version

# Check if tools are in PATH
which samtools tabix bgzip
```

#### Memory Issues
```bash
# Check available memory
free -h

# Monitor memory usage during execution
htop
```

### Getting Help

1. **Check logs**: Look at `ont_amplicon_phase.log` for detailed error messages
2. **Validate inputs**: Use `ont-amplicon-phase validate` to check sample sheets
3. **Test configuration**: Use `ont-amplicon-phase config` to verify settings
4. **GitHub Issues**: Report bugs at https://github.com/your-org/ont-amplicon-phase/issues

## Performance Optimization

### Resource Configuration

```bash
# Set optimal thread count (usually number of CPU cores)
export ONT_THREADS=$(nproc)

# Increase memory for Java-based tools
export JAVA_OPTS="-Xmx8g"
```

### Storage Optimization

```bash
# Use SSD storage for working directory
export TMPDIR=/fast/ssd/tmp

# Clean up intermediate files
ont-amplicon-phase run --clean-intermediate
```

## Uninstallation

### Remove Python Package
```bash
pip uninstall ont-amplicon-phase
```

### Remove Docker Images
```bash
docker rmi ont-amplicon-phase:latest
```

### Remove Conda Environment
```bash
conda env remove -n ont-amplicon
```

### Clean Up Files
```bash
# Remove configuration
rm -rf ~/.ont_amplicon_phase

# Remove downloaded references (optional)
rm -rf /data/reference
```