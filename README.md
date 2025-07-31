# ONT Amplicon Phase Pipeline

A bioinformatics pipeline for analyzing barcoded amplicon sequences from Oxford Nanopore Technology (ONT) data, with capabilities for quality control, variant calling, and haplotype phasing.

## Installation

```bash
# Clone and install
git clone https://github.com/your-org/ont-amplicon-phase.git
cd ont-amplicon-phase
pip install .
```

### Docker

See [Docker Usage](#docker-usage) section below for complete instructions.

## Usage

```bash
# Basic run
ont-amplicon-phase run RUN_001

# Custom paths
ont-amplicon-phase run RUN_001 \
    --input-dir /path/to/data \
    --output-dir /path/to/results

# Validate sample sheet
ont-amplicon-phase validate sample_sheet.csv
```

## Features

- **Quality Control**: Comprehensive QC metrics for ONT amplicon data
- **Variant Calling**: Clair3 integration for SNP and indel detection
- **Haplotype Phasing**: WhatsHap and HapCUT2 for variant phasing
- **Flexible Analysis**: Single-variant QC or two-variant phasing
- **S3 Integration**: Automated upload of results to AWS S3
- **IGV Visualization**: Automatic generation of IGV XML files
- **Docker Support**: Containerized execution

## Architecture

```
src/ont_amplicon_phase/
├── core/           # Pipeline modules
├── utils/          # Utilities
├── data/           # Templates
└── config/         # Configuration
```

## Configuration

The pipeline uses the original default paths from the bash scripts:

```yaml
paths:
  reference_genome: "/EFSGaiaDataDrive/ref/ONT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  clair3_path: "/EBSDataDrive/ONT/Clair3"
  base_data_dir: "/EBSDataDrive/ONT/Runs"
  hapcut2_path: "/EBSDataDrive/ONT/HapCUT2-1.3.4/build"
  software_dir: "/EBSDataDrive/software"

aws:
  s3_bucket: "nswhp-gaia-poc-pl"
  s3_prefix: "ONT"
```

Override with environment variables if needed:

```bash
export ONT_REFERENCE_GENOME="/path/to/reference.fna"
export ONT_CLAIR3_PATH="/path/to/clair3"
export ONT_S3_BUCKET="your-bucket"
```

## Sample Sheet Format

```csv
Batch,Barcode,Episode,Coordinate,Variant1,Variant2,EpisodeWES
1,01,SAMPLE_001,chr1:1000000-1002000,chr1:1001000 A>T,chr1:1001500 G>C,WES_001
1,02,SAMPLE_002,chr2:2000000-2002000,chr2:2001000 T>C,,WES_002
```

**Required columns**: Batch, Barcode, Episode, Coordinate  
**Optional**: Variant1, Variant2 (for phasing), EpisodeWES

## Analysis Types

1. **QC-Only**: No variants → comprehensive quality control
2. **Single Variant**: One variant → QC + variant validation  
3. **Two-Variant Phasing**: Two variants → QC + phasing analysis

## Output Files

For each sample in `results/barcodeXX/`:

- `{episode}.bam` - Amplicon BAM file
- `{episode}_report.txt` - Analysis report
- `{episode}.wf_snp.vcf.gz` - Variant calls
- `{episode}.xml` - IGV visualization file
- `{episode}_phased.vcf` - Phased variants (if applicable)

Results are automatically uploaded to S3 with presigned URLs for visualization.

## Docker Usage

### Build the Image

```bash
# Build image (includes all tools and reference genome)
docker build -t ont-amplicon-phase .
```

### Method 1: Using docker-compose (Recommended)

```bash
# Set environment variables
export DATA_DIR=/path/to/your/data
export RESULTS_DIR=/path/to/results
export AWS_ACCESS_KEY_ID=your_key
export AWS_SECRET_ACCESS_KEY=your_secret

# Run pipeline
docker-compose run ont-amplicon-phase ont-amplicon-phase run RUN_001

# Validate sample sheet
docker-compose run ont-amplicon-phase ont-amplicon-phase validate sample_sheet.csv
```

### Method 2: Direct Docker Run

```bash
# Run pipeline (reference genome is built into container)
docker run --rm \
  -v /EBSDataDrive:/EBSDataDrive \
  -e AWS_ACCESS_KEY_ID=your_key \
  -e AWS_SECRET_ACCESS_KEY=your_secret \
  ont-amplicon-phase \
  ont-amplicon-phase run RUN_001 --input-dir /EBSDataDrive/ONT/Runs/RUN_001

# Custom input/output directories
docker run --rm \
  -v /path/to/input:/input \
  -v /path/to/output:/output \
  -e AWS_ACCESS_KEY_ID=your_key \
  ont-amplicon-phase \
  ont-amplicon-phase run RUN_001 --input-dir /input --output-dir /output
```

### Docker Environment Variables

```bash
# Required for S3 upload
AWS_ACCESS_KEY_ID=your_access_key
AWS_SECRET_ACCESS_KEY=your_secret_key

# Optional overrides
ONT_S3_BUCKET=your-custom-bucket
ONT_S3_PREFIX=your-prefix
ONT_LOG_LEVEL=DEBUG
```

### Docker Volume Mapping

- **Data directory**: `/path/to/data` → `/EBSDataDrive/ONT/Runs`
- **Results directory**: `/path/to/results` → `/EBSDataDrive/ONT/Runs/results`
- **Software directory**: `/path/to/software` → `/EBSDataDrive/software`

## Requirements

- Python 3.9+
- Conda environment with ONT tools
- Clair3, WhatsHap, HapCUT2
- AWS credentials (for S3 upload)

## Development

```bash
pip install -e ".[dev]"
pytest tests/
```

## License

MIT License