# ONT Amplicon Phase Pipeline - Professional Edition

A production-ready, clinically-suitable bioinformatics pipeline for analyzing barcoded amplicon sequences from Oxford Nanopore Technology (ONT) data, with capabilities for quality control, variant calling, and haplotype phasing.

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/your-org/ont-amplicon-phase.git
cd ont-amplicon-phase

# Install using pip
pip install .

# Or install in development mode
pip install -e ".[dev]"
```

### Docker Installation (Recommended for Production)

```bash
# Build the Docker image
docker build -t ont-amplicon-phase:latest .

# Or use Docker Compose
docker-compose up ont-amplicon-phase
```

### Basic Usage

```bash
# Run the pipeline (uses original hardcoded paths by default)
ont-amplicon-phase run RUN_001 \
    --input-dir /EBSDataDrive/ONT/Runs/RUN_001 \
    --output-dir /EBSDataDrive/ONT/Runs/RUN_001/result

# Skip S3 upload
ont-amplicon-phase run RUN_001 --no-upload

# Skip XML generation
ont-amplicon-phase run RUN_001 --no-xml

# Custom paths
ont-amplicon-phase run RUN_001 \
    --input-dir /path/to/data \
    --output-dir /path/to/results \
    --sample-sheet /path/to/sample_sheet.csv

# Validate sample sheet
ont-amplicon-phase validate /path/to/sample_sheet.csv

# View configuration
ont-amplicon-phase config
```

## 📋 Features

### Core Capabilities
- **Automated Quality Control**: Comprehensive QC metrics for ONT amplicon data
- **Variant Calling**: Integration with Clair3 for accurate SNP and indel detection
- **Haplotype Phasing**: WhatsHap and HapCUT2 integration for variant phasing
- **Flexible Analysis**: Supports single-variant QC or two-variant phasing analysis
- **Clinical-Grade**: Designed for clinical genomics workflows

### Professional Features
- **Configuration Management**: YAML-based configuration with environment variable support
- **Input Validation**: Comprehensive validation of sample sheets and input files
- **Error Handling**: Robust error handling and logging
- **S3 Integration**: Automated upload of results to AWS S3
- **IGV Visualization**: Automatic generation of IGV XML files for visualization
- **BaseSpace Integration**: WES data integration via BaseSpace API
- **Containerization**: Docker support for reproducible deployments
- **CI/CD Ready**: GitHub Actions integration for automated testing
- **Extensible**: Modular design for easy customization

## 🏗️ Architecture

```
src/ont_amplicon_phase/
├── core/                   # Core pipeline modules
│   ├── pipeline.py        # Main pipeline orchestrator
│   ├── qc_analyzer.py     # Quality control analysis
│   ├── phasing_analyzer.py # Phasing analysis (exact replica)
│   ├── variant_processor.py # Variant processing
│   └── xml_generator.py   # IGV XML generation & S3 upload
├── utils/                  # Utility modules
│   ├── config.py          # Configuration management
│   ├── validators.py      # Input validation
│   └── file_utils.py      # File handling utilities
├── data/                   # Template files (VCF, XML)
│   ├── dummy.vcf          # VCF template
│   ├── solo_LR.xml        # IGV template (long reads only)
│   └── solo_LR_SR.xml     # IGV template (long + short reads)
└── config/                 # Default configuration
    └── default.yaml       # Default paths and settings
```

## ⚙️ Configuration

### Default Paths (from Original Implementation)

The pipeline uses the same default paths as the original scripts:

```yaml
paths:
  reference_genome: "/EFSGaiaDataDrive/ref/ONT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
  clair3_path: "/EBSDataDrive/ONT/Clair3"
  hapcut2_path: "/EBSDataDrive/ONT/HapCUT2-1.3.4/build"
  base_data_dir: "/EBSDataDrive/ONT/Runs"
  script_dir: "/EBSDataDrive/ONT/script"
  software_dir: "/EBSDataDrive/software"
  clair3_model: "r1041_e82_400bps_sup_v500"

aws:
  s3_bucket: "nswhp-gaia-poc-pl"
  s3_prefix: "ONT"
  presign_expire_seconds: 604800

basespace:
  config_name: "POWH"
  sample_files:
    - "/EBSDataDrive/software/sample_ran.txt"
    - "/EBSDataDrive/software/sample_ran_CRE_BS.txt"
```

### Environment Variables

```bash
# Core paths (defaults from original scripts)
export ONT_REFERENCE_GENOME="/EFSGaiaDataDrive/ref/ONT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
export ONT_CLAIR3_PATH="/EBSDataDrive/ONT/Clair3"
export ONT_HAPCUT2_PATH="/EBSDataDrive/ONT/HapCUT2-1.3.4/build"
export ONT_BASE_DATA_DIR="/EBSDataDrive/ONT/Runs"
export ONT_SOFTWARE_DIR="/EBSDataDrive/software"
export ONT_CLAIR3_MODEL="r1041_e82_400bps_sup_v500"

# AWS settings (defaults from original scripts)
export ONT_S3_BUCKET="nswhp-gaia-poc-pl"
export ONT_S3_PREFIX="ONT"
export ONT_UPLOAD_S3="true"
export ONT_GENERATE_XML="true"

# BaseSpace settings
export ONT_BASESPACE_CONFIG="POWH"

# Pipeline settings
export ONT_THREADS="6"
export ONT_MIN_COVERAGE="20"
export ONT_LOG_LEVEL="INFO"
```

### Configuration File

Create a custom configuration file:

```yaml
# config/production.yaml
qc:
  min_read_length: 1000
  min_mapping_quality: 20
  min_depth: 50
  min_passing_reads: 50

variant_calling:
  clair3:
    threads: 8
    min_coverage: 30
    platform: "ont"
    enable_phasing: true
    use_whatshap: true

paths:
  reference_genome: "/data/reference/GRCh38.fna"
  clair3_path: "/opt/clair3"
  hapcut2_path: "/opt/hapcut2"
  base_data_dir: "/data/runs"
  software_dir: "/data/software"

aws:
  s3_bucket: "my-custom-bucket"
  s3_prefix: "ONT"

basespace:
  config_name: "MY_CONFIG"

output:
  upload_to_s3: true
  generate_xml: true
  cleanup_intermediate: true
```

Use with: `ont-amplicon-phase --config config/production.yaml run RUN_001`

## 📊 Sample Sheet Format

```csv
Batch,Barcode,Episode,Coordinate,Variant1,Variant2,EpisodeWES
1,01,SAMPLE_001,chr1:1000000-1002000,chr1:1001000 A>T,chr1:1001500 G>C,WES_001
1,02,SAMPLE_002,chr2:2000000-2002000,chr2:2001000 T>C,,WES_002
```

### Required Columns
- **Batch**: Run batch identifier
- **Barcode**: Sample barcode (will be formatted as barcode01, barcode02, etc.)
- **Episode**: Unique sample identifier
- **Coordinate**: Amplicon genomic coordinates (chr:start-end)
- **Variant1**: First variant (optional for QC-only analysis)
- **Variant2**: Second variant (optional, required for phasing)
- **EpisodeWES**: WES episode identifier (optional)

## 🔬 Analysis Types

### 1. QC-Only Analysis
When no variants are provided, performs comprehensive quality control:
- Read length distribution and N50
- Mapping quality assessment
- Base quality statistics
- Alignment identity metrics
- Depth coverage analysis

### 2. Single Variant QC
When one variant is provided:
- All QC metrics above
- Variant validation against Clair3 calls
- Variant matching report

### 3. Two-Variant Phasing
When two variants are provided:
- Variant validation and coverage assessment
- Spanning read analysis
- WhatsHap phasing
- HapCUT2 phasing (alternative method)
- Read-based phasing analysis
- Cis/Trans determination

## 📈 Output Files

For each processed sample in `results/barcodeXX/`:

### Core Output Files
- `{episode}.bam` - Amplicon-specific BAM file
- `{episode}_report.txt` - Comprehensive analysis report
- `{episode}.wf_snp.vcf.gz` - Clair3 variant calls
- `{episode}.xml` - IGV visualization file
- `{episode}_coordinate.bed` - Amplicon coordinates

### QC-Specific Files (Single/No Variant)
- `{episode}_QC.bam` - Quality-filtered BAM file

### Phasing-Specific Files (Two Variants)
- `{episode}_phased.bam` - WhatsHap haplotagged BAM
- `{episode}_phased.vcf` - WhatsHap phased VCF
- `{episode}_Phased.vcf` - Final phased VCF
- `whatshap.log` - WhatsHap execution log
- `HapCUT2.log` - HapCUT2 execution log
- `clean-span-hq.bam` - High-quality spanning reads (temporary)

### Cloud Integration
- **S3 Upload**: All results automatically uploaded to S3
- **Presigned URLs**: Generated for IGV visualization
- **BaseSpace Integration**: WES data integration when available

## 🐳 Docker Usage

### Production Deployment

```bash
# Build image
docker build -t ont-amplicon-phase:v1.0.0 .

# Run pipeline with default paths
docker run --rm \
  -v /EBSDataDrive:/EBSDataDrive \
  -v /EFSGaiaDataDrive:/EFSGaiaDataDrive \
  -e AWS_ACCESS_KEY_ID=your_key \
  -e AWS_SECRET_ACCESS_KEY=your_secret \
  ont-amplicon-phase:v1.0.0 \
  ont-amplicon-phase run RUN_001

# Run pipeline with custom paths
docker run --rm \
  -v /path/to/data:/data \
  -v /path/to/results:/results \
  -e ONT_REFERENCE_GENOME=/data/reference.fna \
  -e ONT_S3_BUCKET=my-bucket \
  ont-amplicon-phase:v1.0.0 \
  ont-amplicon-phase run RUN_001 \
    --input-dir /data \
    --output-dir /results
```

### Development

```bash
# Start development container
docker-compose up -d ont-amplicon-dev
docker-compose exec ont-amplicon-dev bash
```

## 🧪 Testing

```bash
# Run all tests
make test

# Run specific test categories
pytest tests/unit/ -v
pytest tests/integration/ -v

# Run with coverage
pytest --cov=ont_amplicon_phase --cov-report=html
```

## 🔧 Development

### Setup Development Environment

```bash
# Install in development mode
make install-dev

# Setup pre-commit hooks
make setup-dev

# Format code
make format

# Run linting
make lint
```

### Adding New Features

1. Create feature branch: `git checkout -b feature/new-feature`
2. Implement changes in appropriate modules
3. Add tests in `tests/` directory
4. Update documentation
5. Run tests and linting: `make test lint`
6. Submit pull request

## 📋 Quality Assurance

### Code Quality
- **Type Hints**: Full type annotation coverage
- **Linting**: flake8 and mypy integration
- **Formatting**: Black code formatting
- **Testing**: pytest with coverage reporting

### Clinical Standards
- **Input Validation**: Comprehensive validation of all inputs
- **Error Handling**: Graceful error handling with detailed logging
- **Reproducibility**: Containerized execution for consistent results
- **Traceability**: Detailed logging and audit trails

## 🚀 Deployment

### Production Checklist
- [ ] Configure reference genome path
- [ ] Set up Clair3 and HapCUT2 installations
- [ ] Configure resource limits (CPU, memory)
- [ ] Set up logging and monitoring
- [ ] Test with sample data
- [ ] Configure backup and archival

### Scaling Considerations
- **Parallel Processing**: Configure thread counts based on available resources
- **Storage**: Ensure adequate storage for intermediate files
- **Memory**: Monitor memory usage for large amplicons
- **Network**: Consider network storage for large datasets

## 📞 Support

### Documentation
- [API Documentation](docs/api.md)
- [Configuration Guide](docs/configuration.md)
- [Troubleshooting Guide](docs/troubleshooting.md)

### Getting Help
- GitHub Issues: Report bugs and feature requests
- Email: support@example.com
- Documentation: https://ont-amplicon-phase.readthedocs.io

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Oxford Nanopore Technologies for ONT sequencing technology
- Clair3 team for variant calling capabilities
- WhatsHap and HapCUT2 teams for phasing algorithms
- The bioinformatics community for open-source tools and libraries