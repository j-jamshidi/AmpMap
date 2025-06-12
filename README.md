# Amplicon Analysis Pipeline for Oxford Nanopore Sequencing

This bioinformatics pipeline analyzes barcoded amplicon sequences from Oxford Nanopore Technology (ONT). [cite_start]It processes raw sequencing data to perform quality control, variant calling, and, if two variants are provided, haplotype phasing[cite: 1].

## Table of Contents
- [What it does](#what-it-does)
- [Prerequisites](#prerequisites)
- [Input Files](#input-files)
  - [Sample Sheet Format](#sample-sheet-format)
- [How to Run](#how-to-run)
- [Pipeline Workflow](#pipeline-workflow)
- [Output Files](#output-files)
- [Core Functionality Breakdown](#core-functionality-breakdown)
  - [Quality Control (QC)](#quality-control-qc)
  - [Variant Calling](#variant-calling)
  - [Haplotype Phasing](#haplotype-phasing)
- [Configuration](#configuration)

---

## What it does
* [cite_start]**Processes ONT amplicon sequencing data** from FASTQ to phased haplotypes[cite: 1].
* [cite_start]**Performs detailed Quality Control (QC)** on raw sequences, assessing metrics like read depth, read length, mapping quality, and base quality[cite: 1].
* [cite_start]**Calls variants using Clair3**, a highly accurate tool for ONT data[cite: 1].
* [cite_start]**Phases specified variants** to determine if they are on the same chromosome (in *cis*) or on opposite chromosomes (in *trans*)[cite: 1].
* [cite_start]**Uses WhatsHap and HapCUT2** for robust phasing analysis[cite: 1].
* [cite_start]**Generates comprehensive QC and phasing reports** for each sample[cite: 1].
* [cite_start]**Detects and reports chimeric reads** as part of the phasing analysis[cite: 1].

---

## Prerequisites
The pipeline relies on several bioinformatics tools and specific file paths that must be configured in the `run.sh` script.

**Software/Tools:**
* `samtools`
* `Clair3`
* `WhatsHap`
* `HapCUT2`
* `conda` (with a configured `ONT` environment)
* `aws-cli` (for uploading results to S3)
* `Python 3` with the following libraries:
    * `pysam`
    * `pandas`
    * `numpy`

**Reference Files:**
* A reference genome in FASTA format (e.g., `GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`).
* A pre-trained Clair3 model.

---

## Input Files
The primary input for the pipeline is a `sample_sheet.csv` file located in a run-specific directory. The pipeline also requires pass `.bam` files for each barcode, which should be located in a `bam_pass` subdirectory.

/EBSDataDrive/ONT/Runs/
└── [RUNID]/
├── sample_sheet.csv
└── bam_pass/
├── barcode01/
│ └── *.bam
├── barcode02/
│ └── *.bam
└── ...


- **[RUNID]/**: Folder name for the specific sequencing run.
- **sample_sheet.csv**: Input sheet containing metadata, variant info, and amplicon coordinates.
- **bam_pass/**: Directory with BAM files separated by barcode.
  - **barcodeXX/**: Subfolder for each barcode containing one or more `.bam` files.


### Sample Sheet Format
The `sample_sheet.csv` must contain a header and the columns are processed by the `run.sh` script.

| Column | Description | Example |
|---|---|---|
| `Batch` | The batch number for the run. | `7` |
| `Barcode` | The barcode identifier. Can be with or without leading zeros. | `barcode01` or `1` |
| `Episode` | A unique identifier for the sample. | `NA24385_FY2` |
| `Coordinate`| The genomic region of the amplicon. | `chr1:236010990-236033398` |
| `Variant1` | The first variant to analyze. Format: `chr:pos REF>ALT`. | `chr1:236011853 T>C` |
| `Variant2` | The second variant for phasing. Leave blank if only doing QC. | `chr1:236033263 A>G` |
| `EpisodeWES`| Corresponding episode identifier in WES data. | `NA24385` |

---

## How to Run
Execute the main pipeline script by providing a `RUNID` as a command-line argument. This `RUNID` must correspond to a directory containing your `sample_sheet.csv` and `bam_pass` subfolder.

```bash
# Usage: ./run.sh [RUNID]
./run.sh MyRun_01

# Pipeline Workflow

The `run.sh` script automates the following workflow for each sample in the sample sheet:

### 1. Input Preparation
- The `sample_sheet.csv` is parsed, and barcode numbers are formatted (e.g., `1` becomes `barcode01`).

### 2. Conditional VCF Creation
- Based on the number of variants in the sample sheet, a VCF file is prepared for analysis:
  - **Two Variants**: A VCF file containing both variants is created for phasing.
  - **One or No Variants**: A VCF file is created if one variant is provided; otherwise, this step is skipped.

### 3. BAM File Processing
- All `.bam` files for a given barcode are merged, sorted, and indexed.
- Reads for the target amplicon are extracted into a new BAM file.

### 4. BED File Creation
- A BED file is generated from the coordinates in the sample sheet to define the target region for variant calling.

### 5. Variant Calling
- **Clair3** is run on the sample BAM file to call variants within the target region.

### 6. Phasing Analysis *(if two variants provided)*
- **HapCUT2** is run to phase the two specified variants.
- The `basecalling_phasing_amplicon.py` script is executed for detailed phasing analysis.

### 7. QC Analysis *(if one or no variants provided)*
- The `basecalling_QC_amplicon.py` script is executed to generate a QC report.

### 8. Output Generation & Cleanup
- The script generates final reports, removes intermediate files, and uploads results to an AWS S3 bucket.

### 9. Final Report
- An XML file summarizing the run is generated.

---

# Output Files

For each sample, the pipeline generates a results directory containing several key files:

- `[Episode]_report.txt`:  
  The primary output file generated by the Python scripts.  
  - For **phasing analysis**: includes validation of variants, allele balance checks, QC stats on spanning reads, a final phasing call (*Cis* or *Trans*), and an estimate of chimeric reads.  
  - For **QC analysis**: includes a QC Pass/Fail summary, variant validation, and detailed statistics on read depth, length (N50), mapping quality, and base quality.

- `[Episode].bam`:  
  The merged, sorted, and indexed BAM file for the target amplicon.

- `[Episode].wf_snp.vcf.gz`:  
  The gzipped VCF file of variants called by Clair3.

- `[Episode]_phased.bam`: *(Phasing only)*  
  A BAM file where reads are tagged by WhatsHap with their haplotype (`HP` tag).

---

# Core Functionality Breakdown

## Quality Control (QC)

The QC script (`basecalling_QC_amplicon.py`) provides a high-level summary and detailed metrics.  
A sample receives a **QC PASSED** status if it meets the following criteria:

- **High-Quality Reads > 50**: More than 50 reads pass the internal quality filter.
- **Minimum Depth > 50×**: Sequencing depth does not drop below 50× across the entire amplicon.

## Variant Calling

Variant calling is performed by **Clair3**.  
The pipeline is configured to use the following model:

- `r1041_e82_400bps_sup_v500`

Clair3's internal phasing functionality (via WhatsHap) is also enabled.

## Haplotype Phasing

When two variants are provided, the pipeline applies a multi-step strategy to determine phase:

- **Variant Validation**:  
  The script (`basecalling_phasing_amplicon.py`) checks that both variants are present in the reads at heterozygous frequencies.

- **Read Counting**:  
  The number of high-quality reads supporting each haplotype is counted:
  - *Cis*: Reads containing `(ref-ref)` or `(alt-alt)`
  - *Trans*: Reads containing `(ref-alt)` or `(alt-ref)`

- **Tool-Based Phasing**:  
  The phase calls from WhatsHap and HapCUT2 are parsed from their VCF outputs and reported as confirmation.

---

# Configuration

Key paths for tools and reference files are hardcoded at the beginning of the `run.sh` script and may need to be modified as required.

# Reference and Tool Paths in `run.sh`

The following key paths are hardcoded in the `run.sh` script:

- `REFERENCE`:  
  `/EFSGaiaDataDrive/ref/ONT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna`  
  Path to the reference genome used for alignment and variant calling.

- `CLAIR3_PATH`:  
  `/EBSDataDrive/ONT/Clair3`  
  Path to the Clair3 variant calling tool directory.

- `HAPCUT2_PATH`:  
  `/EBSDataDrive/ONT/HapCUT2-1.3.4/build`  
  Path to the compiled HapCUT2 binary directory.

- `SCRIPT_PATH`:  
  `/EBSDataDrive/ONT/script`  
  Path to custom scripts used in the pipeline (e.g., QC and phasing analysis).

