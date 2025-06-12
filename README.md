## Amplicon Analysis Pipeline for Oxford Nanopore Sequencing

A bioinformatics pipeline for analysing barcoded amplicon sequences from Oxford Nanopore Technology (ONT). The pipeline performs quality control, variant calling, and haplotype phasing.

### What it does
* Processes ONT amplicon sequencing data
* Performs quality control on raw sequences
* Calls variants using Clair3
* Phases specified variants from sample sheet
* Uses WhatsHap and HapCUT2 for phasing
* Generates QC metrics and phasing reports
* Detects and reports chimeric reads

# Amplicon Analysis Pipeline for Oxford Nanopore Sequencing

[cite_start]This bioinformatics pipeline analyzes barcoded amplicon sequences from Oxford Nanopore Technology (ONT)[cite: 1]. [cite_start]It processes raw sequencing data to perform quality control, variant calling, and, if two variants are provided, haplotype phasing[cite: 1]. The pipeline is designed to work with a `sample_sheet.csv` and automates the analysis for each sample specified.

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

The expected directory structure is:
/EBSDataDrive/ONT/Runs/
└── [RUNID]/
├── sample_sheet.csv
└── bam_pass/
├── barcode01/
│   └── *.bam
├── barcode02/
│   └── *.bam
└── ...


### Sample Sheet Format
The `sample_sheet.csv` must contain a header and the following columns:

| Column | Description | Example |
|---|---|---|
| `Batch` | The batch number for the run. | `7` |
| `Barcode` | The barcode identifier. Can be with or without leading zeros. | `barcode01` or `1` |
| `Episode` | A unique identifier for the sample. | `NA24385_FY2` |
| `Coordinate`| The genomic region of the amplicon. | `chr1:236010990-236033398` |
| `Variant1` | The first variant to analyze. Format: `chr:pos REF>ALT`. | `chr1:236011853 T>C` |
| `Variant2` | The second variant for phasing. Leave blank if only doing QC. | `chr1:236033263 A>G` |
| `EpisodeWES`| Corresponding episode identifier in WES data. | `NA24385` |

**Example `sample_sheet.csv`:**
```csv
Batch,Barcode,Episode,Coordinate,Variant1,Variant2,EpisodeWES
7,barcode01,NA24385_FY2,chr1:236010990-236033398,chr1:236011853 T>C,chr1:236033263 A>G
5,3,NA24385_FAE,chr3:49412078-49426468,chr3:49412205 G>A,chr3:49425854 G>A,NA24385
How to Run
Execute the main pipeline script by providing a RUNID as a command-line argument. This RUNID must correspond to a directory containing your sample_sheet.csv and bam_pass subfolder.

Bash

# Usage: ./run.sh [RUNID]
./run.sh MyRun_01
Pipeline Workflow
The run.sh script automates the following workflow for each sample in the sample sheet:

Input Preparation: The sample_sheet.csv is parsed. Barcode numbers are formatted (e.g., 1 becomes barcode01).

Conditional VCF Creation:

Two Variants: A VCF file containing both variants is created for phasing analysis.
One Variant: A VCF file with the single variant is created for QC and variant validation.
No Variants: No VCF is created; the pipeline proceeds with general QC only.
BAM File Processing: All .bam files for a given barcode are merged, sorted, and indexed. Reads corresponding to the target amplicon coordinate are extracted into a new sample-specific BAM file.

BED File Creation: A BED file is generated from the coordinates in the sample sheet to define the target region for variant calling.

Variant Calling: Clair3 is run on the sample BAM to call variants within the target amplicon region.

Phasing Analysis (if two variants provided):

HapCUT2 is run to phase the two specified variants.
The basecalling_phasing_amplicon.py script is executed to perform a detailed phasing analysis using the results from WhatsHap (run via Clair3) and HapCUT2.
QC Analysis (if one or no variants provided):

The basecalling_QC_amplicon.py script is executed to generate a detailed quality control report.
Output Generation & Cleanup: The script generates final reports, and intermediate files are removed. The results are then uploaded to an AWS S3 bucket.

Final Report: An XML file summarizing the run is generated.

Output Files
For each sample, the pipeline generates a results directory (/EBSDataDrive/ONT/Runs/[RUNID]/result/[Barcode]/) containing several key files:

[Episode]_report.txt: The primary output file.
For phasing analysis, this report includes:
Validation of the user-provided variants against Clair3's calls.
Allele balance checks to ensure the variants appear heterozygous.
QC statistics on reads spanning both variants.
A final phasing call (Cis or Trans) based on read counting, WhatsHap, and HapCUT2.
An estimate of chimeric read percentage.
For QC analysis, this report includes:
A QC Pass/Fail summary based on read count and depth.
Variant validation results (if a variant was provided).
Detailed statistics on read depth, length (N50), mapping quality, base quality, and strand bias.
[Episode].bam: The merged, sorted, and indexed BAM file containing reads for the target amplicon.
[Episode].wf_snp.vcf.gz: The gzipped VCF file of variants called by Clair3.
[Episode]_phased.bam: (Phasing only) A BAM file where reads have been tagged by WhatsHap with their haplotype (HP tag).
Core Functionality Breakdown
Quality Control (QC)
The QC script (basecalling_QC_amplicon.py) provides a high-level summary and detailed metrics. A sample receives a QC PASSED status if it meets the following criteria:

High-Quality Reads > 50: More than 50 reads pass the internal quality filter.
Minimum Depth > 50x: The sequencing depth never drops below 50x across the entire amplicon.
Variant Calling
Variant calling is performed by Clair3. The pipeline is configured to use a specific model (r1041_e82_400bps_sup_v500) optimized for ONT R10.4.1 flow cells. Clair3's internal phasing capabilities (via WhatsHap) are also enabled.

Haplotype Phasing
When two variants are provided, the pipeline uses a multi-faceted approach to determine their phase:

Variant Validation: Before phasing, the script (basecalling_phasing_amplicon.py) verifies that both variants are present in the reads at approximately heterozygous frequencies (i.e., not skewed heavily toward ref or alt).
Read Counting: The script directly counts the number of high-quality reads that support each haplotype:
Cis: Reads containing (ref-ref) or (alt-alt) alleles for the two variants.
Trans: Reads containing (ref-alt) or (alt-ref) alleles.
Tool-Based Phasing: The phase calls from WhatsHap and HapCUT2 are parsed from their VCF outputs and included in the final report, serving as a confirmation of the read-counting result.
Configuration
Key paths for tools and reference files are hardcoded at the beginning of the run.sh script. These may need to be modified depending on the system environment.

Bash

# Reference and tool paths in run.sh
REFERENCE="/EFSGaiaDataDrive/ref/ONT/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
CLAIR3_PATH="/EBSDataDrive/ONT/Clair3"
HAPCUT2_PATH="/EBSDataDrive/ONT/HapCUT2-1.3.4/build"
SCRIPT_PATH="/EBSDataDrive/ONT/script"
