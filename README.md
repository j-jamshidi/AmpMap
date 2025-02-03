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