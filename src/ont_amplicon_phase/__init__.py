"""
ONT Amplicon Phase Pipeline

A professional bioinformatics pipeline for analyzing barcoded amplicon sequences 
from Oxford Nanopore Technology (ONT) data, with capabilities for quality control,
variant calling, and haplotype phasing.
"""

__version__ = "1.0.0"
__author__ = "Clinical Genomics Team"
__email__ = "support@example.com"

from .core.pipeline import AmpliconPipeline
from .core.qc_analyzer import QCAnalyzer
from .core.phasing_analyzer import PhasingAnalyzer

__all__ = ["AmpliconPipeline", "QCAnalyzer", "PhasingAnalyzer"]