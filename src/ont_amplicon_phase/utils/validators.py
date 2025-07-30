"""Input validation utilities."""

import pandas as pd
import pysam
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional
import re
import logging

logger = logging.getLogger(__name__)


def validate_sample_sheet(sample_sheet_path: Path) -> Tuple[bool, List[str], Optional[pd.DataFrame]]:
    """
    Validate sample sheet format and content.
    
    Returns:
        Tuple of (is_valid, error_messages, dataframe)
    """
    errors = []
    
    if not sample_sheet_path.exists():
        return False, [f"Sample sheet not found: {sample_sheet_path}"], None
    
    try:
        # Read CSV
        df = pd.read_csv(sample_sheet_path)
        
        # Check required columns
        required_columns = ['Batch', 'Barcode', 'Episode', 'Coordinate', 'Variant1', 'Variant2', 'EpisodeWES']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            errors.append(f"Missing required columns: {missing_columns}")
        
        # Validate each row
        for idx, row in df.iterrows():
            row_errors = []
            
            # Validate Barcode format
            if not isinstance(row['Barcode'], (int, str)) or str(row['Barcode']).strip() == '':
                row_errors.append("Barcode is required")
            
            # Validate Episode format
            if not isinstance(row['Episode'], str) or row['Episode'].strip() == '':
                row_errors.append("Episode is required")
            
            # Validate Coordinate format
            coord_pattern = r'^chr[0-9XYM]+:\d+-\d+$'
            if not re.match(coord_pattern, str(row['Coordinate'])):
                row_errors.append(f"Invalid coordinate format: {row['Coordinate']}")
            
            # Validate Variant format (if provided)
            variant_pattern = r'^chr[0-9XYM]+:\d+[:\s][ATCG]+[>:][ATCG]+$'
            for var_col in ['Variant1', 'Variant2']:
                variant = str(row[var_col]).strip()
                if variant and variant.lower() not in ['', 'nan', 'none']:
                    if not re.match(variant_pattern, variant):
                        row_errors.append(f"Invalid {var_col} format: {variant}")
            
            if row_errors:
                errors.append(f"Row {idx + 2}: {'; '.join(row_errors)}")
        
        return len(errors) == 0, errors, df
        
    except Exception as e:
        return False, [f"Error reading sample sheet: {str(e)}"], None


def validate_bam_file(bam_path: Path) -> Tuple[bool, List[str]]:
    """
    Validate BAM file format and accessibility.
    
    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []
    
    if not bam_path.exists():
        return False, [f"BAM file not found: {bam_path}"]
    
    try:
        # Try to open BAM file
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            # Check if file has reads
            try:
                next(iter(bam))
            except StopIteration:
                errors.append("BAM file appears to be empty")
            
            # Check if index exists
            index_path = Path(str(bam_path) + ".bai")
            if not index_path.exists():
                errors.append("BAM index file (.bai) not found")
        
        return len(errors) == 0, errors
        
    except Exception as e:
        return False, [f"Error validating BAM file: {str(e)}"]


def validate_vcf_file(vcf_path: Path) -> Tuple[bool, List[str]]:
    """
    Validate VCF file format and accessibility.
    
    Returns:
        Tuple of (is_valid, error_messages)
    """
    errors = []
    
    if not vcf_path.exists():
        return False, [f"VCF file not found: {vcf_path}"]
    
    try:
        # Try to open VCF file
        with pysam.VariantFile(str(vcf_path)) as vcf:
            # Check if file has variants
            try:
                next(iter(vcf))
            except StopIteration:
                errors.append("VCF file appears to be empty")
        
        return len(errors) == 0, errors
        
    except Exception as e:
        return False, [f"Error validating VCF file: {str(e)}"]


def validate_coordinate_format(coordinate: str) -> bool:
    """Validate genomic coordinate format."""
    pattern = r'^chr[0-9XYM]+:\d+-\d+$'
    return bool(re.match(pattern, coordinate))


def validate_variant_format(variant: str) -> bool:
    """Validate variant format."""
    if not variant or variant.lower() in ['', 'nan', 'none']:
        return True  # Empty variants are allowed
    
    # Support both formats: chr1:123 A>T and chr1:123:A:T
    pattern1 = r'^chr[0-9XYM]+:\d+\s+[ATCG]+>[ATCG]+$'
    pattern2 = r'^chr[0-9XYM]+:\d+:[ATCG]+:[ATCG]+$'
    
    return bool(re.match(pattern1, variant) or re.match(pattern2, variant))