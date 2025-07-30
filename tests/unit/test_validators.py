"""Unit tests for validators module."""

import pytest
import pandas as pd
from pathlib import Path
import tempfile

from ont_amplicon_phase.utils.validators import (
    validate_sample_sheet,
    validate_coordinate_format,
    validate_variant_format
)


class TestValidators:
    """Test validation functions."""
    
    def test_validate_coordinate_format_valid(self):
        """Test valid coordinate formats."""
        valid_coords = [
            "chr1:1000-2000",
            "chr22:123456-789012",
            "chrX:1-1000000",
            "chrY:500-1500",
            "chrM:1-16569"
        ]
        
        for coord in valid_coords:
            assert validate_coordinate_format(coord), f"Should be valid: {coord}"
    
    def test_validate_coordinate_format_invalid(self):
        """Test invalid coordinate formats."""
        invalid_coords = [
            "1:1000-2000",  # Missing chr
            "chr1:1000",    # Missing end position
            "chr1-1000-2000",  # Wrong separator
            "chr1:abc-def",  # Non-numeric positions
            ""  # Empty string
        ]
        
        for coord in invalid_coords:
            assert not validate_coordinate_format(coord), f"Should be invalid: {coord}"
    
    def test_validate_variant_format_valid(self):
        """Test valid variant formats."""
        valid_variants = [
            "chr1:123 A>T",
            "chr1:123:A:T",
            "chr22:456789 G>C",
            "chrX:1000:T:G",
            "",  # Empty is allowed
            "nan"  # NaN is allowed
        ]
        
        for variant in valid_variants:
            assert validate_variant_format(variant), f"Should be valid: {variant}"
    
    def test_validate_variant_format_invalid(self):
        """Test invalid variant formats."""
        invalid_variants = [
            "1:123 A>T",  # Missing chr
            "chr1:abc A>T",  # Non-numeric position
            "chr1:123 A>",  # Missing alt allele
            "chr1:123 >T",  # Missing ref allele
        ]
        
        for variant in invalid_variants:
            assert not validate_variant_format(variant), f"Should be invalid: {variant}"
    
    def test_validate_sample_sheet_valid(self):
        """Test valid sample sheet."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("Batch,Barcode,Episode,Coordinate,Variant1,Variant2,EpisodeWES\n")
            f.write("1,01,TEST_SAMPLE,chr1:1000-2000,chr1:1500 A>T,chr1:1800 G>C,TEST_WES\n")
            f.flush()
            
            is_valid, errors, df = validate_sample_sheet(Path(f.name))
            
            assert is_valid
            assert len(errors) == 0
            assert df is not None
            assert len(df) == 1
        
        Path(f.name).unlink()  # Cleanup
    
    def test_validate_sample_sheet_missing_columns(self):
        """Test sample sheet with missing columns."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("Batch,Barcode,Episode\n")
            f.write("1,01,TEST_SAMPLE\n")
            f.flush()
            
            is_valid, errors, df = validate_sample_sheet(Path(f.name))
            
            assert not is_valid
            assert "Missing required columns" in errors[0]
        
        Path(f.name).unlink()  # Cleanup
    
    def test_validate_sample_sheet_invalid_data(self):
        """Test sample sheet with invalid data."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write("Batch,Barcode,Episode,Coordinate,Variant1,Variant2,EpisodeWES\n")
            f.write("1,,TEST_SAMPLE,invalid_coord,invalid_variant,,TEST_WES\n")
            f.flush()
            
            is_valid, errors, df = validate_sample_sheet(Path(f.name))
            
            assert not is_valid
            assert len(errors) > 0
        
        Path(f.name).unlink()  # Cleanup