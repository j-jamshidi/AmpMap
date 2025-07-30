"""Variant processing and comparison utilities."""

import pysam
from pathlib import Path
from typing import List, Tuple, Dict, Any, Optional
import re
import logging

logger = logging.getLogger(__name__)


class VariantProcessor:
    """Handles variant processing, comparison, and VCF operations."""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
    
    def create_variant_vcf(self, template_vcf: Path, output_vcf: Path, 
                          episode: str, variant1: str, variant2: str = "") -> None:
        """
        Create a VCF file from template with specified variants.
        
        Args:
            template_vcf: Path to template VCF file
            output_vcf: Path for output VCF file
            episode: Sample episode identifier
            variant1: First variant specification
            variant2: Second variant specification (optional)
        """
        logger.info(f"Creating VCF for {episode} with variants: {variant1}, {variant2}")
        
        # Read template
        with open(template_vcf, 'r') as f:
            content = f.read()
        
        # Replace sample name
        content = content.replace('sample', episode)
        
        # Process first variant
        if variant1 and self._is_valid_variant(variant1):
            var1_parts = self._parse_variant(variant1)
            content = content.replace('chrA', var1_parts['chrom'])
            content = content.replace('POS1', str(var1_parts['pos']))
            content = content.replace('REF1', var1_parts['ref'])
            content = content.replace('ALT1', var1_parts['alt'])
        
        # Process second variant or remove line
        if variant2 and self._is_valid_variant(variant2):
            var2_parts = self._parse_variant(variant2)
            content = content.replace('chrB', var2_parts['chrom'])
            content = content.replace('POS2', str(var2_parts['pos']))
            content = content.replace('REF2', var2_parts['ref'])
            content = content.replace('ALT2', var2_parts['alt'])
        else:
            # Remove second variant line
            lines = content.split('\n')
            content = '\n'.join([line for line in lines if 'chrB' not in line])
        
        # Write output VCF
        with open(output_vcf, 'w') as f:
            f.write(content)
        
        # Sort VCF if multiple variants
        if variant2 and self._is_valid_variant(variant2):
            self._sort_vcf(output_vcf)
        
        logger.info(f"Created VCF file: {output_vcf}")
    
    def _is_valid_variant(self, variant: str) -> bool:
        """Check if variant string is valid and not empty."""
        if not variant or variant.lower().strip() in ['', 'nan', 'none']:
            return False
        return 'chr' in variant.lower()
    
    def _parse_variant(self, variant: str) -> Dict[str, str]:
        """
        Parse variant string into components.
        
        Supports formats:
        - chr1:123 A>T
        - chr1:123:A:T
        """
        variant = variant.strip()
        
        # Format: chr1:123 A>T
        if ' ' in variant and '>' in variant:
            pos_part, allele_part = variant.split(' ', 1)
            chrom, pos = pos_part.split(':')
            ref, alt = allele_part.split('>')
        # Format: chr1:123:A:T
        elif variant.count(':') >= 3:
            parts = variant.split(':')
            chrom = parts[0]
            pos = parts[1]
            ref = parts[2]
            alt = parts[3]
        else:
            raise ValueError(f"Invalid variant format: {variant}")
        
        return {
            'chrom': chrom.strip(),
            'pos': pos.strip(),
            'ref': ref.strip().upper(),
            'alt': alt.strip().upper()
        }
    
    def _sort_vcf(self, vcf_path: Path) -> None:
        """Sort VCF file by chromosome and position."""
        with open(vcf_path, 'r') as f:
            lines = f.readlines()
        
        header_lines = [line for line in lines if line.startswith('#')]
        variant_lines = [line for line in lines if not line.startswith('#')]
        
        # Sort variant lines
        def sort_key(line):
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom = parts[0]
                pos = int(parts[1])
                # Extract numeric part from chromosome for sorting
                chrom_num = re.findall(r'\d+', chrom)
                chrom_num = int(chrom_num[0]) if chrom_num else 999
                return (chrom_num, pos)
            return (999, 0)
        
        variant_lines.sort(key=sort_key)
        
        # Write sorted VCF
        with open(vcf_path, 'w') as f:
            f.writelines(header_lines)
            f.writelines(variant_lines)
    
    def compare_variants(self, user_vcf: Path, called_vcf: Path) -> List[Dict[str, Any]]:
        """
        Compare user-provided variants with called variants.
        
        Args:
            user_vcf: Path to user-provided VCF
            called_vcf: Path to variant calling VCF
        
        Returns:
            List of comparison results
        """
        logger.info(f"Comparing variants: {user_vcf} vs {called_vcf}")
        
        # Read user variants
        user_variants = self._read_variants_from_vcf(user_vcf)
        
        # Read called variants
        called_variants = self._read_variants_from_called_vcf(called_vcf)
        
        # Compare variants
        comparison_results = []
        for i, user_var in enumerate(user_variants, 1):
            found = False
            for called_var in called_variants:
                if self._variants_match(user_var, called_var):
                    found = True
                    comparison_results.append({
                        'variant_num': i,
                        'variant': user_var,
                        'found': True,
                        'matched_variant': called_var
                    })
                    break
            
            if not found:
                comparison_results.append({
                    'variant_num': i,
                    'variant': user_var,
                    'found': False,
                    'matched_variant': None
                })
        
        return comparison_results
    
    def _read_variants_from_vcf(self, vcf_path: Path) -> List[Dict[str, str]]:
        """Read variants from a standard VCF file."""
        variants = []
        
        with open(vcf_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        variants.append({
                            'chrom': parts[0],
                            'pos': int(parts[1]),
                            'ref': parts[3],
                            'alt': parts[4]
                        })
        
        return variants
    
    def _read_variants_from_called_vcf(self, vcf_path: Path) -> List[Dict[str, Any]]:
        """Read variants from a variant calling VCF file."""
        variants = []
        
        try:
            with pysam.VariantFile(str(vcf_path)) as vcf:
                for record in vcf.fetch():
                    if record.alts and record.alts[0] and record.alts[0] != ".":
                        variant_info = {
                            'chrom': record.chrom,
                            'pos': record.pos,
                            'ref': record.ref,
                            'alt': record.alts[0],
                            'qual': record.qual,
                            'gt': None
                        }
                        
                        # Extract genotype if available
                        if record.samples:
                            sample = list(record.samples.values())[0]
                            if 'GT' in sample:
                                gt = sample['GT']
                                variant_info['gt'] = "/".join(str(x) for x in gt) if gt else "."
                        
                        variants.append(variant_info)
        except Exception as e:
            logger.error(f"Error reading called variants: {e}")
        
        return variants
    
    def _variants_match(self, user_var: Dict[str, Any], called_var: Dict[str, Any]) -> bool:
        """Check if two variants match."""
        return (
            user_var['chrom'] == called_var['chrom'] and
            user_var['pos'] == called_var['pos'] and
            user_var['ref'] == called_var['ref'] and
            user_var['alt'] == called_var['alt']
        )
    
    def write_variant_comparison_report(self, comparison_results: List[Dict[str, Any]], 
                                      called_variants: List[Dict[str, Any]], 
                                      output_file: Path, user_variants: List[Dict[str, str]]) -> None:
        """Write variant comparison results to report file."""
        with open(output_file, 'a') as f:
            f.write("=== Variant Calling ===\n")
            
            # Write user-provided variants
            f.write("User provided variant(s):\n")
            for var in user_variants:
                f.write(f"{var['chrom']}\t{var['pos']}\t{var['ref']}\t{var['alt']}\n")
            
            # Calculate distance if two variants
            if len(user_variants) == 2:
                distance = abs(user_variants[1]['pos'] - user_variants[0]['pos'])
                f.write(f"Distance between the variants is {distance:,} bp\n")
            
            f.write("\nVariants called from the amplicon by Clair3:\n")
            f.write("CHROM\tPOS\tREF\tALT\tQUAL\tGT\n")
            
            # Create set of user variants for matching
            user_variant_set = {
                (var['chrom'], var['pos'], var['ref'], var['alt']) 
                for var in user_variants
            }
            
            # Write called variants with indicators
            for var in called_variants:
                qual_str = f"{var['qual']:.1f}" if var['qual'] is not None else "."
                gt_str = var.get('gt', '.')
                
                # Check if this variant matches user variants
                is_matched = (var['chrom'], var['pos'], var['ref'], var['alt']) in user_variant_set
                indicator = "<-" if is_matched else ""
                
                f.write(f"{var['chrom']}\t{var['pos']}\t{var['ref']}\t{var['alt']}\t{qual_str}\t{gt_str} {indicator}\n")
            
            f.write("\nVariant Matching Results:\n")
            for result in comparison_results:
                var = result['variant']
                if result['found']:
                    f.write(f"Variant {result['variant_num']}: {var['chrom']} {var['pos']} {var['ref']} > {var['alt']} - FOUND\n")
                else:
                    f.write(f"Variant {result['variant_num']}: {var['chrom']} {var['pos']} {var['ref']} > {var['alt']} - NOT FOUND\n")
        
        logger.info(f"Variant comparison report written to {output_file}")