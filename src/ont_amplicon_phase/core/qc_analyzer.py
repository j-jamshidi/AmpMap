"""Quality control analysis module."""

import pysam
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)


class QCAnalyzer:
    """Performs quality control analysis on amplicon BAM files."""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.qc_thresholds = config.get('qc', {})
    
    def analyze_amplicon(self, bam_file: Path, bed_file: Path, episode: str) -> Dict[str, Any]:
        """
        Perform comprehensive QC analysis on amplicon BAM file.
        
        Args:
            bam_file: Path to BAM file
            bed_file: Path to BED file with amplicon coordinates
            episode: Sample episode identifier
        
        Returns:
            Dictionary containing QC metrics and results
        """
        logger.info(f"Starting QC analysis for {episode}")
        
        # Read amplicon region from BED file
        region = self._read_bed_file(bed_file)
        
        # Open BAM file
        with pysam.AlignmentFile(str(bam_file), "rb") as bam:
            # Collect metrics
            metrics = self._collect_read_metrics(bam, region)
            depth_stats = self._calculate_depth_stats(bam, region)
            
            # Generate QC report
            qc_results = {
                'episode': episode,
                'amplicon_region': f"{region[0]}:{region[1]}-{region[2]}",
                'amplicon_length': region[2] - region[1],
                'total_reads': metrics['total_reads'],
                'passing_qc_reads': metrics['passing_qc_reads'],
                'qc_pass_rate': metrics['passing_qc_reads'] / metrics['total_reads'] * 100 if metrics['total_reads'] > 0 else 0,
                'depth_stats': depth_stats,
                'read_length_stats': metrics['read_length_stats'],
                'mapping_quality_stats': metrics['mapping_quality_stats'],
                'base_quality_stats': metrics['base_quality_stats'],
                'identity_stats': metrics['identity_stats'],
                'strand_distribution': metrics['strand_distribution'],
                'qc_status': self._determine_qc_status(metrics, depth_stats)
            }
        
        logger.info(f"QC analysis completed for {episode}")
        return qc_results
    
    def _read_bed_file(self, bed_file: Path) -> Tuple[str, int, int]:
        """Read amplicon coordinates from BED file."""
        with open(bed_file, 'r') as f:
            line = f.readline().strip()
            parts = line.split('\t')
            return parts[0], int(parts[1]), int(parts[2])
    
    def _collect_read_metrics(self, bam: pysam.AlignmentFile, region: Tuple[str, int, int]) -> Dict[str, Any]:
        """Collect comprehensive read metrics."""
        read_lengths = []
        mapping_quals = []
        base_qualities = []
        read_identities = []
        strand_counts = defaultdict(int)
        passing_reads = 0
        
        for read in bam.fetch(region[0], region[1], region[2]):
            # Basic metrics
            read_lengths.append(read.query_length)
            mapping_quals.append(read.mapping_quality)
            
            # Strand information
            strand_counts['forward' if not read.is_reverse else 'reverse'] += 1
            
            # Base qualities
            if read.query_qualities is not None:
                mean_qual = np.mean(read.query_qualities)
                base_qualities.append(mean_qual)
            
            # Calculate read identity
            identity = self._calculate_read_identity(read)
            read_identities.append(identity)
            
            # Check if read passes QC
            if self._passes_qc_criteria(read, identity, mean_qual if read.query_qualities else 0):
                passing_reads += 1
        
        total_reads = len(read_lengths)
        
        return {
            'total_reads': total_reads,
            'passing_qc_reads': passing_reads,
            'read_length_stats': self._calculate_stats(read_lengths),
            'mapping_quality_stats': self._calculate_mapping_quality_stats(mapping_quals),
            'base_quality_stats': self._calculate_stats(base_qualities) if base_qualities else {},
            'identity_stats': self._calculate_identity_stats(read_identities),
            'strand_distribution': dict(strand_counts)
        }
    
    def _calculate_read_identity(self, read: pysam.AlignedSegment) -> float:
        """Calculate read identity percentage."""
        aligned_pairs = read.get_aligned_pairs(with_seq=True)
        if not aligned_pairs:
            return 0.0
        
        matches = sum(1 for _, _, base in aligned_pairs if base and base.isupper())
        return matches / len(aligned_pairs) if aligned_pairs else 0.0
    
    def _passes_qc_criteria(self, read: pysam.AlignedSegment, identity: float, mean_qual: float) -> bool:
        """Check if read passes QC criteria."""
        return (
            read.query_length >= self.qc_thresholds.get('min_read_length', 1000) and
            read.mapping_quality >= self.qc_thresholds.get('min_mapping_quality', 20) and
            identity >= self.qc_thresholds.get('min_identity', 0.8) and
            mean_qual >= self.qc_thresholds.get('min_base_quality', 7)
        )
    
    def _calculate_depth_stats(self, bam: pysam.AlignmentFile, region: Tuple[str, int, int]) -> Dict[str, float]:
        """Calculate depth statistics for the amplicon region."""
        depths = []
        for pileup_column in bam.pileup(region[0], region[1], region[2], truncate=True, min_base_quality=0):
            depths.append(pileup_column.n)
        
        if depths:
            return {
                'median': float(np.median(depths)),
                'mean': float(np.mean(depths)),
                'min': float(min(depths)),
                'max': float(max(depths))
            }
        return {'median': 0, 'mean': 0, 'min': 0, 'max': 0}
    
    def _calculate_stats(self, values: List[float]) -> Dict[str, float]:
        """Calculate basic statistics for a list of values."""
        if not values:
            return {}
        
        return {
            'mean': float(np.mean(values)),
            'median': float(np.median(values)),
            'min': float(min(values)),
            'max': float(max(values)),
            'n50': float(self._calculate_n50(values))
        }
    
    def _calculate_mapping_quality_stats(self, mapping_quals: List[int]) -> Dict[str, Any]:
        """Calculate mapping quality specific statistics."""
        if not mapping_quals:
            return {}
        
        stats = self._calculate_stats(mapping_quals)
        stats['above_q20_count'] = sum(1 for q in mapping_quals if q >= 20)
        stats['above_q20_percent'] = stats['above_q20_count'] / len(mapping_quals) * 100
        
        return stats
    
    def _calculate_identity_stats(self, identities: List[float]) -> Dict[str, Any]:
        """Calculate identity specific statistics."""
        if not identities:
            return {}
        
        stats = self._calculate_stats([i * 100 for i in identities])  # Convert to percentage
        stats['above_80_count'] = sum(1 for i in identities if i >= 0.8)
        stats['above_80_percent'] = stats['above_80_count'] / len(identities) * 100
        
        return stats
    
    def _calculate_n50(self, lengths: List[float]) -> float:
        """Calculate N50 for a list of lengths."""
        if not lengths:
            return 0
        
        sorted_lengths = sorted(lengths, reverse=True)
        total_length = sum(sorted_lengths)
        running_sum = 0
        
        for length in sorted_lengths:
            running_sum += length
            if running_sum >= total_length / 2:
                return length
        
        return 0
    
    def _determine_qc_status(self, metrics: Dict[str, Any], depth_stats: Dict[str, float]) -> str:
        """Determine overall QC status."""
        min_passing_reads = self.qc_thresholds.get('min_passing_reads', 50)
        min_depth = self.qc_thresholds.get('min_depth', 50)
        
        passing_reads = metrics['passing_qc_reads']
        min_depth_value = depth_stats.get('min', 0)
        
        if passing_reads >= min_passing_reads and min_depth_value >= min_depth:
            return "PASS"
        elif passing_reads >= min_passing_reads:
            return "FAIL - Insufficient depth"
        elif min_depth_value >= min_depth:
            return "FAIL - Insufficient high-quality reads"
        else:
            return "FAIL - Insufficient reads and depth"
    
    def write_qc_report(self, qc_results: Dict[str, Any], output_file: Path) -> None:
        """Write QC results to a report file."""
        with open(output_file, 'w') as f:
            f.write(f"Report for {qc_results['episode']}\n")
            f.write("-" * 30 + "\n")
            
            # Summary information
            f.write(f"Amplicon region: {qc_results['amplicon_region']}\n")
            f.write(f"Amplicon length: {qc_results['amplicon_length']:,} bp\n")
            f.write(f"Total reads: {qc_results['total_reads']}\n")
            f.write(f"Passing QC reads: {qc_results['passing_qc_reads']} ({qc_results['qc_pass_rate']:.2f}%)\n")
            
            # QC status
            f.write(f"\n*QC {qc_results['qc_status']}*\n\n")
            
            # Detailed metrics
            f.write("==Quality control details==\n")
            
            # Depth statistics
            depth = qc_results['depth_stats']
            f.write(f"Depth Statistics:\n")
            f.write(f"Median depth: {depth['median']:.0f}X\n")
            f.write(f"Minimum depth: {depth['min']:.0f}X\n")
            f.write(f"Maximum depth: {depth['max']:.0f}X\n\n")
            
            # Read length statistics
            length_stats = qc_results['read_length_stats']
            f.write(f"Read Length Distribution:\n")
            f.write(f"Median length: {length_stats['median']:.0f}\n")
            f.write(f"Mean length: {length_stats['mean']:.0f}\n")
            f.write(f"Length N50: {length_stats['n50']:.0f}\n\n")
            
            # Mapping quality
            mapq_stats = qc_results['mapping_quality_stats']
            f.write(f"Mapping Quality:\n")
            f.write(f"Mean MAPQ: {mapq_stats['mean']:.0f}\n")
            f.write(f"Reads with MAPQ≥20: {mapq_stats['above_q20_count']} ({mapq_stats['above_q20_percent']:.0f}%)\n\n")
            
            # Base quality
            if qc_results['base_quality_stats']:
                base_stats = qc_results['base_quality_stats']
                f.write(f"Base Quality:\n")
                f.write(f"Mean base quality: {base_stats['mean']:.0f}\n")
                f.write(f"Median base quality: {base_stats['median']:.0f}\n\n")
            
            # Identity statistics
            identity_stats = qc_results['identity_stats']
            f.write(f"Alignment Statistics:\n")
            f.write(f"Mean identity: {identity_stats['mean']:.0f}%\n")
            f.write(f"Reads ≥80% identity: {identity_stats['above_80_count']} ({identity_stats['above_80_percent']:.0f}%)\n\n")
            
            # Strand distribution
            strand_dist = qc_results['strand_distribution']
            total_reads = qc_results['total_reads']
            f.write(f"Strand Distribution:\n")
            f.write(f"Forward: {strand_dist.get('forward', 0)} ({strand_dist.get('forward', 0)/total_reads*100:.0f}%)\n")
            f.write(f"Reverse: {strand_dist.get('reverse', 0)} ({strand_dist.get('reverse', 0)/total_reads*100:.0f}%)\n")
        
        logger.info(f"QC report written to {output_file}")
    
    def create_filtered_bam(self, input_bam: Path, output_bam: Path, region: Tuple[str, int, int]) -> int:
        """Create a BAM file with only reads passing QC criteria."""
        passing_count = 0
        
        with pysam.AlignmentFile(str(input_bam), "rb") as infile:
            with pysam.AlignmentFile(str(output_bam), "wb", template=infile) as outfile:
                for read in infile.fetch(region[0], region[1], region[2]):
                    identity = self._calculate_read_identity(read)
                    mean_qual = np.mean(read.query_qualities) if read.query_qualities else 0
                    
                    if self._passes_qc_criteria(read, identity, mean_qual):
                        outfile.write(read)
                        passing_count += 1
        
        # Index the output BAM
        pysam.index(str(output_bam))
        logger.info(f"Created filtered BAM with {passing_count} reads: {output_bam}")
        
        return passing_count