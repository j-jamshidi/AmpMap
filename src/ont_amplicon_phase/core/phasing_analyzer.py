"""Phasing analysis module """

import pysam
import subprocess
import os
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from collections import Counter
import logging

logger = logging.getLogger(__name__)


class PhasingAnalyzer:
    """Performs phasing analysis on amplicon data with two variants."""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.reference_fasta = config['paths']['reference_genome']
    
    def analyze_phasing(self, bam_file: Path, vcf_file: Path, episode: str, output_dir: Path) -> Dict[str, Any]:
        """Main phasing analysis function - replicates original main() function."""
        logger.info(f"Starting phasing analysis for {episode}")
        
        # Create output names
        output_bam = output_dir / f"{episode}_phased.bam"
        output_txt = output_dir / f"{episode}_report.txt"
        bed_file = output_dir / f"{episode}_coordinate.bed"
        
        # Read amplicon coordinates from BED file
        try:
            with open(bed_file, 'r') as bed:
                line = bed.readline().strip()
                chrom, start, end = line.split()[:3]
                start, end = int(start), int(end)
                amplicon_length = end - start
                region = (chrom, start, end)
        except FileNotFoundError:
            logger.error(f"BED file '{bed_file}' not found")
            raise
        
        # Capture all output to report file
        with open(output_txt, 'w') as output_file:
            # Write header
            print(f"Report for: {episode}", file=output_file)
            print(f"{'-'*30}", file=output_file)
            print(f"Amplicon region: {chrom}:{start}-{end}", file=output_file)
            print(f"Amplicon length: {amplicon_length:,} bp", file=output_file)
        
        # Look for the variant calling VCF
        variant_calling_vcf = output_dir / f"{episode}.wf_snp.vcf.gz"
        
        # Write variant comparison results
        if vcf_file.exists() and variant_calling_vcf.exists():
            self.write_variant_comparison_results(str(vcf_file), str(variant_calling_vcf), str(output_txt))
        
        # Perform quality control and get clean spanning reads
        clean_span_bam = output_dir / "clean-span-hq.bam"
        high_quality_reads = self.quality_control(str(bam_file), str(vcf_file), str(clean_span_bam), str(output_txt))
        
        results = {
            'episode': episode,
            'high_quality_reads': high_quality_reads,
            'output_files': {
                'report': output_txt,
                'phased_bam': output_bam,
                'clean_span_bam': clean_span_bam
            }
        }
        
        if high_quality_reads > 0:
            # Run WhatsHap phasing
            phased_vcf_gz = self.run_whatshap(str(clean_span_bam), str(vcf_file), str(output_bam))
            
            # Analyze reads for phasing
            self.analyze_reads(str(output_bam), str(vcf_file), str(output_txt))
            
            # Unzip and analyze WhatsHap phased VCF
            if phased_vcf_gz:
                unzip_command = f"gunzip -f {phased_vcf_gz}"
                subprocess.run(unzip_command, shell=True, check=True)
                
                phased_vcf = os.path.splitext(phased_vcf_gz)[0]
                whatshap_phasing = self.parse_vcf_phasing(phased_vcf)
                
                with open(output_txt, 'a') as f:
                    print(" "*40, file=f)
                    print(f"WhatsHap analysis determined the phase as {whatshap_phasing}", file=f)
                
                # Look for and analyze HapCUT2 phased VCF
                hapcut2_vcf = output_dir / f"hap2cut_{episode}.phased.VCF"
                if hapcut2_vcf.exists():
                    try:
                        hapcut2_phasing = self.parse_vcf_phasing(str(hapcut2_vcf))
                        with open(output_txt, 'a') as f:
                            print(f"\nHapCUT2 analysis determined the phase as {hapcut2_phasing}\n", file=f)
                    except Exception as e:
                        with open(output_txt, 'a') as f:
                            print(f"Error analyzing HapCUT2 VCF: {e}", file=f)
                
                results['whatshap_phasing'] = whatshap_phasing
                results['phased_vcf'] = phased_vcf
        else:
            with open(output_txt, 'a') as f:
                print("No high quality spanning reads found. Analysis cannot proceed.", file=f)
        
        logger.info(f"Phasing analysis completed for {episode}")
        return results
    
    def normalize_allele(self, allele: str) -> str:
        """Normalize allele to uppercase and replace empty alleles with <DEL>."""
        if allele == "":
            return "<DEL>"
        return allele.upper()
    
    def get_allele(self, read, variant) -> Optional[str]:
        """Determine which allele is present in a read at a given variant position - exact replica."""
        pos = variant.pos - 1  # Convert to 0-based position
        
        # Check if the variant position is covered by the read
        if pos < read.reference_start or pos >= read.reference_end:
            return None
        
        # Get the reference and alternate alleles
        ref_allele = variant.ref.upper()
        alt_allele = variant.alts[0].upper() if variant.alts else ""
        
        # Get the aligned pairs with sequence
        aligned_pairs = read.get_aligned_pairs(with_seq=True)
        
        # Extract the variant region with context
        context_size = max(len(ref_allele), len(alt_allele)) + 10  # Increased context size
        variant_region = []
        ref_pos_idx = None
        
        # Find the variant position and collect context
        for idx, (qpos, rpos, base) in enumerate(aligned_pairs):
            if rpos == pos:
                ref_pos_idx = idx
                # Collect preceding context
                start_idx = max(0, idx - context_size//2)
                variant_region.extend(aligned_pairs[start_idx:idx])
                break
        
        if ref_pos_idx is None:
            return None
        
        # Collect following context
        end_idx = min(len(aligned_pairs), ref_pos_idx + context_size//2)
        variant_region.extend(aligned_pairs[ref_pos_idx:end_idx])
        
        # Extract sequences for comparison
        read_seq = ""
        ref_seq = ""
        read_positions = []
        ref_positions = []
        
        for qpos, rpos, base in variant_region:
            if qpos is not None and rpos is not None:
                read_seq += read.query_sequence[qpos]
                ref_seq += base.upper() if base else 'N'
                read_positions.append(qpos)
                ref_positions.append(rpos)
            elif qpos is None and rpos is not None:
                # Deletion in read
                read_seq += '-'
                ref_seq += base.upper() if base else 'N'
                read_positions.append(None)
                ref_positions.append(rpos)
            elif qpos is not None and rpos is None:
                # Insertion in read
                read_seq += read.query_sequence[qpos]
                ref_seq += '-'
                read_positions.append(qpos)
                ref_positions.append(None)
        
        # Handle different variant types
        if len(ref_allele) != len(alt_allele):
            # Indel case
            variant_start_idx = ref_positions.index(pos) if pos in ref_positions else -1
            if variant_start_idx == -1:
                return None
                
            # For insertions 
            if len(alt_allele) > len(ref_allele):
                # Extract the sequence around the insertion point
                insertion_length = len(alt_allele) - len(ref_allele)
                sequence_window = 5  # Consider surrounding context
                
                # Get the read sequence at and after the variant position
                read_variant_seq = ""
                current_idx = variant_start_idx
                insertion_bases = 0
                context_bases = 0
                total_bases_needed = insertion_length + sequence_window
                
                while current_idx < len(read_positions) and len(read_variant_seq) < total_bases_needed:
                    if read_positions[current_idx] is not None:
                        read_variant_seq += read_seq[current_idx]
                        
                        # Count insertion bases (those without corresponding reference position)
                        if ref_positions[current_idx] is None:
                            insertion_bases += 1
                        else:
                            context_bases += 1
                    current_idx += 1
                
                # Check if we have enough bases for comparison
                if len(read_variant_seq) >= insertion_length:
                    # Calculate alignment scores against both alleles
                    alt_score = 0
                    ref_score = 0
                    
                    # Compare with alt allele
                    for i in range(min(len(read_variant_seq), len(alt_allele))):
                        if read_variant_seq[i] == alt_allele[i]:
                            alt_score += 1
                    
                    # Compare with ref allele
                    for i in range(min(len(read_variant_seq), len(ref_allele))):
                        if read_variant_seq[i] == ref_allele[i]:
                            ref_score += 1
                    
                    # Add weight to the insertion detection
                    if insertion_bases >= insertion_length * 0.8:  # Allow for some sequencing errors
                        alt_score += 2  # Bonus for having the right number of inserted bases
                    
                    # Calculate match percentages
                    alt_match_percent = alt_score / len(alt_allele)
                    ref_match_percent = ref_score / len(ref_allele)
                    
                    # Make the final call
                    if alt_match_percent >= 0.8 and insertion_bases >= insertion_length * 0.8:
                        return 'alt'
                    elif ref_match_percent >= 0.8 and insertion_bases == 0:
                        return 'ref'
                    
            # For deletions 
            else:
                deletion_length = len(ref_allele) - len(alt_allele)
                expected_ref_positions = list(range(pos, pos + len(ref_allele)))
                
                # Count missing bases in the read
                missing_bases = 0
                for exp_pos in expected_ref_positions:
                    if exp_pos not in ref_positions or read_positions[ref_positions.index(exp_pos)] is None:
                        missing_bases += 1
                
                # Calculate match scores for both alleles
                ref_match_score = len(ref_allele) - missing_bases
                alt_match_score = missing_bases
                
                # Add additional check for partial deletions
                if missing_bases >= deletion_length * 0.8:  # Allow some mismatches
                    return 'alt'
                elif ref_match_score >= len(ref_allele) * 0.8:  # Allow some mismatches
                    return 'ref'
        
        # Handle SNVs 
        else:
            query_pos = aligned_pairs[ref_pos_idx][0]
            if query_pos is None:
                return None
            
            read_base = read.query_sequence[query_pos]
            if read_base == ref_allele:
                return 'ref'
            elif read_base == alt_allele:
                return 'alt'
        
        return None
    
    def run_whatshap(self, bam_file: str, vcf_file: str, output_bam: str) -> Optional[str]:
        """Run WhatsHap for phasing variants and tagging reads - exact replica."""
        # Generate output file names based on input VCF
        vcf_prefix = os.path.splitext(os.path.basename(vcf_file))[0]
        output_dir = Path(vcf_file).parent
        phased_vcf = os.path.join(output_dir, f"{vcf_prefix}_Phased.vcf")
        phased_vcf_gz = f"{phased_vcf}.gz"
        whatshap_log = os.path.join(output_dir, "whatshap.log")

        # Open log file for WhatsHap output
        with open(whatshap_log, 'w') as log_file:
            # Index the input BAM file
            log_file.write(f"\nIndexing input BAM file: {bam_file}\n")
            index_command = f"samtools index {bam_file}"
            try:
                result = subprocess.run(index_command, shell=True, check=True, capture_output=True, text=True)
                log_file.write("BAM indexing completed successfully\n")
                if result.stderr:
                    log_file.write(f"Index stderr: {result.stderr}\n")
            except subprocess.CalledProcessError as e:
                log_file.write(f"Error indexing BAM file: {e.stderr}\n")
                raise

            # Run WhatsHap phase
            phase_command = f"whatshap phase -o {phased_vcf} --reference {self.reference_fasta} --internal-downsampling 23 --ignore-read-groups {vcf_file} {bam_file}"
            log_file.write("\nRunning WhatsHap phase command:\n")
            log_file.write(f"{phase_command}\n")
            
            try:
                phase_result = subprocess.run(phase_command, shell=True, check=True, capture_output=True, text=True)
                log_file.write("\nWhatsHap Phase stdout:\n")
                log_file.write(phase_result.stdout)
                
                if phase_result.stderr:
                    log_file.write("\nWhatsHap Phase stderr:\n")
                    log_file.write(phase_result.stderr)
            except subprocess.CalledProcessError as e:
                log_file.write("\nError running WhatsHap phase:\n")
                log_file.write(f"stdout: {e.stdout}\n")
                log_file.write(f"stderr: {e.stderr}\n")
                raise

            # Compress and index the phased VCF
            bgzip_command = f"bgzip -f {phased_vcf}"
            subprocess.run(bgzip_command, shell=True, check=True)
            
            index_command = f"tabix -p vcf {phased_vcf_gz}"
            subprocess.run(index_command, shell=True, check=True)

            # Run WhatsHap haplotag
            haplotag_command = f"whatshap haplotag --tag-supplementary -o {output_bam} --reference {self.reference_fasta} {phased_vcf_gz} {bam_file} --ignore-read-groups"
            log_file.write("\nRunning WhatsHap haplotag command:\n")
            log_file.write(f"{haplotag_command}\n")
            
            try:
                haplotag_result = subprocess.run(haplotag_command, shell=True, check=True, capture_output=True, text=True)
                log_file.write("\nWhatsHap Haplotag stdout:\n")
                log_file.write(haplotag_result.stdout)
                
                if haplotag_result.stderr:
                    log_file.write("\nWhatsHap Haplotag stderr:\n")
                    log_file.write(haplotag_result.stderr)
            except subprocess.CalledProcessError as e:
                log_file.write("\nError running WhatsHap haplotag:\n")
                log_file.write(f"stdout: {e.stdout}\n")
                log_file.write(f"stderr: {e.stderr}\n")
                raise

            # Index the output BAM
            log_file.write(f"\nIndexing output BAM file: {output_bam}\n")
            index_bam_command = f"samtools index {output_bam}"
            try:
                subprocess.run(index_bam_command, shell=True, check=True)
                log_file.write("Output BAM indexing completed successfully\n")
            except subprocess.CalledProcessError as e:
                log_file.write(f"Error indexing output BAM file: {e.stderr}\n")
                raise

        return phased_vcf_gz
    
    def validate_variants(self, input_bam, variants):
        """Validate that the specified variants match the expected patterns in the reads - exact replica."""
        validation_results = []
        warnings = []
        errors = []
        
        for variant in variants:
            ref_count = 0
            alt_count = 0
            unexpected_alleles = Counter()  # Track specific unexpected alleles
            total_covering_reads = 0
            
            # Expected alleles from VCF
            expected_ref = variant.ref.upper()
            expected_alt = variant.alts[0].upper() if variant.alts else ""
            
            # Fetch reads covering this variant position
            for read in input_bam.fetch(variant.chrom, variant.pos - 1, variant.pos):
                allele = self.get_allele(read, variant)
                if allele == 'ref':
                    ref_count += 1
                    total_covering_reads += 1
                elif allele == 'alt':
                    alt_count += 1
                    total_covering_reads += 1
                elif allele:  # Track unexpected alleles
                    query_pos = None
                    for qpos, rpos in read.get_aligned_pairs():
                        if rpos == variant.pos - 1 and qpos is not None:
                            query_pos = qpos
                            break
                    if query_pos is not None:
                        unexpected_base = read.query_sequence[query_pos]
                        unexpected_alleles[unexpected_base] += 1
                        total_covering_reads += 1
            
            if total_covering_reads == 0:
                errors.append(f"No reads cover variant at position {variant.pos}")
                continue
                
            # Calculate percentages
            ref_percentage = (ref_count / total_covering_reads * 100) if total_covering_reads > 0 else 0
            alt_percentage = (alt_count / total_covering_reads * 100) if total_covering_reads > 0 else 0
            unexpected_percentage = sum(unexpected_alleles.values()) / total_covering_reads * 100 if total_covering_reads > 0 else 0
            
            result = {
                'position': variant.pos,
                'total_reads': total_covering_reads,
                'ref_count': ref_count,
                'alt_count': alt_count,
                'ref_percentage': ref_percentage,
                'alt_percentage': alt_percentage,
                'unexpected_alleles': unexpected_alleles,
                'unexpected_percentage': unexpected_percentage,
                'expected_ref': expected_ref,
                'expected_alt': expected_alt
            }
            validation_results.append(result)
            
            # Validation checks
            if total_covering_reads < 30:
                errors.append(
                    f"Insufficient coverage at position {variant.pos}:\n"
                    f"Total reads: {total_covering_reads} (minimum 30 required)\n"
                    "Please verify that the correct variant position was provided."
                )
                continue
            
            # Check for unexpected alleles
            if unexpected_percentage > 0:
                unexpected_details = ", ".join(f"{base}: {count} reads ({count/total_covering_reads*100:.1f}%)" 
                                            for base, count in unexpected_alleles.most_common())
                errors.append(
                    f"Unexpected alleles found at position {variant.pos}:\n"
                    f"Expected ref: {expected_ref}, Expected alt: {expected_alt}\n"
                    f"Unexpected alleles: {unexpected_details}\n"
                    "Please verify that the correct variant position and alleles were provided."
                )
            
            # Check allele balance
            if ref_percentage > 95 or alt_percentage > 95:
                # Severe skew - error
                errors.append(
                    f"Variant at position {variant.pos} appears to be homozygous:\n"
                    f"Reference allele ({expected_ref}): {ref_count} reads ({ref_percentage:.1f}%)\n"
                    f"Alternate allele ({expected_alt}): {alt_count} reads ({alt_percentage:.1f}%)\n"
                    "This level of skew suggests this position is homozygous. A wrong variant is provided or allele dropout has occured.\n"
                    "Please verify that the correct variant position was provided."
                )
            elif ref_percentage > 80 or alt_percentage > 80:
                # Moderate skew - warning
                warnings.append(
                    f"Warning: Variant at position {variant.pos} shows skewed allele frequencies:\n"
                    f"Reference allele ({expected_ref}): {ref_count} reads ({ref_percentage:.1f}%)\n"
                    f"Alternate allele ({expected_alt}): {alt_count} reads ({alt_percentage:.1f}%)\n"
                    "While this is acceptable, please verify that the correct variant position was provided."
                )
        
        # Determine if validation passed
        validation_passed = len(errors) == 0
        
        return validation_passed, warnings, errors
    
    def quality_control(self, input_bam: str, vcf_file: str, output_bam: str, output_txt: str) -> int:
        """Perform quality control on input BAM file and filter reads - exact replica."""
        if not os.path.exists(input_bam + '.bai'):
            pysam.index(input_bam)
        vcf = pysam.VariantFile(vcf_file)
        variants = list(vcf.fetch())
        if len(variants) != 2:
            raise ValueError("VCF file should contain exactly two variants")

        var1, var2 = variants
        
        input_bam_file = pysam.AlignmentFile(input_bam, "rb")
        
        # Validate variants before proceeding
        variants_valid, warnings, errors = self.validate_variants(input_bam_file, variants)
        
        with open(output_txt, 'a') as f:
            print("\n=== Variant Validation ===", file=f)
            if variants_valid and not warnings:
                print("Variants were validated successfully", file=f)
                print(f"Both variants show expected heterozygous patterns in the reads.", file=f)
            elif variants_valid and warnings:
                print("Variants were validated successfully but with warnings:", file=f)
                for warning in warnings:
                    print(warning, file=f)
                    print(file=f)
            
            if errors:
                print("\nErrors:", file=f)
                for error in errors:
                    print(error, file=f)
                    print(file=f)
                input_bam_file.close()
                raise ValueError("Variant validation failed. Analysis cannot proceed.")
        
        output_bam_file = pysam.AlignmentFile(output_bam, "wb", template=input_bam_file)

        total_reads = input_bam_file.count()
        spanning_reads = 0
        clean_spanning_reads = 0
        high_quality_reads = 0

        with open(output_txt, 'a') as f:
            print("\n=== Quality Control ===", file=f)

        for read in input_bam_file.fetch():
            if read.reference_start <= var1.pos - 1 and read.reference_end >= var2.pos:
                spanning_reads += 1

                alleles = [self.get_allele(read, variant) for variant in variants]

                if alleles[0] in ['ref', 'alt'] and alleles[1] in ['ref', 'alt']:
                    clean_spanning_reads += 1

                    if read.mapping_quality >= 20:
                        high_quality_reads += 1
                        output_bam_file.write(read)

        input_bam_file.close()
        output_bam_file.close()

        with open(output_txt, 'a') as f:
            print(f"Total reads: {total_reads}", file=f)
            print(f"Spanning reads: {spanning_reads} ({spanning_reads/total_reads*100:.2f}%)", file=f)
            print(f"Clean spanning reads: {clean_spanning_reads} ({clean_spanning_reads/spanning_reads*100:.2f}% of spanning reads)", file=f)
            print(f"High quality spanning reads (MAPQ >= 20): {high_quality_reads} ({high_quality_reads/clean_spanning_reads*100:.2f}% of clean spanning reads)", file=f)

            if high_quality_reads > 50:
                print("\n*QC PASSED (>50 high quality spanning reads)*", file=f)
            else:
                print("*QC FAILED (<50 high quality spanning reads)*", file=f)

        return high_quality_reads
    
    def analyze_reads(self, bam_file: str, vcf_file: str, output_txt: str) -> None:
        """Analyze phased reads to determine variant relationships (Cis/Trans) - exact replica."""
        vcf = pysam.VariantFile(vcf_file)
        variants = list(vcf.fetch())
        if len(variants) != 2:
            raise ValueError("VCF file should contain exactly two variants")

        var1, var2 = variants

        bam = pysam.AlignmentFile(bam_file, "rb")

        total_reads = bam.count()
        ref_reads = 0
        alt_reads = 0
        ref_alt_reads = 0
        alt_ref_reads = 0

        for read in bam.fetch():
            alleles = [self.get_allele(read, variant) for variant in variants]

            if alleles[0] == 'ref' and alleles[1] == 'ref':
                ref_reads += 1
            elif alleles[0] == 'alt' and alleles[1] == 'alt':
                alt_reads += 1
            elif alleles[0] == 'ref' and alleles[1] == 'alt':
                ref_alt_reads += 1
            elif alleles[0] == 'alt' and alleles[1] == 'ref':
                alt_ref_reads += 1

        cis_reads = ref_reads + alt_reads
        trans_reads = ref_alt_reads + alt_ref_reads
        cis_percentage = cis_reads / total_reads * 100 if total_reads > 0 else 0
        trans_percentage = trans_reads / total_reads * 100 if total_reads > 0 else 0
        
        with open(output_txt, 'a') as f:
            print("\n=== Result ===", file=f)
            print(f"Total high quality spanning reads: {total_reads}", file=f)
            print("\nDetailed categorisation of reads:", file=f)
            print(f"Reads with ref allele for both variants (Cis): {ref_reads} ({ref_reads/total_reads*100:.2f}%)", file=f)
            print(f"Reads with alt allele for both variants (Cis): {alt_reads} ({alt_reads/total_reads*100:.2f}%)", file=f)
            print(f"Reads with ref for first, alt for second (Trans): {ref_alt_reads} ({ref_alt_reads/total_reads*100:.2f}%)", file=f)
            print(f"Reads with alt for first, ref for second (Trans): {alt_ref_reads} ({alt_ref_reads/total_reads*100:.2f}%)", file=f)
            print(f"\nCis reads (both ref or both alt): {cis_reads} ({cis_percentage:.2f}%)", file=f)
            print(f"Trans reads (one ref, one alt): {trans_reads} ({trans_percentage:.2f}%)", file=f)

            if cis_percentage > trans_percentage:
                print(f"\nCounting the reads determined the phase as Cis \n\nChimeric reads percentage: {trans_percentage:.2f}% ", file=f)
            else:
                print(f"\nCounting reads determined the phase as Trans \n\nChimeric reads percentage: {cis_percentage:.2f}% ", file=f)

        bam.close()
    
    def parse_vcf_phasing(self, phased_vcf_file: str) -> str:
        """Read the phased VCF file and determine if variants are in Cis or Trans - exact replica."""
        vcf = pysam.VariantFile(phased_vcf_file)
        variants = list(vcf.fetch())
        
        if len(variants) != 2:
            raise ValueError("VCF file should contain exactly two variants")
        
        # Extract genotype information
        var1_genotype = variants[0].samples[0]['GT']
        var2_genotype = variants[1].samples[0]['GT']
        
        # Determine phasing
        if var1_genotype == var2_genotype:
            return "Cis"
        else:
            return "Trans"
    
    def write_variant_comparison_results(self, dummy_vcf: str, variant_calling_vcf: str, output_file: str):
        """Write variant comparison results - exact replica from original."""
        # Ensure dummy VCF exists
        if not os.path.exists(dummy_vcf):
            print(f"Error: Dummy VCF file does not exist at {dummy_vcf}")
            return
        
        # Ensure variant calling VCF exists
        if not os.path.exists(variant_calling_vcf):
            print(f"Error: Variant calling VCF file does not exist at {variant_calling_vcf}")
            return

        try:
            matched_variants = self.compare_variants(dummy_vcf, variant_calling_vcf)
            
            with open(output_file, 'a') as f:  # Changed to append mode

                f.write("Variants to be phased:\n")
                # Read dummy VCF variants and store them
                variants = []
                dummy_vcf_content = []
                with open(dummy_vcf, 'r') as dummy_f:
                    for line in dummy_f:
                        if not line.startswith('#'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 5:
                                chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                                variants.append((chrom, int(pos), ref, alt))
                                dummy_vcf_content.append(f"{chrom}\t{pos}\t{ref}\t{alt}\n")

                # Create a set of them for matching
                dummy_variants_set = set()
                with open(dummy_vcf, 'r') as dummy_f:
                    for line in dummy_f:
                        if not line.startswith('#'):
                            parts = line.strip().split('\t')
                            if len(parts) >= 10:
                                chrom, pos, _, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
                                dummy_variants_set.add((chrom, int(pos), ref, alt))
                
                # Write the variants and distance information
                for line in dummy_vcf_content:
                    f.write(line)
                
                if len(variants) >= 2:
                    f.write(f"Distance between the variants is {variants[1][1] - variants[0][1]:,} bp\n")

                f.write("\n===Variant Calling===\n")           
                f.write("Variants called from the amplicon by Clair3:\n")
                f.write("CHROM\tPOS\tREF\tALT\tQUAL\tGT\n")
                # Read variant calling VCF variants
                vcf = pysam.VariantFile(variant_calling_vcf)
                for record in vcf.fetch():
                    if record.alts and record.alts[0] and record.alts[0] != ".":
                        qual = f"{record.qual:.1f}" if record.qual is not None else "."
                        gt = record.samples[0]['GT']
                        gt_str = "/".join(str(x) for x in gt) if gt is not None else "."
                        
                        # Check if this variant matches any in the dummy VCF
                        is_matched = (record.chrom, record.pos, record.ref, record.alts[0]) in dummy_variants_set
                        
                        # Add ">" for matched variants
                        indicator = "<-" if is_matched else " "
                        f.write(f"{record.chrom}\t{record.pos}\t{record.ref}\t{record.alts[0]}\t{qual}\t{gt_str} {indicator} \n")
                
                f.write("\nVariant Matching Results:\n")
                for i, variant in enumerate(matched_variants, 1):
                    if variant:
                        f.write(f"Variant {i}: {variant[0]} {variant[1]} {variant[2]} > {variant[3]} - FOUND\n")
                    else:
                        f.write(f"Variant {i}: NOT FOUND in variant calling VCF\n")
                
                # Summary statistics
                found_count = sum(1 for v in matched_variants if v)
                not_found_count = sum(1 for v in matched_variants if v is None)
        
        except Exception as e:
            print(f"Error during variant comparison: {e}")
            import traceback
            traceback.print_exc()
    
    def compare_variants(self, dummy_vcf: str, variant_calling_vcf: str) -> List[Tuple[str, int, str, str]]:
        """Compare variants between dummy VCF and variant calling VCF - exact replica."""
        # Read dummy VCF variants
        dummy_variants = []
        with open(dummy_vcf, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 5:
                        chrom, pos, _, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
                        dummy_variants.append((chrom, pos, ref, alt))
                    else:
                        print(f"Skipping malformed line in dummy VCF: {line.strip()}")

        # Read variant calling VCF variants
        variant_calling_variants = []
        try:
            vcf = pysam.VariantFile(variant_calling_vcf)
            for record in vcf.fetch():
                # Only add variants with valid alternate alleles (not "." or None)
                if record.alts and record.alts[0] and record.alts[0] != ".":
                    variant_calling_variants.append((record.chrom, record.pos, record.ref, record.alts[0]))
            
            print(f"Number of variants called from the amplicon: {len(variant_calling_variants)}")
        except Exception as e:
            print(f"Error reading variant calling VCF: {e}")
            return []

        # Compare variants
        matched_variants = []
        for dummy_var in dummy_variants:
            found = False
            for vc_var in variant_calling_variants:
                if (dummy_var[0] == vc_var[0] and  # chromosome
                    dummy_var[1] == vc_var[1] and   # position
                    dummy_var[2] == vc_var[2] and   # reference allele
                    dummy_var[3] == vc_var[3]):     # alternate allele
                    found = True
                    matched_variants.append(dummy_var)
                    break
            
            if not found:
                matched_variants.append(None)
        
        return matched_variants