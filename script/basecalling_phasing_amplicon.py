import pysam
import os
import sys
from typing import List, Tuple
from pathlib import Path
from quality_control import get_allele

def compare_variants(dummy_vcf: str, variant_calling_vcf: str) -> List[Tuple[str, int, str, str]]:
    """
    Compare variants between dummy VCF and variant calling VCF.
    
    Args:
        dummy_vcf (str): Path to the dummy VCF file
        variant_calling_vcf (str): Path to the variant calling VCF file
    
    Returns:
        List of variants found in both files, with details about matching status
    """
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

def write_variant_comparison_results(dummy_vcf: str, variant_calling_vcf: str, output_file: str):
    """
    Write variant comparison results to the report file.
    
    Args:
        report_file: File handle for the report
        dummy_vcf (str): Path to the dummy VCF file
        variant_calling_vcf (str): Path to the variant calling VCF file
    """    
    # Ensure dummy VCF exists
    if not os.path.exists(dummy_vcf):
        print(f"Error: Dummy VCF file does not exist at {dummy_vcf}")
        return
    
    # Ensure variant calling VCF exists
    if not os.path.exists(variant_calling_vcf):
        print(f"Error: Variant calling VCF file does not exist at {variant_calling_vcf}")
        return

    try:
        matched_variants = compare_variants(dummy_vcf, variant_calling_vcf)
        
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

def analyze_reads(bam_file, vcf_file):
    """
    Analyze phased reads to determine variant relationships (Cis/Trans).
    """
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
        alleles = [get_allele(read, variant) for variant in variants]

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
    print("\n=== Result ===")
    print(f"Total high quality spanning reads: {total_reads}")
    print("\nDetailed categorisation of reads:")
    print(f"Reads with ref allele for both variants (Cis): {ref_reads} ({ref_reads/total_reads*100:.2f}%)")
    print(f"Reads with alt allele for both variants (Cis): {alt_reads} ({alt_reads/total_reads*100:.2f}%)")
    print(f"Reads with ref for first, alt for second (Trans): {ref_alt_reads} ({ref_alt_reads/total_reads*100:.2f}%)")
    print(f"Reads with alt for first, ref for second (Trans): {alt_ref_reads} ({alt_ref_reads/total_reads*100:.2f}%)")
    print(f"\nCis reads (both ref or both alt): {cis_reads} ({cis_percentage:.2f}%)")
    print(f"Trans reads (one ref, one alt): {trans_reads} ({trans_percentage:.2f}%)")

    if cis_percentage > trans_percentage:
        print(f"\nCounting the reads determined the phase as Cis \n\nChimeric reads percentage: {trans_percentage:.2f}% ")
    else:
        print(f"\nCounting reads determined the phase as Trans \n\nChimeric reads percentage: {cis_percentage:.2f}% ")

    bam.close()

def parse_vcf_phasing(phased_vcf_file):
    """
    Read the phased VCF file and determine if variants are in Cis or Trans.
    
    Args:
        phased_vcf_file (str): Path to the phased VCF file (can be gzipped)
    
    Returns:
        str: Interpretation of variant phasing (Cis or Trans)
    """
    vcf = pysam.VariantFile(phased_vcf_file)  # pysam automatically handles gzipped files
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

def main():
    """
    Main function to run the complete analysis pipeline.
    """
    if len(sys.argv) != 3:
        print("Usage: python script.py <input.bam> <input.vcf>")
        sys.exit(1)
    
    bam_file = sys.argv[1]
    vcf_file = sys.argv[2]
    episode = os.path.splitext(os.path.basename(bam_file))[0]
    
    # Create output names
    output_bam = bam_file.removesuffix('.bam') + "_phased.bam"
    output_txt = bam_file.removesuffix('.bam') + "_report.txt"

    with open(output_txt, 'a') as output_file:
        original_stdout = sys.stdout
        original_stderr = sys.stderr
        sys.stdout = output_file
        sys.stderr = output_file

        try:
            dummy_vcf = vcf_file
            
            # Look for the variant calling VCF
            variant_calling_vcf = os.path.join(Path(bam_file).parent._str, f"{episode}.wf_snp.vcf.gz")
            
            if os.path.exists(dummy_vcf) and os.path.exists(variant_calling_vcf):
                write_variant_comparison_results(dummy_vcf, variant_calling_vcf, output_txt)
            else:
                print("Skipping variant comparison due to missing files.")
                print(f"Dummy VCF exists: {os.path.exists(dummy_vcf)}")
                print(f"Variant calling VCF exists: {os.path.exists(variant_calling_vcf)}")
                print(f"Looking for variant calling VCF at: {variant_calling_vcf}")

            # Quality control is now handled by separate script
            clean_span_bam = os.path.join(Path(bam_file).parent._str, f"clean-span-hq.bam")
            
            if os.path.exists(clean_span_bam):
                # WhatsHap is now run from the shell script
                analyze_reads(output_bam, vcf_file)

                # Analyze WhatsHap phased VCF (gzipped)
                phased_vcf_gz = os.path.join(Path(bam_file).parent._str, f"{episode}_Phased.vcf.gz")
                if os.path.exists(phased_vcf_gz):
                    whatshap_phasing = parse_vcf_phasing(phased_vcf_gz)
                    print(" "*40)
                    print(f"WhatsHap analysis determined the phase as {whatshap_phasing}")
                else:
                    print("WhatsHap phased VCF not found.")

                # Look for and analyze HapCUT2 phased VCF
                hapcut2_vcf = os.path.join(Path(bam_file).parent._str, f"hap2cut_{episode}.phased.VCF")
                if os.path.exists(hapcut2_vcf):
                    try:
                        hapcut2_phasing = parse_vcf_phasing(hapcut2_vcf)
                        print(f"\nHapCUT2 analysis determined the phase as {hapcut2_phasing}\n")
                    except Exception as e:
                        print(f"Error analyzing HapCUT2 VCF: {e}")
                else:
                    print(f"HapCUT2 phased VCF not found at: {hapcut2_vcf}")
            else:
                print("No high quality spanning reads found. Analysis cannot proceed.")

        except Exception as e:
            print(f"An error occurred: {e}")
            import traceback
            traceback.print_exc()
        finally:
            sys.stdout = original_stdout
            sys.stderr = original_stderr

    print(f"Analysis complete. Results written to {output_txt}")

if __name__ == "__main__":
    main()