"""Main pipeline orchestrator."""

import pandas as pd
import subprocess
import pysam
import os
from pathlib import Path
from typing import Dict, List, Any, Optional
import logging
import shutil

from ..utils.config import ConfigManager
from ..utils.file_utils import setup_output_directory, cleanup_files, copy_template_file
from ..utils.validators import validate_bam_file, validate_vcf_file
from .qc_analyzer import QCAnalyzer
from .phasing_analyzer import PhasingAnalyzer
from .variant_processor import VariantProcessor
from .xml_generator import XMLGenerator

logger = logging.getLogger(__name__)


class AmpliconPipeline:
    """Main pipeline orchestrator for ONT amplicon analysis."""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.qc_analyzer = QCAnalyzer(config)
        self.phasing_analyzer = PhasingAnalyzer(config)
        self.variant_processor = VariantProcessor(config)
        self.xml_generator = XMLGenerator(config)
    
    def run(self, run_id: str, input_dir: Path, output_dir: Path, 
            sample_sheet: Path, clean: bool = False) -> List[Dict[str, Any]]:
        """
        Run the complete pipeline for all samples.
        
        Args:
            run_id: Run identifier
            input_dir: Input directory containing BAM files
            output_dir: Output directory for results
            sample_sheet: Path to sample sheet CSV
            clean: Whether to clean existing output directories
        
        Returns:
            List of processing results for each sample
        """
        logger.info(f"Starting pipeline run: {run_id}")
        
        # Read and process sample sheet
        df = self._prepare_sample_sheet(sample_sheet, input_dir)
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        results = []
        
        for _, row in df.iterrows():
            try:
                result = self._process_sample(
                    row=row,
                    run_id=run_id,
                    input_dir=input_dir,
                    output_dir=output_dir,
                    clean=clean
                )
                results.append(result)
                
            except Exception as e:
                logger.error(f"Failed to process {row['Episode']}: {e}")
                results.append({
                    'episode': row['Episode'],
                    'barcode': row['Barcode'],
                    'success': False,
                    'error': str(e)
                })
        
        # Generate XML files after all samples are processed
        if self.config.get('output', {}).get('generate_xml', True):
            self.xml_generator.generate_xml_files(run_id, input_dir, output_dir)
        
        logger.info(f"Pipeline completed. Processed {len(results)} samples")
        return results
    
    def _prepare_sample_sheet(self, sample_sheet: Path, input_dir: Path) -> pd.DataFrame:
        """Prepare and validate sample sheet."""
        df = pd.read_csv(sample_sheet)
        
        # Process barcode format
        df['Barcode'] = df['Barcode'].apply(self._format_barcode)
        
        # Convert Episode and EpisodeWES to uppercase
        df['Episode'] = df['Episode'].str.upper()
        df['EpisodeWES'] = df['EpisodeWES'].str.upper()
        
        # Clean whitespace
        for col in df.columns:
            if df[col].dtype == 'object':
                df[col] = df[col].astype(str).str.strip()
        
        return df
    
    def _format_barcode(self, barcode: Any) -> str:
        """Format barcode with leading zeros."""
        if isinstance(barcode, (int, float)):
            return f"barcode{int(barcode):02d}"
        elif isinstance(barcode, str) and barcode.isdigit():
            return f"barcode{int(barcode):02d}"
        else:
            return str(barcode)
    
    def _process_sample(self, row: pd.Series, run_id: str, input_dir: Path, 
                       output_dir: Path, clean: bool) -> Dict[str, Any]:
        """Process a single sample."""
        episode = row['Episode']
        barcode = row['Barcode']
        coordinate = row['Coordinate']
        variant1 = row['Variant1']
        variant2 = row['Variant2']
        
        logger.info(f"Processing {episode} ({barcode})")
        
        # Setup output directory
        sample_output_dir = setup_output_directory(output_dir, barcode, clean)
        
        # Determine analysis type
        has_variant1 = self._has_valid_variant(variant1)
        has_variant2 = self._has_valid_variant(variant2)
        
        if has_variant1 and has_variant2:
            analysis_type = "phasing"
        elif has_variant1:
            analysis_type = "single_variant_qc"
        else:
            analysis_type = "qc_only"
        
        logger.info(f"Analysis type for {episode}: {analysis_type}")
        
        try:
            # Step 1: Merge and extract BAM files
            bam_file = self._process_bam_files(
                input_dir, barcode, episode, coordinate, sample_output_dir
            )
            
            # Step 2: Create BED file
            bed_file = self._create_bed_file(coordinate, episode, sample_output_dir)
            
            # Step 3: Create VCF file if variants provided
            vcf_file = None
            if has_variant1:
                vcf_file = self._create_variant_vcf(
                    episode, variant1, variant2, sample_output_dir
                )
            
            # Step 4: Run variant calling
            try:
                called_vcf = self._run_variant_calling(
                    bam_file, bed_file, episode, sample_output_dir
                )
            except Exception as e:
                logger.error(f"Variant calling failed for {episode}: {e}")
                # Continue with analysis even if variant calling fails
                called_vcf = None
            
            # Step 5: Run analysis based on type
            if analysis_type == "phasing":
                # Run phasing analysis
                self._run_hapcut2(bam_file, vcf_file, episode, sample_output_dir)
                phasing_results = self.phasing_analyzer.analyze_phasing(
                    bam_file, vcf_file, episode, sample_output_dir
                )
                self.phasing_analyzer.write_phasing_report(
                    phasing_results, sample_output_dir / f"{episode}_report.txt"
                )
            else:
                # Run QC analysis
                qc_results = self.qc_analyzer.analyze_amplicon(bam_file, bed_file, episode)
                
                # Add variant comparison if applicable
                report_file = sample_output_dir / f"{episode}_report.txt"
                self.qc_analyzer.write_qc_report(qc_results, report_file)
                
                if vcf_file and called_vcf and called_vcf.exists():
                    try:
                        comparison_results = self.variant_processor.compare_variants(vcf_file, called_vcf)
                        user_variants = self.variant_processor._read_variants_from_vcf(vcf_file)
                        called_variants = self.variant_processor._read_variants_from_called_vcf(called_vcf)
                        
                        # Write comparison to report
                        self.variant_processor.write_variant_comparison_report(
                            comparison_results, called_variants, report_file, user_variants
                        )
                    except Exception as e:
                        logger.warning(f"Variant comparison failed for {episode}: {e}")
                        with open(report_file, 'a') as f:
                            f.write(f"\nVariant comparison failed: {e}\n")
                elif vcf_file:
                    with open(report_file, 'a') as f:
                        f.write("\nVariant calling failed - no comparison available\n")
            
            # Step 6: Upload to S3
            if self.config.get('output', {}).get('upload_to_s3', True):
                self.xml_generator.upload_results_to_s3(run_id, barcode, sample_output_dir)
            
            # Step 7: Cleanup
            self._cleanup_intermediate_files(sample_output_dir, episode, analysis_type)
            
            return {
                'episode': episode,
                'barcode': barcode,
                'analysis_type': analysis_type,
                'success': True,
                'output_dir': sample_output_dir
            }
            
        except Exception as e:
            logger.error(f"Error processing {episode}: {e}")
            raise
    
    def _has_valid_variant(self, variant: str) -> bool:
        """Check if variant string is valid."""
        if pd.isna(variant) or not variant or str(variant).strip().lower() in ['', 'nan', 'none']:
            return False
        return 'chr' in str(variant).lower()
    
    def _process_bam_files(self, input_dir: Path, barcode: str, episode: str, 
                          coordinate: str, output_dir: Path) -> Path:
        """Merge BAM files and extract amplicon region."""
        logger.info(f"Processing BAM files for {episode}")
        
        bam_dir = input_dir / "bam_pass" / barcode
        if not bam_dir.exists():
            raise FileNotFoundError(f"BAM directory not found: {bam_dir}")
        
        bam_files = list(bam_dir.glob("*.bam"))
        if not bam_files:
            raise FileNotFoundError(f"No BAM files found in {bam_dir}")
        
        temp_bam = output_dir / "temp.bam"
        output_bam = output_dir / f"{episode}.bam"
        
        # Merge BAM files
        if len(bam_files) == 1:
            # Single file, just copy
            shutil.copy2(bam_files[0], temp_bam)
        else:
            # Multiple files, merge
            pysam.merge("-@", "6", "-u", str(temp_bam), *[str(f) for f in bam_files])
        
        # Sort merged BAM
        pysam.sort("-@", "6", "-o", str(temp_bam), str(temp_bam))
        pysam.index(str(temp_bam))
        
        # Extract amplicon region
        with pysam.AlignmentFile(str(temp_bam), "rb") as infile:
            with pysam.AlignmentFile(str(output_bam), "wb", template=infile) as outfile:
                chrom, start_end = coordinate.split(':')
                start, end = map(int, start_end.split('-'))
                
                for read in infile.fetch(chrom, start, end):
                    outfile.write(read)
        
        # Index output BAM
        pysam.index(str(output_bam))
        
        # Cleanup temp files
        cleanup_files([temp_bam, Path(str(temp_bam) + ".bai")])
        
        logger.info(f"Created amplicon BAM: {output_bam}")
        return output_bam
    
    def _create_bed_file(self, coordinate: str, episode: str, output_dir: Path) -> Path:
        """Create BED file from coordinate string."""
        bed_file = output_dir / f"{episode}_coordinate.bed"
        
        chrom, start_end = coordinate.split(':')
        start, end = start_end.split('-')
        
        with open(bed_file, 'w') as f:
            f.write(f"{chrom}\t{start}\t{end}\n")
        
        logger.debug(f"Created BED file: {bed_file}")
        return bed_file
    
    def _create_variant_vcf(self, episode: str, variant1: str, variant2: str, 
                           output_dir: Path) -> Path:
        """Create VCF file with specified variants."""
        template_vcf = Path(__file__).parent.parent / "data" / "dummy.vcf"
        output_vcf = output_dir / f"{episode}.vcf"
        
        self.variant_processor.create_variant_vcf(
            template_vcf, output_vcf, episode, variant1, variant2
        )
        
        return output_vcf
    
    def _run_variant_calling(self, bam_file: Path, bed_file: Path, 
                           episode: str, output_dir: Path) -> Path:
        """Run Clair3 variant calling."""
        logger.info(f"Running variant calling for {episode}")
        
        clair3_path = self.config['paths']['clair3_path']
        reference = self.config['paths']['reference_genome']
        clair3_config = self.config['variant_calling']['clair3']
        
        # Validate required files exist
        clair3_script = Path(clair3_path) / "run_clair3.sh"
        if not clair3_script.exists():
            raise FileNotFoundError(f"Clair3 script not found: {clair3_script}")
        
        if not Path(reference).exists():
            raise FileNotFoundError(f"Reference genome not found: {reference}")
        
        model_path = Path(clair3_path) / "models" / self.config['paths']['clair3_model']
        if not model_path.exists():
            raise FileNotFoundError(f"Clair3 model not found: {model_path}")
        
        clair3_output = output_dir / "variant_calling_output"
        log_file = output_dir / "clair3.log"
        
        # Build command
        cmd = [
            f"{clair3_path}/run_clair3.sh",
            f"--bam_fn={bam_file}",
            f"--bed_fn={bed_file}",
            f"--ref_fn={reference}",
            f"--threads={clair3_config['threads']}",
            f"--platform={clair3_config['platform']}",
            f"--model_path={clair3_path}/models/{self.config['paths']['clair3_model']}",
            f"--sample_name={episode}",
            f"--output={clair3_output}",
            f"--min_coverage={clair3_config['min_coverage']}",
            f"--var_pct_full={clair3_config['var_pct_full']}",
            f"--ref_pct_full={clair3_config['ref_pct_full']}",
            f"--var_pct_phasing={clair3_config['var_pct_phasing']}"
        ]
        
        if clair3_config.get('enable_phasing'):
            cmd.append("--enable_phasing")
        if clair3_config.get('use_whatshap'):
            cmd.append("--use_whatshap_for_final_output_phasing")
        if clair3_config.get('remove_intermediate'):
            cmd.append("--remove_intermediate_dir")
        
        logger.info(f"Running Clair3 command: {' '.join(cmd)}")
        
        # Prepare environment - activate conda environment if needed
        env = os.environ.copy()
        
        # Run Clair3 with detailed logging
        with open(log_file, 'w') as log:
            log.write(f"Clair3 command: {' '.join(cmd)}\n")
            log.write(f"Working directory: {output_dir}\n")
            log.write(f"Environment PATH: {env.get('PATH', 'Not set')}\n\n")
            log.flush()
            
            # Activate Clair3 conda environment and run command
            bash_cmd = f"source $(conda info --base)/etc/profile.d/conda.sh && conda activate clair3 && cd {output_dir} && {' '.join(cmd)}"
            result = subprocess.run(['bash', '-c', bash_cmd], stdout=log, stderr=subprocess.STDOUT, text=True, env=env)
            
        # Check result and log details
        if result.returncode != 0:
            with open(log_file, 'r') as log:
                log_content = log.read()
            logger.error(f"Clair3 failed with return code {result.returncode}")
            logger.error(f"Clair3 log content:\n{log_content}")
            raise RuntimeError(f"Clair3 variant calling failed with return code {result.returncode}. Check {log_file} for details.")
        
        # Check if output files exist
        source_vcf = clair3_output / "merge_output.vcf.gz"
        source_tbi = clair3_output / "merge_output.vcf.gz.tbi"
        
        if not source_vcf.exists():
            with open(log_file, 'r') as log:
                log_content = log.read()
            logger.error(f"Clair3 output file not found: {source_vcf}")
            logger.error(f"Clair3 log content:\n{log_content}")
            raise FileNotFoundError(f"Clair3 output file not found: {source_vcf}. Check {log_file} for details.")
        
        # Copy output files
        dest_vcf = output_dir / f"{episode}.wf_snp.vcf.gz"
        dest_tbi = output_dir / f"{episode}.wf_snp.vcf.gz.tbi"
        
        shutil.copy2(source_vcf, dest_vcf)
        if source_tbi.exists():
            shutil.copy2(source_tbi, dest_tbi)
        
        # Cleanup Clair3 output directory
        if clair3_config.get('remove_intermediate'):
            cleanup_files([clair3_output])
        
        logger.info(f"Variant calling completed: {dest_vcf}")
        return dest_vcf
    
    def _run_hapcut2(self, bam_file: Path, vcf_file: Path, episode: str, output_dir: Path) -> None:
        """Run HapCUT2 for phasing."""
        logger.info(f"Running HapCUT2 for {episode}")
        
        hapcut2_path = self.config['paths']['hapcut2_path']
        reference = self.config['paths']['reference_genome']
        
        fragment_file = output_dir / f"fragment_{episode}"
        hapcut2_output = output_dir / f"hap2cut_{episode}"
        log_file = output_dir / "HapCUT2.log"
        
        # Extract HAIRS
        extract_cmd = [
            f"{hapcut2_path}/extractHAIRS",
            "--ont", "1",
            "--bam", str(bam_file),
            "--VCF", str(vcf_file),
            "--out", str(fragment_file),
            "--indels", "1",
            "--ref", reference
        ]
        
        with open(log_file, 'w') as log:
            result = subprocess.run(extract_cmd, stdout=log, stderr=log)
            if result.returncode != 0:
                raise RuntimeError(f"extractHAIRS failed for {episode}")
        
        # Run HAPCUT2
        hapcut2_cmd = [
            f"{hapcut2_path}/HAPCUT2",
            "--fragments", str(fragment_file),
            "--VCF", str(vcf_file),
            "--output", str(hapcut2_output)
        ]
        
        with open(log_file, 'a') as log:
            result = subprocess.run(hapcut2_cmd, stdout=log, stderr=log)
            if result.returncode != 0:
                raise RuntimeError(f"HAPCUT2 failed for {episode}")
        
        logger.info(f"HapCUT2 completed for {episode}")
    
    def _cleanup_intermediate_files(self, output_dir: Path, episode: str, analysis_type: str) -> None:
        """Clean up intermediate files."""
        if not self.config.get('output', {}).get('cleanup_intermediate', True):
            return
        
        cleanup_files_list = []
        
        if analysis_type == "phasing":
            cleanup_files_list.extend([
                output_dir / f"{episode}.vcf",
                output_dir / f"fragment_{episode}",
                output_dir / f"hap2cut_{episode}",
                output_dir / "clean-span-hq.bam",
                output_dir / "clean-span-hq.bam.bai"
            ])
        
        cleanup_files(cleanup_files_list, ignore_errors=True)
        logger.debug(f"Cleaned up intermediate files for {episode}")