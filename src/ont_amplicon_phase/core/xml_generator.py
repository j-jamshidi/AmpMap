"""XML generator for IGV visualization files."""

import subprocess
import os
from pathlib import Path
from typing import Dict, Any, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


class XMLGenerator:
    """Generates IGV XML files for visualization of results."""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.s3_bucket = config.get('aws', {}).get('s3_bucket', 'nswhp-gaia-poc-pl')
        self.s3_prefix = config.get('aws', {}).get('s3_prefix', 'ONT')
    
    def generate_xml_files(self, run_id: str, input_dir: Path, output_dir: Path) -> None:
        """Generate XML files for all samples in the run."""
        logger.info(f"Generating XML files for run: {run_id}")
        
        info_file = input_dir / f"{run_id}.info"
        if not info_file.exists():
            logger.error(f"Info file not found: {info_file}")
            return
        
        # Process each sample
        with open(info_file, 'r') as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split(',')
                    if len(parts) >= 7:
                        batch, barcode, episode, coordinate, variant1, variant2, episode_wes = parts[:7]
                        
                        try:
                            self._generate_sample_xml(
                                run_id, barcode, episode, coordinate, 
                                variant1, variant2, episode_wes, output_dir
                            )
                        except Exception as e:
                            logger.error(f"Failed to generate XML for {episode}: {e}")
    
    def _generate_sample_xml(self, run_id: str, barcode: str, episode: str, 
                           coordinate: str, variant1: str, variant2: str, 
                           episode_wes: str, output_dir: Path) -> None:
        """Generate XML file for a single sample."""
        logger.debug(f"Generating XML for {episode} ({barcode})")
        
        sample_output_dir = output_dir / barcode
        xml_file = sample_output_dir / f"{episode}.xml"
        
        # Determine template based on WES availability
        if episode_wes == "NA" or not episode_wes.strip():
            template_file = Path(__file__).parent.parent / "data" / "solo_LR.xml"
        else:
            # Check if WES data is available
            wes_available = self._check_wes_availability(episode_wes)
            if wes_available:
                template_file = Path(__file__).parent.parent / "data" / "solo_LR_SR.xml"
            else:
                template_file = Path(__file__).parent.parent / "data" / "solo_LR.xml"
        
        # Copy template
        with open(template_file, 'r') as src:
            xml_content = src.read()
        
        # Replace placeholders
        xml_content = self._replace_placeholders(
            xml_content, run_id, barcode, episode, coordinate, 
            variant1, variant2, episode_wes
        )
        
        # Write XML file
        with open(xml_file, 'w') as f:
            f.write(xml_content)
        
        logger.debug(f"Generated XML file: {xml_file}")
    
    def _check_wes_availability(self, episode_wes: str) -> bool:
        """Check if WES data is available for the episode."""
        if not episode_wes or episode_wes.upper() == "NA":
            return False
        
        try:
            # Combine sample files as in original script
            software_dir = self.config.get('paths', {}).get('software_dir', '/EBSDataDrive/software')
            sample_files = self.config.get('basespace', {}).get('sample_files', [
                f"{software_dir}/sample_ran.txt",
                f"{software_dir}/sample_ran_CRE_BS.txt"
            ])
            combined_file = self.config.get('basespace', {}).get('combined_file', f"{software_dir}/sample_SR.txt")
            
            # Combine files
            with open(combined_file, 'w') as outfile:
                for sample_file in sample_files:
                    try:
                        with open(sample_file, 'r') as infile:
                            outfile.write(infile.read())
                    except FileNotFoundError:
                        logger.warning(f"Sample file not found: {sample_file}")
            
            # Check if episode exists in combined file
            with open(combined_file, 'r') as f:
                for line in f:
                    if episode_wes.upper() in line.upper():
                        parts = line.strip().split('\t')
                        if len(parts) >= 6 and parts[5].strip():  # wesrun field
                            return True
            
            return False
            
        except Exception as e:
            logger.error(f"Error checking WES availability for {episode_wes}: {e}")
            return False
    
    def _replace_placeholders(self, xml_content: str, run_id: str, barcode: str, 
                            episode: str, coordinate: str, variant1: str, 
                            variant2: str, episode_wes: str) -> str:
        """Replace placeholders in XML template."""
        
        # Replace episode ID
        xml_content = xml_content.replace("C1_IID", episode)
        
        # Determine coordinate for visualization
        viz_coordinate = self._calculate_visualization_coordinate(
            coordinate, variant1, variant2
        )
        xml_content = xml_content.replace("CHRSTARTEND", viz_coordinate)
        
        # Generate S3 presigned URLs
        urls = self._generate_s3_urls(run_id, barcode, episode, variant2)
        
        # Replace URL placeholders
        xml_content = xml_content.replace("WF_SNP_VCF", urls['wf_vcf'])
        xml_content = xml_content.replace("PHASED_VCF", urls['phased_vcf'])
        xml_content = xml_content.replace("C1_LR_BAM_URL", urls['lr_bam'])
        xml_content = xml_content.replace("C1_LR_index_URL", urls['lr_bai'])
        
        # Replace WES URLs if available
        if "C1_SR_BAM_URL" in xml_content:
            wes_urls = self._get_wes_urls(episode_wes)
            if wes_urls:
                xml_content = xml_content.replace("C1_SR_BAM_URL", wes_urls['bam'])
                xml_content = xml_content.replace("C1_SR_index_URL", wes_urls['bai'])
        
        return xml_content
    
    def _calculate_visualization_coordinate(self, coordinate: str, 
                                          variant1: str, variant2: str) -> str:
        """Calculate coordinate for visualization based on variants."""
        # If variants are provided, use them to define the region
        if self._has_variant(variant1) or self._has_variant(variant2):
            if not self._has_variant(variant2):
                variant2 = variant1
            
            # Extract positions
            chr1, pos1 = self._parse_variant_position(variant1)
            chr2, pos2 = self._parse_variant_position(variant2)
            
            if chr1 != chr2:
                logger.warning("Variants are on different chromosomes")
                return coordinate
            
            # Create region around variants
            start = min(pos1, pos2) - 100
            end = max(pos1, pos2) + 100
            return f"{chr1}:{start}-{end}"
        
        return coordinate
    
    def _has_variant(self, variant: str) -> bool:
        """Check if variant string contains a valid variant."""
        return variant and "chr" in variant.lower()
    
    def _parse_variant_position(self, variant: str) -> Tuple[str, int]:
        """Parse chromosome and position from variant string."""
        if ":" in variant:
            parts = variant.split(":")
            chrom = parts[0]
            pos_part = parts[1].split()[0] if " " in parts[1] else parts[1]
            pos = int(pos_part)
            return chrom, pos
        return "chr1", 0
    
    def _generate_s3_urls(self, run_id: str, barcode: str, episode: str, 
                         variant2: str) -> Dict[str, str]:
        """Generate S3 presigned URLs for files."""
        base_path = f"s3://{self.s3_bucket}/{self.s3_prefix}/{run_id}/{barcode}"
        
        # Determine file names based on analysis type
        if self._has_variant(variant2):
            # Phasing analysis
            lr_bam = f"{base_path}/{episode}_phased.bam"
            lr_bai = f"{base_path}/{episode}_phased.bam.bai"
            phased_vcf = f"{base_path}/{episode}_Phased.vcf"
        else:
            # QC analysis
            lr_bam = f"{base_path}/{episode}_QC.bam"
            lr_bai = f"{base_path}/{episode}_QC.bam.bai"
            phased_vcf = f"{base_path}/{episode}.vcf"
        
        wf_vcf = f"{base_path}/{episode}.wf_snp.vcf.gz"
        
        # Generate presigned URLs
        urls = {
            'lr_bam': self._generate_presigned_url(lr_bam),
            'lr_bai': self._generate_presigned_url(lr_bai),
            'phased_vcf': self._generate_presigned_url(phased_vcf),
            'wf_vcf': self._generate_presigned_url(wf_vcf)
        }
        
        return urls
    
    def _generate_presigned_url(self, s3_path: str, expire_seconds: int = 604800) -> str:
        """Generate presigned URL for S3 object."""
        try:
            cmd = ["aws", "s3", "presign", s3_path, "--expires-in", str(expire_seconds)]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            url = result.stdout.strip()
            
            # Escape URL for XML
            url = url.replace("/", "\\/")
            url = url.replace("&", "&amp;")
            
            return url
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to generate presigned URL for {s3_path}: {e}")
            return s3_path
    
    def _get_wes_urls(self, episode_wes: str) -> Optional[Dict[str, str]]:
        """Get WES BAM URLs from BaseSpace."""
        if not episode_wes or episode_wes.upper() == "NA":
            return None
        
        try:
            # Get BaseSpace configuration
            config_name = self.config.get('basespace', {}).get('config_name', 'POWH')
            
            # Get sample ID from combined file
            software_dir = self.config.get('paths', {}).get('software_dir', '/EBSDataDrive/software')
            combined_file = f"{software_dir}/sample_SR.txt"
            
            sid = None
            with open(combined_file, 'r') as f:
                for line in f:
                    if episode_wes.upper() in line.upper():
                        parts = line.strip().split('\t')
                        if len(parts) >= 1:
                            sid = parts[0].upper()
                            break
            
            if not sid:
                logger.warning(f"Sample ID not found for {episode_wes}")
                return None
            
            # Get dataset ID
            cmd = ["bs", "-c", config_name, "list", "datasets", "--input-biosample", sid, 
                   "--not-type", "illumina.fastq.v1.8", "--sort-by", "AppSession.DateCreated", "--terse"]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            dataset_ids = result.stdout.strip().split('\n')
            if not dataset_ids or not dataset_ids[-1]:
                return None
            
            datid = dataset_ids[-1]
            
            # Get BAM and BAI file IDs
            bam_cmd = ["bs", "dataset", "-c", config_name, "content", "--id", datid, "--extension", "bam", "--terse"]
            bai_cmd = ["bs", "dataset", "-c", config_name, "content", "--id", datid, "--extension", "bam.bai", "--terse"]
            
            bam_result = subprocess.run(bam_cmd, capture_output=True, text=True, check=True)
            bai_result = subprocess.run(bai_cmd, capture_output=True, text=True, check=True)
            
            bam_id = bam_result.stdout.strip()
            bai_id = bai_result.stdout.strip()
            
            if not bam_id or not bai_id:
                return None
            
            # Get presigned URLs
            bam_url = self._get_basespace_presigned_url(bam_id, config_name)
            bai_url = self._get_basespace_presigned_url(bai_id, config_name)
            
            if bam_url and bai_url:
                return {
                    'bam': bam_url,
                    'bai': bai_url
                }
            
            return None
            
        except Exception as e:
            logger.error(f"Error getting WES URLs for {episode_wes}: {e}")
            return None
    
    def _get_basespace_presigned_url(self, file_id: str, config_name: str) -> Optional[str]:
        """Get presigned URL from BaseSpace file ID."""
        try:
            # Get AWS link
            cmd = ["bs", "-c", config_name, "file", "link", "-i", file_id]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            aws_link = result.stdout.strip()
            
            # Get presigned URL using wget
            wget_cmd = ["wget", "--save-headers", "--max-redirect=0", "-O", "-", aws_link]
            wget_result = subprocess.run(wget_cmd, capture_output=True, text=True)
            
            # Extract Location header
            for line in wget_result.stderr.split('\n'):
                if 'Location' in line:
                    url = line.split()[-1]
                    # Escape URL for XML
                    url = url.replace("/", "\\/")
                    url = url.replace("&", "&amp;")
                    return url
            
            return None
            
        except Exception as e:
            logger.error(f"Error getting BaseSpace presigned URL for {file_id}: {e}")
            return None
    
    def upload_results_to_s3(self, run_id: str, barcode: str, 
                            local_dir: Path) -> bool:
        """Upload results to S3 bucket."""
        s3_path = f"s3://{self.s3_bucket}/{self.s3_prefix}/{run_id}/{barcode}/"
        
        try:
            cmd = [
                "aws", "s3", "cp", str(local_dir), s3_path,
                "--recursive", "--quiet"
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info(f"Successfully uploaded {barcode} results to S3")
            return True
            
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to upload {barcode} to S3: {e.stderr}")
            return False
        except Exception as e:
            logger.error(f"Error uploading {barcode} to S3: {e}")
            return False