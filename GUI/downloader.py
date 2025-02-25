import paramiko
import time
import os
from pathlib import Path
import logging
import subprocess
from datetime import datetime

class RemoteFileMonitor:
    def __init__(self, 
                 hostname='ubuntu@3.24.162.88',
                 base_path='/EBSDataDrive/ONT/Runs',
                 local_path='/Users/javadjamshidi/Desktop/Runs',
                 pem_path='/Volumes/NanoDisk2/gaia-ec2-chromwell.pem'):
        """
        Initialize the Remote File Monitor
        
        Args:
            hostname (str): Remote hostname (username@ip)
            base_path (str): Base path to monitor
            local_path (str): Local path to save downloaded files
            pem_path (str): Path to .pem key file
        """
        self.hostname = hostname
        self.base_path = base_path
        self.local_path = Path(local_path)
        self.pem_path = Path(pem_path)
        self.processed_signals = set()
        
        # Set up logging
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler('remote_monitor.log'),
                logging.StreamHandler()
            ]
        )
        self.logger = logging.getLogger(__name__)
        
        # Create local directory if it doesn't exist
        self.local_path.mkdir(parents=True, exist_ok=True)
        
        # Set up SSH client
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    def connect(self):
        """Establish SSH connection"""
        try:
            self.ssh.connect(
                hostname=self.hostname.split('@')[1],
                username=self.hostname.split('@')[0],
                key_filename=str(self.pem_path)
            )
            self.logger.info(f"Connected to {self.hostname}")
        except Exception as e:
            self.logger.error(f"Connection failed: {str(e)}")
            raise

    def check_new_signal_files(self):
        """Check for new .analysed signal files"""
        try:
            cmd = f'find {self.base_path} -maxdepth 1 -type f -name "*.analysed"'
            stdin, stdout, stderr = self.ssh.exec_command(cmd)
            
            new_signals = set(stdout.read().decode().splitlines())
            signals_to_process = new_signals - self.processed_signals
            
            for signal_file in signals_to_process:
                folder_name = os.path.basename(signal_file).replace('.analysed', '')
                self.process_data_folder(folder_name)
                self.processed_signals.add(signal_file)
                
        except Exception as e:
            self.logger.error(f"Error checking for new signal files: {str(e)}")
            raise

    def process_data_folder(self, folder_name):
        """
        Process a data folder when its signal file is detected
        
        Args:
            folder_name (str): Name of the folder to processcd 
        """
        try:
            # Path to the data folder
            data_folder = f"{self.base_path}/{folder_name}"
            
            # Check if the data folder exists
            stdin, stdout, stderr = self.ssh.exec_command(f'test -d "{data_folder}" && echo "exists"')
            if not stdout.read().decode().strip():
                self.logger.error(f"Data folder {data_folder} does not exist")
                return
            
            # Download .html, .csv and .txt files from the main folder
            self.download_files_by_pattern(data_folder, ['*.html', '*.csv'])
            
            # Download .log files from barcode directories (not subdirectories)
            for barcode in range(1, 93):
                barcode_dir = f"{data_folder}/result/barcode{barcode:02d}"
                self.download_files_by_pattern(barcode_dir, ['*.log', '*.txt', '*.xml'])
            
            # Rename .analysed file to .done if successful
            analysed_file = f"{self.base_path}/{folder_name}.analysed"
            done_file = f"{self.base_path}/{folder_name}.done"
            self.ssh.exec_command(f'mv "{analysed_file}" "{done_file}"')
            self.logger.info(f"Renamed {analysed_file} to {done_file}")
            
        except Exception as e:
            self.logger.error(f"Error processing folder {folder_name}: {str(e)}")
            raise

    def download_files_by_pattern(self, remote_dir, patterns):
        """
        Download files matching specific patterns from a remote directory
        
        Args:
            remote_dir (str): Remote directory path
            patterns (list): List of file patterns to download
        """
        try:
            # Check if remote directory exists
            stdin, stdout, stderr = self.ssh.exec_command(f'test -d "{remote_dir}" && echo "exists"')
            if not stdout.read().decode().strip():
                self.logger.info(f"Directory {remote_dir} does not exist, skipping")
                return
            
            # Find files matching patterns with maxdepth 1
            pattern_str = ' -o '.join([f'-name "{p}"' for p in patterns])
            cmd = f'find "{remote_dir}" -maxdepth 1 -type f \( {pattern_str} \)'
            stdin, stdout, stderr = self.ssh.exec_command(cmd)
            
            files = stdout.read().decode().splitlines()
            for remote_file in files:
                try:
                    # Create local directory structure
                    relative_path = os.path.relpath(remote_file, self.base_path)
                    local_file = self.local_path / relative_path
                    local_file.parent.mkdir(parents=True, exist_ok=True)
                    
                    # Download file using scp
                    scp_command = [
                        'scp',
                        '-i', str(self.pem_path),
                        f'{self.hostname}:{remote_file}',
                        str(local_file)
                    ]
                    
                    self.logger.info(f"Downloading {remote_file}")
                    result = subprocess.run(scp_command, capture_output=True, text=True)
                    
                    if result.returncode == 0:
                        self.logger.info(f"Successfully downloaded {remote_file}")
                    else:
                        self.logger.error(f"Download failed: {result.stderr}")
                        
                except Exception as e:
                    self.logger.error(f"Error downloading {remote_file}: {str(e)}")
                    continue
                    
        except Exception as e:
            self.logger.error(f"Error processing directory {remote_dir}: {str(e)}")
            raise

    def run(self, check_interval=60):
        """
        Start monitoring for new files
        
        Args:
            check_interval (int): Time in seconds between checks
        """
        self.logger.info(f"Starting monitoring of {self.hostname}:{self.base_path}")
        self.logger.info(f"Downloads will be saved to {self.local_path}")
        
        while True:
            try:
                self.connect()
                self.check_new_signal_files()
                self.ssh.close()
                time.sleep(check_interval)
            except KeyboardInterrupt:
                self.logger.info("Monitoring stopped by user")
                break
            except Exception as e:
                self.logger.error(f"Unexpected error: {str(e)}")
                time.sleep(check_interval)  # Continue monitoring despite errors

if __name__ == "__main__":
    monitor = RemoteFileMonitor(
        hostname='ubuntu@3.24.162.88',
        base_path='/EBSDataDrive/ONT/Runs',
        local_path='/Users/javadjamshidi/Desktop/Runs',
        pem_path='/Volumes/NanoDisk2/gaia-ec2-chromwell.pem'
    )
    monitor.run()