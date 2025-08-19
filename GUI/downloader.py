from pathlib import Path
from time import sleep
from os.path import basename, relpath
from logging import getLogger, basicConfig, INFO, FileHandler, StreamHandler, Formatter
from subprocess import run, PIPE
from shlex import quote
import paramiko
import re
import logging

class RemoteFileMonitor:
    def __init__(self, 
                 hostname='ubuntu@3.24.162.88',
                 base_path='/EBSDataDrive/ONT/Runs',
                 local_path='/Users/javadjamshidi/Desktop/Runs',
                 pem_path='/Volumes/1-ONT_Refs/gaia-ec2-chromwell.pem'):
        """
        Initialize the Remote File Monitor
        
        Args:
            hostname (str): Remote hostname (username@ip)
            base_path (str): Base path to monitor
            local_path (str): Local path to save downloaded files
            pem_path (str): Path to .pem key file
        """
        self.hostname = hostname
        self.base_path = self._sanitize_path(base_path)
        self.local_path = Path(local_path).resolve()
        self.pem_path = Path(pem_path).resolve()
        self.processed_signals = set()
        
        # Parse hostname once
        self.username, self.host_ip = hostname.split('@')
        
        # Set up logging with instance-specific logger
        self.logger = getLogger(f"{__name__}.{id(self)}")
        if not self.logger.handlers:
            handler = FileHandler('remote_monitor.log')
            handler.setFormatter(Formatter('%(asctime)s - %(levelname)s - %(message)s'))
            self.logger.addHandler(handler)
            self.logger.addHandler(StreamHandler())
            self.logger.setLevel(INFO)
        
        # Validate paths
        if not self.local_path.parent.exists():
            raise ValueError(f"Parent directory does not exist: {self.local_path.parent}")
        if not self.pem_path.exists():
            raise ValueError(f"PEM file does not exist: {self.pem_path}")
            
        # Create local directory if it doesn't exist
        self.local_path.mkdir(parents=True, exist_ok=True)
        
        # Set up SSH client with proper host key policy
        self.ssh = paramiko.SSHClient()
        # Load system host keys for security
        self.ssh.load_system_host_keys()
        # Only add policy for known development environments
        if '3.24.162.88' in hostname:  # Your specific dev server
            self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        else:
            self.ssh.set_missing_host_key_policy(paramiko.RejectPolicy())
        
        self.sftp = None

    def _sanitize_path(self, path):
        """Sanitize path to prevent traversal attacks"""
        # Remove any path traversal attempts
        sanitized = re.sub(r'\.\.[\/\\]', '', path)
        return sanitized
    
    def _sanitize_log_input(self, text):
        """Sanitize input for logging to prevent log injection"""
        if isinstance(text, str):
            # Remove newlines and control characters
            return re.sub(r'[\r\n\x00-\x1f\x7f-\x9f]', '', text)
        return str(text)

    def connect(self):
        """Establish SSH connection"""
        try:
            self.ssh.connect(
                hostname=self.host_ip,
                username=self.username,
                key_filename=str(self.pem_path)
            )
            self.sftp = self.ssh.open_sftp()
            self.logger.info(f"Connected to {self._sanitize_log_input(self.hostname)}")
        except Exception as e:
            self.logger.error(f"Connection failed: {self._sanitize_log_input(str(e))}")
            raise

    def check_new_signal_files(self):
        """Check for new .analysed signal files"""
        try:
            # Use parameterized command to prevent injection
            cmd = f'find {quote(self.base_path)} -maxdepth 1 -type f -name "*.analysed"'
            self.logger.info(f"Executing command: {cmd}")
            _, stdout, stderr = self.ssh.exec_command(cmd)
            
            new_signals = set(stdout.read().decode().splitlines())
            self.logger.info(f"Found {len(new_signals)} .analysed files: {new_signals}")
            self.logger.info(f"Already processed: {self.processed_signals}")
            
            signals_to_process = new_signals - self.processed_signals
            self.logger.info(f"New signals to process: {signals_to_process}")
            
            if not signals_to_process:
                self.logger.info("No new .analysed files to process")
            
            for signal_file in signals_to_process:
                self.logger.info(f"Processing signal file: {signal_file}")
                folder_name = basename(signal_file).replace('.analysed', '')
                # Sanitize folder name
                folder_name = re.sub(r'[^\w\-_.]', '', folder_name)
                self.logger.info(f"Extracted folder name: {folder_name}")
                self.process_data_folder(folder_name)
                self.processed_signals.add(signal_file)
                
        except Exception as e:
            self.logger.error(f"Error checking for new signal files: {self._sanitize_log_input(str(e))}")
            raise

    def process_data_folder(self, folder_name):
        """
        Process a data folder when its signal file is detected
        
        Args:
            folder_name (str): Name of the folder to process
        """
        try:
            self.logger.info(f"Starting to process data folder: {folder_name}")
            # Sanitize and construct paths safely
            safe_folder_name = quote(folder_name)
            data_folder = f"{self.base_path}/{safe_folder_name}"
            self.logger.info(f"Data folder path: {data_folder}")
            
            # Check if the data folder exists
            test_cmd = f'test -d {quote(data_folder)} && echo "exists"'
            self.logger.info(f"Testing folder existence with: {test_cmd}")
            _, stdout, stderr = self.ssh.exec_command(test_cmd)
            result = stdout.read().decode().strip()
            self.logger.info(f"Folder test result: '{result}'")
            
            if not result:
                self.logger.error(f"Data folder {self._sanitize_log_input(data_folder)} does not exist")
                return
            
            self.logger.info(f"Data folder exists, starting downloads...")
            
            # Download .html, .csv files from the main folder
            self.download_files_by_pattern(data_folder, ['*.html', '*.csv'])
            
            # Download .log files from barcode directories
            for barcode in range(1, 93):
                barcode_dir = f"{data_folder}/result/barcode{barcode:02d}"
                self.download_files_by_pattern(barcode_dir, ['*.log', '*.txt', '*.xml'])
            
            # Rename .analysed file to .done if successful
            analysed_file = f"{self.base_path}/{safe_folder_name}.analysed"
            done_file = f"{self.base_path}/{safe_folder_name}.done"
            mv_cmd = f'mv {quote(analysed_file)} {quote(done_file)}'
            self.logger.info(f"Renaming with command: {mv_cmd}")
            self.ssh.exec_command(mv_cmd)
            self.logger.info(f"Renamed {self._sanitize_log_input(analysed_file)} to {self._sanitize_log_input(done_file)}")
            
        except Exception as e:
            self.logger.error(f"Error processing folder {self._sanitize_log_input(folder_name)}: {self._sanitize_log_input(str(e))}")
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
            _, stdout, stderr = self.ssh.exec_command(f'test -d {quote(remote_dir)} && echo "exists"')
            if not stdout.read().decode().strip():
                self.logger.info(f"Directory {self._sanitize_log_input(remote_dir)} does not exist, skipping")
                return
            
            # Sanitize patterns and find files
            safe_patterns = [quote(p) for p in patterns if re.match(r'^[\w*.-]+$', p)]
            if not safe_patterns:
                return
                
            pattern_str = ' -o '.join([f'-name {p}' for p in safe_patterns])
            cmd = f'find {quote(remote_dir)} -maxdepth 1 -type f \\( {pattern_str} \\)'
            _, stdout, stderr = self.ssh.exec_command(cmd)
            
            files = stdout.read().decode().splitlines()
            for remote_file in files:
                try:
                    # Validate file path
                    if '..' in remote_file or not remote_file.startswith(self.base_path):
                        self.logger.warning(f"Skipping suspicious file path: {self._sanitize_log_input(remote_file)}")
                        continue
                        
                    # Create local directory structure
                    relative_path = relpath(remote_file, self.base_path)
                    local_file = self.local_path / relative_path
                    local_file.parent.mkdir(parents=True, exist_ok=True)
                    
                    # Download file using SFTP (more efficient than scp subprocess)
                    self.logger.info(f"Downloading {self._sanitize_log_input(remote_file)}")
                    self.sftp.get(remote_file, str(local_file))
                    self.logger.info(f"Successfully downloaded {self._sanitize_log_input(remote_file)}")
                        
                except Exception as e:
                    self.logger.error(f"Error downloading {self._sanitize_log_input(remote_file)}: {self._sanitize_log_input(str(e))}")
                    continue
                    
        except Exception as e:
            self.logger.error(f"Error processing directory {self._sanitize_log_input(remote_dir)}: {self._sanitize_log_input(str(e))}")
            raise

    def run(self, check_interval=60):
        """
        Start monitoring for new files
        
        Args:
            check_interval (int): Time in seconds between checks
        """
        self.logger.info(f"Starting monitoring of {self._sanitize_log_input(self.hostname)}:{self._sanitize_log_input(self.base_path)}")
        self.logger.info(f"Downloads will be saved to {self._sanitize_log_input(str(self.local_path))}")
        
        # Maintain persistent connection for better performance
        connected = False
        retry_count = 0
        max_retries = 3
        
        while True:
            try:
                if not connected:
                    self.connect()
                    connected = True
                    retry_count = 0
                    
                self.check_new_signal_files()
                sleep(check_interval)
                
            except KeyboardInterrupt:
                self.logger.info("Monitoring stopped by user")
                break
            except Exception as e:
                self.logger.error(f"Unexpected error: {self._sanitize_log_input(str(e))}")
                connected = False
                retry_count += 1
                
                if self.ssh:
                    try:
                        self.ssh.close()
                    except:
                        pass
                        
                if self.sftp:
                    try:
                        self.sftp.close()
                    except:
                        pass
                        
                # Exponential backoff for retries
                wait_time = min(check_interval * (2 ** min(retry_count, 3)), 300)
                sleep(wait_time)
                
        # Cleanup on exit
        try:
            if self.sftp:
                self.sftp.close()
            if self.ssh:
                self.ssh.close()
        except:
            pass

if __name__ == "__main__":
    monitor = RemoteFileMonitor(
        hostname='ubuntu@3.24.162.88',
        base_path='/EBSDataDrive/ONT/Runs',
        local_path='/Users/javadjamshidi/Desktop/Runs',
        pem_path='/Volumes/1-ONT_Refs/gaia-ec2-chromwell.pem'
    )
    monitor.run()