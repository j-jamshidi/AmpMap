from flask import Flask, render_template, jsonify, send_file, abort
import os
import glob
import re
import paramiko
import shlex
from pathlib import Path
from time import sleep
from os.path import basename, relpath
from threading import Thread, Lock
import logging

# Import configuration
from config import HOSTNAME, USERNAME, BASE_PATH, LOCAL_PATH, PEM_PATH, FLASK_HOST, FLASK_PORT

class RemoteFileMonitor:
    def __init__(self, hostname=HOSTNAME, base_path=BASE_PATH, local_path=LOCAL_PATH, pem_path=PEM_PATH):
        self.hostname = hostname
        self.base_path = base_path
        self.local_path = Path(local_path)
        self.pem_path = Path(pem_path)
        self.host_ip = hostname
        self.username = "ubuntu"
        
        # Set up logging
        self.logger = logging.getLogger('RemoteMonitor')
        self.logger.setLevel(logging.INFO)
        
        if not self.pem_path.exists():
            raise FileNotFoundError(f"PEM file not found: {self.pem_path}")
            
        # Create local directory if it doesn't exist
        self.local_path.mkdir(parents=True, exist_ok=True)
        
        # Set up SSH client with proper host key policy
        self.ssh = paramiko.SSHClient()
        self.ssh.load_system_host_keys()
        if '3.24.162.88' in hostname:
            self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        
        self.sftp = None
        self.processed_signals = set()
        self._lock = Lock()

    def connect(self):
        """Establish SSH connection"""
        try:
            if not os.path.exists(self.pem_path):
                raise FileNotFoundError(f"PEM file not found at {self.pem_path}")
                
            # Close existing connections if any
            if self.sftp:
                try:
                    self.sftp.close()
                except:
                    pass
            if self.ssh:
                try:
                    self.ssh.close()
                except:
                    pass
                    
            self.ssh.connect(
                hostname=self.host_ip,
                username=USERNAME,
                key_filename=str(self.pem_path),
                timeout=10
            )
            self.sftp = self.ssh.open_sftp()
            self.logger.info(f"Connected to {self._sanitize_log_input(self.hostname)}")
        except Exception as e:
            self.logger.error(f"Connection failed: {self._sanitize_log_input(str(e))}")
            self.sftp = None
            raise

    def _sanitize_path(self, path):
        """Sanitize path to prevent traversal attacks"""
        return re.sub(r'\.\.[\/\\]', '', path)
    
    def _sanitize_log_input(self, text):
        """Sanitize input for logging to prevent log injection"""
        if isinstance(text, str):
            return text.replace('\n', ' ').replace('\r', ' ')
        return str(text)

    def check_new_signal_files(self):
        """Check for new .analysed signal files"""
        try:
            cmd = f'find {shlex.quote(self.base_path)} -maxdepth 1 -type f -name "*.analysed"'
            _, stdout, stderr = self.ssh.exec_command(cmd)
            
            new_signals = set(stdout.read().decode().splitlines())
            signals_to_process = new_signals - self.processed_signals
            
            if not signals_to_process:
                self.logger.info("No new .analysed files to process")
                return
            
            for signal_file in signals_to_process:
                folder_name = basename(signal_file).replace('.analysed', '')
                folder_name = re.sub(r'[^\w\-_.]', '', folder_name)
                self.process_data_folder(folder_name)
                with self._lock:
                    self.processed_signals.add(signal_file)
                
        except Exception as e:
            self.logger.error(f"Error checking for new signal files: {self._sanitize_log_input(str(e))}")
            raise

    def process_data_folder(self, folder_name):
        """Process a data folder when its signal file is detected"""
        try:
            safe_folder_name = shlex.quote(folder_name)
            data_folder = f"{self.base_path}/{safe_folder_name}"
            
            # Check if the data folder exists
            test_cmd = f'test -d {shlex.quote(data_folder)} && echo "exists"'
            _, stdout, stderr = self.ssh.exec_command(test_cmd)
            if not stdout.read().decode().strip():
                self.logger.error(f"Data folder {self._sanitize_log_input(data_folder)} does not exist")
                return
            
            # Download .html, .csv files from the main folder
            self.download_files_by_pattern(data_folder, ['*.html', '*.csv'])
            
            # Download .log files from barcode directories
            for barcode in range(1, 93):
                barcode_dir = f"{data_folder}/result/barcode{barcode:02d}"
                self.download_files_by_pattern(barcode_dir, ['*.log', '*.txt', '*.xml'])
            
            # Rename .analysed file to .done
            analysed_file = f"{self.base_path}/{safe_folder_name}.analysed"
            done_file = f"{self.base_path}/{safe_folder_name}.done"
            mv_cmd = f'mv {shlex.quote(analysed_file)} {shlex.quote(done_file)}'
            self.ssh.exec_command(mv_cmd)
            
        except Exception as e:
            self.logger.error(f"Error processing folder {self._sanitize_log_input(folder_name)}: {self._sanitize_log_input(str(e))}")
            raise

    def download_files_by_pattern(self, remote_dir, patterns):
        """Download files matching specific patterns from a remote directory"""
        try:
            _, stdout, stderr = self.ssh.exec_command(f'test -d {shlex.quote(remote_dir)} && echo "exists"')
            if not stdout.read().decode().strip():
                self.logger.info(f"Directory {self._sanitize_log_input(remote_dir)} does not exist, skipping")
                return
            
            safe_patterns = [shlex.quote(p) for p in patterns if re.match(r'^[\w*.-]+$', p)]
            if not safe_patterns:
                return
                
            pattern_str = ' -o '.join([f'-name {p}' for p in safe_patterns])
            cmd = f'find {shlex.quote(remote_dir)} -maxdepth 1 -type f \\( {pattern_str} \\)'
            _, stdout, stderr = self.ssh.exec_command(cmd)
            
            files = stdout.read().decode().splitlines()
            for remote_file in files:
                try:
                    if '..' in remote_file or not remote_file.startswith(self.base_path):
                        self.logger.warning(f"Skipping suspicious file path: {self._sanitize_log_input(remote_file)}")
                        continue
                        
                    relative_path = relpath(remote_file, self.base_path)
                    local_file = self.local_path / relative_path
                    local_file.parent.mkdir(parents=True, exist_ok=True)
                    
                    if self.sftp is None:
                        self.sftp = self.ssh.open_sftp()
                    self.sftp.get(remote_file, str(local_file))
                    self.logger.info(f"Successfully downloaded {self._sanitize_log_input(remote_file)}")
                    
                except Exception as e:
                    self.logger.error(f"Error downloading {self._sanitize_log_input(remote_file)}: {self._sanitize_log_input(str(e))}")
                    # Try to reconnect if SFTP is failing
                    try:
                        self.sftp = self.ssh.open_sftp()
                    except:
                        pass
                    continue
                    
        except Exception as e:
            self.logger.error(f"Error processing directory {self._sanitize_log_input(remote_dir)}: {self._sanitize_log_input(str(e))}")
            raise

    def run(self, check_interval=300):

        self.logger.info(f"Starting monitoring of {self._sanitize_log_input(self.hostname)}:{self._sanitize_log_input(self.base_path)}")
        
        connected = False
        retry_count = 0
        
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
                        
                wait_time = min(check_interval * (2 ** min(retry_count, 3)), 600)
                sleep(wait_time)

class GUI_AmpMap:
    def __init__(self):
        # Set up logging first
        self.logger = logging.getLogger('GUI_AmpMap')
        
        self.logger.info("Initializing GUI AmpMap")
        self.logger.info(f"Local path: {LOCAL_PATH}")
        self.logger.info(f"Remote host: {HOSTNAME}")
        
        # Initialize Flask app
        self.app = Flask(__name__)
        self.init_routes()
        self.logger.info("Flask app initialized")
        
        # Initialize remote monitor
        self.remote_monitor = RemoteFileMonitor()
        self.logger.info("Remote monitor initialized")

    def init_routes(self):
        """Initialize Flask routes"""
        
        @self.app.route('/')
        def index():
            return render_template('index.html')

        @self.app.route('/api/runs')
        def get_runs():
            try:
                runs = [d for d in os.listdir(str(LOCAL_PATH)) 
                        if os.path.isdir(os.path.join(str(LOCAL_PATH), d))]
                return jsonify(sorted(runs))
            except Exception as e:
                return jsonify({'error': str(e)}), 500

        @self.app.route('/api/runs/<run_id>/html')
        def get_run_html(run_id):
            try:
                run_dir = os.path.join(str(LOCAL_PATH), run_id)
                html_files = glob.glob(os.path.join(run_dir, '*.html'))
                if not html_files:
                    return jsonify({'error': 'HTML file not found'}), 404
                
                with open(html_files[0], 'r') as f:
                    content = f.read()
                return content
            except Exception as e:
                return jsonify({'error': str(e)}), 500

        @self.app.route('/api/runs/<run_id>/samples')
        def get_samples(run_id):
            try:
                samples = set()
                result_dir = os.path.join(str(LOCAL_PATH), run_id, 'result')
                
                for i in range(1, 93):
                    barcode_dir = os.path.join(result_dir, f'barcode{i:02d}')
                    if not os.path.exists(barcode_dir):
                        continue
                        
                    report_files = glob.glob(os.path.join(barcode_dir, '*_report.txt'))
                    for report_file in report_files:
                        episode = os.path.basename(report_file).replace('_report.txt', '')
                        if episode:
                            samples.add(f"{episode}_b{i:02d}")
                
                return jsonify(sorted(list(samples)))
            except Exception as e:
                return jsonify({'error': str(e)}), 500

        @self.app.route('/api/runs/<run_id>/samples/<sample_id>/report')
        def get_sample_report(run_id, sample_id):
            try:
                result_dir = os.path.join(str(LOCAL_PATH), run_id, 'result')
                
                match = re.match(r"(.+)_b(\d{2})$", sample_id)
                if not match:
                    return jsonify({'error': 'Invalid sample ID format'}), 400
                
                original_sample_id, barcode = match.groups()
                barcode_dir = os.path.join(result_dir, f'barcode{barcode}')
                report_path = os.path.join(barcode_dir, f'{original_sample_id}_report.txt')
                
                if not os.path.exists(report_path):
                    return jsonify({'error': 'Report not found'}), 404
                    
                with open(report_path, 'r') as f:
                    report_content = f.read()
                    
                return report_content
            except Exception as e:
                return jsonify({'error': str(e)}), 500

        @self.app.route('/api/runs/<run_id>/samples/<sample_id>/download/<file_type>')
        def download_file(run_id, sample_id, file_type):
            try:
                result_dir = os.path.join(str(LOCAL_PATH), run_id, 'result')
                
                match = re.match(r"(.+)_b(\d{2})$", sample_id)
                if not match:
                    return jsonify({'error': 'Invalid sample ID format'}), 400
                
                original_sample_id, barcode = match.groups()
                
                file_names = {
                    'xml': f'{original_sample_id}.xml',
                    'HapCUT2': f'HapCUT2.log',
                    'WhatsHap': f'whatshap.log',
                    'Pipeline': f'pipeline.log'
                }
                
                if file_type not in file_names:
                    return jsonify({'error': 'Invalid file type'}), 400
                    
                file_name = file_names[file_type]
                barcode_dir = os.path.join(result_dir, f'barcode{barcode}')
                file_path = os.path.join(barcode_dir, file_name)
                
                if not os.path.exists(file_path):
                    return jsonify({'error': 'File not found'}), 404
                    
                return send_file(file_path, as_attachment=True)
            except Exception as e:
                return jsonify({'error': str(e)}), 500



    def start(self):
        """Start both the monitoring thread and Flask app"""
        try:
            self.logger.info("Starting GUI AmpMap components")
            
            # Start remote monitoring in a separate thread
            monitoring_thread = Thread(target=self.remote_monitor.run)
            monitoring_thread.daemon = True
            monitoring_thread.start()
            self.logger.info("Remote monitoring thread started")
            
            # Start Flask app
            self.logger.info(f"Starting web interface on {FLASK_HOST}:{FLASK_PORT}")
            self.app.run(host=FLASK_HOST, port=FLASK_PORT)
            
        except Exception as e:
            self.logger.error(f"Error starting GUI AmpMap: {str(e)}")
            raise

if __name__ == '__main__':
    # Set up logging configuration
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    
    # Ensure log directory exists
    log_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Create logs directory
    logs_dir = os.path.join(log_dir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)
    
    # Set up common file handler for all logs
    file_handler = logging.FileHandler(os.path.join(logs_dir, 'gui_ampmap.log'))
    file_handler.setFormatter(logging.Formatter(log_format))
    
    # Configure root logger
    logging.basicConfig(
        level=logging.INFO,
        format=log_format,
        handlers=[
            logging.StreamHandler(),  # Console output
            file_handler  # File output
        ]
    )
    
    # Configure component-specific loggers
    remote_monitor_logger = logging.getLogger('RemoteMonitor')
    remote_monitor_logger.setLevel(logging.INFO)
    remote_monitor_logger.propagate = False  # Prevent propagation to root logger
    remote_monitor_handler = logging.FileHandler(os.path.join(logs_dir, 'remotemonitor.log'))
    remote_monitor_handler.setFormatter(logging.Formatter(log_format))
    remote_monitor_logger.addHandler(remote_monitor_handler)
    
    # GUI_AmpMap logger - separate from root logger
    gui_logger = logging.getLogger('GUI_AmpMap')
    gui_logger.setLevel(logging.INFO)
    gui_logger.propagate = False  # Prevent propagation to root logger
    gui_handler = logging.FileHandler(os.path.join(logs_dir, 'gui_ampmap.log'))
    gui_handler.setFormatter(logging.Formatter(log_format))
    gui_logger.addHandler(gui_handler)
    gui_logger.addHandler(logging.StreamHandler())  # Console output for GUI
    
    root_logger = logging.getLogger()
    root_logger.info("=== GUI AmpMap Starting ===")
    root_logger.info(f"Log directory: {logs_dir}")
    
    try:
        # Start the application
        app = GUI_AmpMap()
        app.start()
    except Exception as e:
        root_logger.error(f"Fatal error: {str(e)}", exc_info=True)
        raise
