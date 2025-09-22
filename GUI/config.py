import os
from pathlib import Path

# Configuration - Modify these paths as needed
HOSTNAME = os.getenv('HOSTNAME', '3.24.162.88')
USERNAME = os.getenv('USERNAME', 'ubuntu')
BASE_PATH = os.getenv('BASE_PATH', '/EBSDataDrive/ONT/Runs')
LOCAL_PATH = Path(os.getenv('LOCAL_PATH', '/app/data'))
PEM_PATH = Path(os.getenv('PEM_PATH', '/app/config/ssh_key.pem'))

# Network Configuration
FLASK_HOST = os.getenv('FLASK_HOST', '0.0.0.0')
FLASK_PORT = int(os.getenv('FLASK_PORT', '5001'))