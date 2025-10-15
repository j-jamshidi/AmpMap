# AmpMap GUI

![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![Docker](https://img.shields.io/badge/docker-%230db7ed.svg?style=for-the-badge&logo=docker&logoColor=white)
![Flask](https://img.shields.io/badge/flask-%23000.svg?style=for-the-badge&logo=flask&logoColor=white)

A containerized web interface for accessing and visualizing ONT amplicon analysis results.

## Features

- **Automated Data Monitoring**: Continuously monitors remote server for new analysis results
- **Web Interface**: Clean, responsive web UI for browsing runs and samples
- **Report Visualization**: Interactive display of analysis reports with tabbed sections
- **File Downloads**: Direct download of XML files and log files
- **Cross-Platform**: Works on Linux and macOS via Docker

## Quick Start

### Prerequisites

- Docker and Docker Compose installed
- SSH key file for remote server access

### 1. Pull the Docker Image

```bash
docker pull javadj/ampmap-gui:latest
```

### 2. Setup Configuration

Create a directory for your configuration:

```bash
mkdir -p ampmap-gui/config
mkdir -p ampmap-gui/data
mkdir -p ampmap-gui/logs
```

Copy your SSH key to the config directory:

```bash
cp /path/to/your/ssh_key.pem ampmap-gui/config/ssh_key.pem
chmod 600 ampmap-gui/config/ssh_key.pem
```

### 3. Create Environment File

Create `ampmap-gui/.env`:

```bash
HOSTNAME=your.server.ip
USERNAME=ubuntu
BASE_PATH=/EBSDataDrive/ONT/Runs
```

### 4. Run with Docker Compose

Create `ampmap-gui/docker-compose.yml`:

```yaml
version: '3.8'

services:
  ampmap-gui:
    image: javadj/ampmap-gui:latest
    ports:
      - "5001:5001"
    volumes:
      - ./data:/app/data
      - ./logs:/app/logs
      - ./config:/app/config
    env_file:
      - .env
    environment:
      - FLASK_ENV=production
      - FLASK_HOST=0.0.0.0
      - FLASK_PORT=5001
    restart: unless-stopped
```

Start the application:

```bash
cd ampmap-gui
docker-compose up -d
```

### 5. Access the Interface

Open your browser and navigate to:
```
http://localhost:5001
```

## Alternative: Run with Docker

```bash
docker run -d \
  --name ampmap-gui \
  --network host \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/logs:/app/logs \
  -v $(pwd)/config:/app/config \
  -e HOSTNAME=your.server.ip \
  -e USERNAME=ubuntu \
  -e BASE_PATH=/EBSDataDrive/ONT/Runs \
  javadj/ampmap-gui:latest
```

## Configuration

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `HOSTNAME` | Remote server IP/hostname | `3.24.162.88` |
| `USERNAME` | SSH username | `ubuntu` |
| `BASE_PATH` | Remote path to ONT runs | `/EBSDataDrive/ONT/Runs` |
| `LOCAL_PATH` | Local data storage path | `/app/data` |
| `PEM_PATH` | Path to SSH private key | `/app/config/ssh_key.pem` |
| `FLASK_HOST` | Flask bind address | `0.0.0.0` |
| `FLASK_PORT` | Flask port number | `5001` |

### Directory Structure

```
ampmap-gui/
├── config/
│   └── ssh_key.pem          # SSH private key
├── data/                    # Downloaded analysis results
├── logs/                    # Application logs
├── .env                     # Environment configuration
└── docker-compose.yml       # Docker Compose configuration
```

## How It Works

1. **Monitoring**: The application continuously monitors the remote server for `.analysed` files. These files are generated when a new run uploaded from GridION (indicated by a `.complete` file) is detected and analysed. For a detailed explanation of the monitoring process, see [ampmap_watchdog README](../helper/ampmap_watchdog/README.md).
2. **Download**: When new results are detected, it downloads HTML reports, CSV files, and log files
3. **Web Interface**: Provides a clean web interface to browse runs and samples
4. **Visualization**: Displays analysis reports with interactive tabs and formatting

## Network Access

**✅ Network access is enabled by default!** Other computers can access the GUI immediately.

### Find Your IP Address
```bash
# macOS/Linux
ifconfig | grep "inet " | grep -v 127.0.0.1
```

### Access from Other Computers
```
http://YOUR_COMPUTER_IP:5001
```
Example: `http://192.168.1.100:5001`

### Security Note
The GUI is accessible from any computer on your network. Only use on trusted networks.

## Troubleshooting

### Check Container Status
```bash
docker-compose ps
docker-compose logs ampmap-gui
```

### Test Access
```bash
# Local access
curl http://localhost:5001

# Network access
curl http://YOUR_IP:5001
```

### SSH Connection Issues
- Ensure SSH key has correct permissions: `chmod 600 config/ssh_key.pem`
- Verify server hostname/IP is accessible
- Check username is correct

## Development

To build from source:

```bash
git clone https://github.com/j-jamshidi/AmpMap.git
cd AmpMap/GUI
docker build -t ampmap-gui .
```

## License

MIT License

## Author
Javad Jamshidi j.jamshidi@neura.edu.au
