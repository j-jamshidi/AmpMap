# ONT AmpMap GUI

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
docker pull javadj/ontampip_gui:latest
```

### 2. Setup Configuration

Create a directory for your configuration:

```bash
mkdir -p ontampip-gui/config
mkdir -p ontampip-gui/data
mkdir -p ontampip-gui/logs
```

Copy your SSH key to the config directory:

```bash
cp /path/to/your/ssh_key.pem ontampip-gui/config/ssh_key.pem
chmod 600 ontampip-gui/config/ssh_key.pem
```

### 3. Create Environment File

Create `ontampip-gui/.env`:

```bash
# Server Configuration
HOSTNAME=your.server.ip
USERNAME=ubuntu
BASE_PATH=/EBSDataDrive/ONT/Runs

# Local Configuration (use defaults for Docker)
LOCAL_PATH=/app/data
PEM_PATH=/app/config/ssh_key.pem
```

### 4. Run with Docker Compose

Create `ontampip-gui/docker-compose.yml`:

```yaml
services:
  ontampip-gui:
    image: javadj/ontampip_gui:latest
    ports:
      - "5001:5001"
    volumes:
      - ./data:/app/data
      - ./logs:/app/logs
      - ./config:/app/config
    env_file:
      - .env
    restart: unless-stopped
```

Start the application:

```bash
cd ontampip-gui
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
  --name ontampip-gui \
  -p 5001:5001 \
  -v $(pwd)/data:/app/data \
  -v $(pwd)/logs:/app/logs \
  -v $(pwd)/config:/app/config \
  -e HOSTNAME=your.server.ip \
  -e USERNAME=ubuntu \
  -e BASE_PATH=/EBSDataDrive/ONT/Runs \
  javadj/ontampip_gui:latest
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
ontampip-gui/
├── config/
│   └── ssh_key.pem          # SSH private key
├── data/                    # Downloaded analysis results
├── logs/                    # Application logs
├── .env                     # Environment configuration
└── docker-compose.yml       # Docker Compose configuration
```

## How It Works

1. **Monitoring**: The application continuously monitors the remote server for `.analysed` files
2. **Download**: When new results are detected, it downloads HTML reports, CSV files, and log files
3. **Web Interface**: Provides a clean web interface to browse runs and samples
4. **Visualization**: Displays analysis reports with interactive tabs and formatting

## Network Access

### Access from Other Computers

To allow other computers on your network to access the GUI:

1. **Set Flask host to bind to all interfaces** (default in Docker):
   ```bash
   FLASK_HOST=0.0.0.0
   ```

2. **Find your computer's IP address**:
   ```bash
   # On Linux/macOS
   ip addr show | grep inet
   # or
   ifconfig | grep inet
   ```

3. **Access from other computers**:
   ```
   http://YOUR_COMPUTER_IP:5001
   ```
   Example: `http://192.168.1.100:5001`

4. **Firewall considerations**:
   - Ensure port 5001 is open in your firewall
   - On macOS: System Preferences → Security & Privacy → Firewall
   - On Linux: `sudo ufw allow 5001`

### Security Note

When `FLASK_HOST=0.0.0.0`, the GUI will be accessible from any computer on your network. For security:
- Only use this on trusted networks
- Consider using a VPN for remote access
- Set `FLASK_HOST=127.0.0.1` to restrict to localhost only

## Troubleshooting

### Check Logs

```bash
docker-compose logs -f ontampip-gui
```

### SSH Connection Issues

- Ensure SSH key has correct permissions (600)
- Verify server hostname/IP is accessible
- Check username is correct

### Port Conflicts

If port 5001 is in use, change it in docker-compose.yml:

```yaml
ports:
  - "8080:5001"  # Use port 8080 instead
```

## Development

To build from source:

```bash
git clone <repository>
cd ONT_amplicon_phase/GUI
docker build -t ontampip-gui .
```

## License

MIT License