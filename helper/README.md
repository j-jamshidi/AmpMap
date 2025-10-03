# AmpMap Watchdog

Automated file processing daemon for GridION uploads and AmpMap analysis.

## Overview

The watchdog monitors for new sequencing runs uploaded from GridION devices and automatically processes them through the AmpMap pipeline.

## Workflow
```mermaid
flowchart TD
       subgraph GridION Device
              A[Upload sequencing data]
              B[Creates: run.tar.gz & run.complete]
       end

       subgraph Server
              C[Watchdog Monitoring]
              D[Detects .complete file]
              E[Extract archive & rename to .analysing]
              F[AmpMap Analysis Pipeline]
              G{Analysis Result}
              H[Rename to .analysed]
              I[Rename to .failed]
       end

       subgraph Local GUI
              J[Detects .analysed file]
              K[Downloads results]
              L[Rename to .done]
              M[Results available in GUI]
       end

       A --> B
       B --> C
       C --> D
       D --> E
       E --> F
       F --> G
       G -- Success --> H
       G -- Failure --> I
       H --> J
       J --> K
       K --> L
       L --> M
```

### Process Steps

1. **Upload Detection**: Monitors `/EBSDataDrive/ONT/Runs/` for `.complete` files
2. **File Processing**: 
   - Extracts corresponding `.tar.gz` archives
   - Changes status to `.analysing`
3. **Analysis**: Runs AmpMap pipeline automatically
4. **Status Updates**:
   - Success: `.analysing` → `.analysed` 
   - Failure: `.analysing` → `.failed`
5. **Download Ready**: GUI downloader detects `.analysed` files and changes to `.done` after download

## File States

- `.complete` - New upload ready for processing
- `.analysing` - Currently being processed
- `.analysed` - Analysis complete, ready for download
- `.failed` - Processing failed
- `.done` - Downloaded by GUI

## Usage

```bash
# Start watchdog
./ampmap_watchdog

# Check status
tail -f helper/watchdog.log
```

## Configuration

- **BASEDIR**: `/EBSDataDrive/ONT/Runs` - Monitor directory
- **ANALYSIS_SCRIPT**: `/EBSDataDrive/ONT/AmpMap/ampmap` - Pipeline script
- **Check Interval**: 30 minutes

## Requirements

- Docker (for AmpMap pipeline)
- Write permissions to BASEDIR
- AmpMap installation

## Troubleshooting
For issues, check `helper/watchdog.log` for detailed logs.