#!/bin/bash
set -euo pipefail

# Get script directory for relative paths
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# Configuration
BASEDIR="/EBSDataDrive/ONT/Runs"
ANALYSIS_SCRIPT="/EBSDataDrive/ONT/script/run.sh"
LOCKFILE="/tmp/ont_watchdog.lock"
LOGFILE="${SCRIPT_DIR}/watchdog.log"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >> "$LOGFILE"
}

cleanup() {
    if [ "${LOCKFILE_FD:-}" != "" ]; then
        flock -u "$LOCKFILE_FD"
    fi
    rm -f "$LOCKFILE"
    log "Watchdog stopped"
    exit 0
}

# Set up trap for cleanup
trap cleanup SIGTERM SIGINT EXIT

# Ensure exclusive execution using flock
exec {LOCKFILE_FD}>"$LOCKFILE"
if ! flock -n "$LOCKFILE_FD"; then
    echo "Another instance is running. Exiting."
    exit 1
fi

# Write PID to lockfile
echo $$ > "$LOCKFILE"

log "Watchdog started with PID $$"

while true; do
    # Use find instead of ls
    while IFS= read -r -d $'\0' complete_file; do
        runid=$(basename "$complete_file" .complete)
        log "Processing run: $runid"
        
        # Create target directory
        target_dir="$BASEDIR/$runid"
        if ! mkdir -p "$target_dir"; then
            log "ERROR: Failed to create directory $target_dir"
            continue
        fi
        
        # Check if tar.gz exists
        if [ ! -f "$BASEDIR/${runid}.tar.gz" ]; then
            log "ERROR: ${runid}.tar.gz not found"
            continue
        fi
        
        # Move to analyzing state
        if ! mv "$complete_file" "$BASEDIR/${runid}.analysing"; then
            log "ERROR: Failed to move to analyzing state for $runid"
            continue
        fi
        
        # Extract archive
        if ! tar -zxf "$BASEDIR/${runid}.tar.gz" -C "$target_dir"; then
            log "ERROR: Failed to extract ${runid}.tar.gz"
            mv "$BASEDIR/${runid}.analysing" "$complete_file"
            continue
        fi
        
        # Remove archive after successful extraction
        rm "$BASEDIR/${runid}.tar.gz"
        
        # Run analysis
        if ! bash "$ANALYSIS_SCRIPT" "$runid"; then
            log "ERROR: Analysis failed for $runid"
            mv "$BASEDIR/${runid}.analysing" "$BASEDIR/${runid}.failed"
            continue
        fi
        
        # Mark as analyzed
        mv "$BASEDIR/${runid}.analysing" "$BASEDIR/${runid}.analysed"
        log "Successfully processed $runid"
        
    done < <(find "$BASEDIR" -maxdepth 1 -name "*.complete" -print0)
    
    # Sleep for 2 hours if no files found
    sleep 7200
done