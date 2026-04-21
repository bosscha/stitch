#!/bin/bash
# Auto-resuming loop for Gaia DR3 full sky ingestion
# This script will continuously restart gaia_bulk.py if it exits with an error 
# (e.g., due to a network drop).
# If gaia_bulk.py completes successfully (exit code 0), the loop will break.

# Get the directory where this bash script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/gaia_bulk.py"

echo "Starting auto-resuming Gaia ingestion loop..."

while true; do
    echo "=================================================="
    echo "Running: python3 $PYTHON_SCRIPT"
    echo "=================================================="
    
    python3 "$PYTHON_SCRIPT"
    
    # Store the exit code of the python script
    EXIT_CODE=$?
    
    if [ $EXIT_CODE -eq 0 ]; then
        echo "=================================================="
        echo "Ingestion completed successfully!"
        echo "=================================================="
        break
    else
        echo "=================================================="
        echo "Ingestion stopped with an error (exit code $EXIT_CODE)."
        echo "Likely a network failure or database timeout."
        echo "Restarting in 10 seconds... (Press Ctrl+C to abort)"
        echo "=================================================="
        sleep 10
    fi
done
