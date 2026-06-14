#!/bin/bash

# Script to clear the Gaia clustering results and progress from PostgreSQL
# Tables: clusters, clusters_metadata, clusters_processed_pixels

DB_NAME="${DB_NAME:-gaiadb}"
DB_USER="${DB_USER:-stephane}"
DB_HOST="${DB_HOST:-127.0.0.1}"  # Force TCP connection to avoid Peer authentication issues

echo "⚠️  WARNING: This will permanently delete all extracted clusters, metadata, and pixel progress tracking from the database."
read -p "Are you sure you want to proceed? (y/N): " confirm

if [[ "$confirm" == [yY] || "$confirm" == [yY][eE][sS] ]]; then
    echo "Clearing tables in database '$DB_NAME' on $DB_HOST..."
    
    # Set the password for psql
    export PGPASSWORD="${DB_PASS:-tallis}"
    psql -h $DB_HOST -U $DB_USER -d $DB_NAME -c "TRUNCATE TABLE clusters; TRUNCATE TABLE clusters_metadata; TRUNCATE TABLE clusters_processed_pixels;"
    
    if [ $? -eq 0 ]; then
        echo "✅ Tables successfully cleared."
    else
        echo "❌ Error: Failed to clear tables. Make sure the database is running and you have the correct password."
    fi
else
    echo "Operation cancelled. No data was deleted."
fi
