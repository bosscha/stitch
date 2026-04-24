-- PostgreSQL Optimization Script for Gaia HEALPix Pipeline
-- Note: This requires PostgreSQL superuser privileges to alter system parameters.
-- ONLY RUN THIS AFTER YOU HAVE FINISHED INGESTING ALL GAIA CSVs.

-- 1. System Configuration Tuning for Data Warehousing (requires PostgreSQL restart)
-- These settings assume a modern workstation with at least 32GB RAM and SSD storage.
ALTER SYSTEM SET shared_buffers = '8GB';
ALTER SYSTEM SET work_mem = '1GB';
ALTER SYSTEM SET maintenance_work_mem = '4GB';
ALTER SYSTEM SET random_page_cost = '1.1';
ALTER SYSTEM SET max_worker_processes = '16';
ALTER SYSTEM SET max_parallel_workers_per_gather = '8';
ALTER SYSTEM SET max_parallel_maintenance_workers = '4';

-- 2. Create the Functional Index for HEALPix bit-shifting
-- This index explicitly calculates the Level 5 nested pixel for all 1.8 billion rows
CREATE INDEX IF NOT EXISTS idx_healpix_level5 ON gaia_source ((source_id >> 49));

-- 3. Cluster the Data
-- WARNING: This will physically rewrite the 1.8 billion rows on your disk.
-- It will take several hours and require free disk space equivalent to the table size.
-- However, once completed, it turns the extraction from a random-seek nightmare into a sequential lightning-fast read.
CLUSTER gaia_source USING idx_healpix_level5;

-- 4. Update Table Statistics
-- This ensures the query planner knows how the new index and clustering are distributed.
ANALYZE gaia_source;
