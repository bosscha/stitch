# Gaia Full Sky Processing

This directory contains the scripts required to ingest the Gaia DR3 full sky dataset into a local PostgreSQL database and execute distributed DBSCAN clustering utilizing HEALPix parallelization.

## Prerequisites

Ensure you have your PostgreSQL database initialized and the dataset fully ingested before running the clustering algorithm.

1. **Database Creation**: Run `create_db.sql` in PostgreSQL to create the `gaiadb` schema.
2. **Ingestion Pipeline**: Run the auto-resuming ingestion wrapper `./run_gaia_bulk.sh` to download and stream all Gaia DR3 CSV pieces into your local database.

## Running the Full Sky Extraction

You can run the extraction pipeline natively via Julia.

To make running the pipeline convenient from any directory, you can add the following alias to your `~/.bashrc`:

```bash
alias run_fullsky='julia /home/stephane/Science_0/GAIA/run/scripts/fullSky/build_hp.jl'
```

After adding it, simply run:
```bash
source ~/.bashrc
```

### Execution
Once the alias is set up, you can execute the full sky extraction by simply providing the `.bld` (TOML configuration) file:

```bash
run_fullsky /home/stephane/Science_0/GAIA/stitch/template/fullsky.bld
```

*(Note: The clustering process is fully resumable. Progress is tracked automatically using a `_done_hp.csv` file created in your working directory. If your machine crashes or shuts down, simply run the command again and it will resume from the exact HEALPix pixel it left off at.)*
