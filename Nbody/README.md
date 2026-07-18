# GAIA Hypercluster N-Body Simulations

This repository contains N-Body gravitational simulations designed to model hyperclusters across various spatial dimensions (ranging from 1D to 100D). The project features multiple computational backends and a robust PostgreSQL database integration for downstream analysis.

## Implementations

### 1. Rust CPU-Based Integrators
Highly optimized Rust implementations are provided for each dimensionality (e.g., `rust_nbody_1D`, `rust_nbody_3D`, `rust_nbody_100D`). These are compiled into executables and symlinked into this directory (e.g., `nbody-rust-1D`). 

### 2. PyTorch GPU-Based Integrators (ROCm)
Python scripts utilizing PyTorch on the AMD ROCm backend (`nbody-hypercluster-rocm-*D.py`) to compute heavy $O(N^2)$ pairwise forces in a vectorized manner. 
- Features Velocity Verlet Integration and Salpeter Initial Mass Function (IMF) sampling.
- Designed to run efficiently on AMD GPUs inside the `astro_env` Python environment.

## Running the Simulations

You can run batches of the Rust simulations using the provided bash scripts:

- **Batch 1 (1D, 2D, 3D):**
  ```bash
  ./run_simulations_batch1.sh
  ```
- **Batch 2 (4D, 5D, 6D):**
  ```bash
  ./run_simulations_batch2.sh
  ```
- **Batch 3 (25D, 50D, 100D):**
  ```bash
  ./run_simulations_batch3.sh
  ```

For the ROCm-based PyTorch simulations, execute the scripts directly via python:
```bash
conda activate astro_env
python nbody-hypercluster-rocm-3D.py
```

## Database Integration (PostgreSQL)

The simulations heavily utilize a local PostgreSQL database to store cluster states, trajectory data, and snapshots.

- **Database Name**: `hypercluster`
- **Default User**: `stephane`
- **Host**: `localhost`

### Data Ingestion
Both Rust and Python simulations log outputs to the database. The tables store coordinate spaces, velocities, dimensions, and integration steps natively. 

## Analysis Tools

The `Analysis/` subdirectory contains Python and Julia scripts to query the `hypercluster` database and produce physical plots. 

**Examples of Analysis scripts:**
- `plot_radial_distance_time_series_from_db.py`: Plots time series of Lagrangian radii/radial distances.
- `plot_virial_evolution_from_db.py` / `.jl`: Tracks the virial ratio over the simulation.
- `plot_radial_density_from_db.py`: Calculates and plots projected/radial density profiles.
- `plot_radial_mass_function_from_db.py`: Examines the evolution of mass functions dynamically.

These scripts usually accept arguments like `--dim` for filtering spatial dimensions and standard PostgreSQL connection flags (`--dbname`, `--user`, `--password`, `--host`).
