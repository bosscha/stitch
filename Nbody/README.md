# N-Body Simulation for AMD GPU (ROCm)

This is a PyTorch-based N-Body gravitational simulation designed to run efficiently on AMD GPUs (e.g., AMD 8060S) using ROCm. By leveraging PyTorch tensor operations on the `cuda` backend, it computes heavy $O(N^2)$ pairwise forces in a fully vectorized manner.

## Features
- **GPU Accelerated**: Executes primarily on the GPU to drastically reduce integration time compared to CPU equivalents.
- **Velocity Verlet Integrator**: Provides stable and energy-conserving orbital integration over long periods.
- **Salpeter Initial Mass Function (IMF)**: Natively samples realistic stellar masses on the GPU using Inverse Transform Sampling.
- **Physical Time Conversion**: Internally translates N-body time scaling into Megayears (Myr) by assuming Length=1 Parsec and actual cluster masses (in $M_\odot$).
- **Softened Gravity**: Uses a smoothing parameter to prevent numerical explosions during extremely close stellar encounters.

## Requirements & Environment

The script is intended to run inside the `astro_env` python environment. 

### Dependencies:
- `torch` (compiled with ROCm support)
- `numpy`

## Usage

First, enter the Distrobox development environment and activate the python environment:
```bash
distrobox enter dev-env
conda activate astro_env  # or source path/to/astro_env/bin/activate if using venv
```

Then, simply execute the python script:
```bash
python nbody-rocm.py
```

The terminal will output the progress every 1,000 steps, indicating the current step and the total elapsed physical time (in Myr). Once complete, the final simulation arrays (positions, velocities, and masses) are returned to the CPU and are ready for downstream analysis or plotting.
