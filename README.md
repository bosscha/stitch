# gaia-stitch

Pipeline to extract stellar clusters from GAIA data as votable. The method is using a HDBSCAN-optimized clustering method in the 8-fold dimension with photometry and astrometry data. The optimization the clustering parameters is performed with an Approximate Bayesian Computation method (Weyant et al. 2013).
The extraction process is currently with a stable version using the Julia language for speed constraint. An (experimental) capability is offered from version 1.7 to fit an isochrone and to extract more isolated stars in a second step ("tails", subclumps, etc).
Be aware that some stellar clusters are not physically real and a final cleaning can be done with some Jupyter NB but are still experimental. A visual inspection of the final stellar cluster is still recommended.
The current version is working with Gaia DR3 data.

A postgresql backend is added to deal locally with Gaia sources for the full sky analysis 

**Documentation and information** can be found in the [wiki pages](https://github.com/bosscha/gaia-shock/wiki).

## Running Full Sky Analysis (`build_hp.jl`)

To run the full-sky clustering extraction using HEALPix indexing, use the `build_hp.jl` script. It takes a configuration file (in TOML format) as an argument.

### Example Usage
```bash
julia scripts/fullSky/build_hp.jl fullsky.bld
```

### Configuration
The pipeline relies on parameter files (`.ext`) referenced inside your TOML configuration file (e.g., `template/optimal_galaxy.ext`). 

#### Batch Processing Mode
You can speed up the extraction process by analyzing multiple HEALPix pixels simultaneously. To do this, edit your `.ext` file and set the `nbatch` parameter to your desired concurrency level:
```text
nbatch = 4      ## number of simultaneous batches (pixels) to process
```

#### Lock File Management
The script utilizes an `active_pixels.lock` file in the working directory to track which pixels are currently being processed across all processes. This robust PID-based locking mechanism safely allows you to run multiple instances of the script simultaneously, and ensures no duplicate work is performed. 

If you interrupt the script via `Ctrl-C`, your terminal might suspend the background processes rather than killing them. To completely terminate all running or suspended instances, use the provided helper script:
```bash
./kill_build_hp.sh
```
The pipeline will automatically identify dead PIDs and clean up any stale locks the next time you start it.
