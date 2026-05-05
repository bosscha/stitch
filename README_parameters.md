# Gaia Clustering Parameters Guide

This guide explains the parameters used in the `.ext` configuration files (like `fullsky.ext`) for the Gaia clustering pipeline.

## 📂 Directories & Infrastructure
| Parameter | Description |
| :--- | :--- |
| `rootdir` | Base directory of the project (usually set via `GAIA_ROOT` env). |
| `wdir` | Working directory where processing happens and products are stored. |
| `prefile` | Prefix for the generated result files (e.g., `ocres-isofit`). |
| `isomodel` | Path to the isochrone models (e.g., MIST serial models). |
| `savedb` | Save cluster members and metadata to PostgreSQL (`yes`/`no`). |
| `dbhost`, `dbuser`, `dbpass`, `dbname`, `dbtable` | PostgreSQL connection details. |

## ⚙️ Workflow Control
| Parameter | Description |
| :--- | :--- |
| `optim` | Perform ABC/MCMC optimization for DBSCAN parameters (`yes`/`no`). |
| `algo` | Clustering algorithm selection (`dbscan` or `hdbscan`). |
| `tail` | Perform a 2nd step extraction to find isolated members (tails) (`yes`/`no`). |
| `iso` | Perform isochrone fitting to estimate age, metallicity, and mass (`yes`/`no`). |
| `zpt` | Apply Zero Point offset correction for parallax (Lindegren 2020) (`yes`/`no`). |
| `pca` | Add Principal Component (PC) components and PC vector files to output (`yes`/`no`). |

## 🔍 Data Filtering (Star Selection)
| Parameter | Description |
| :--- | :--- |
| `mindist` / `maxdist` | Range for star distances in parsecs (pc). |
| `minvra` / `maxvra` | Range for Proper Motion in RA (km/s). |
| `minvdec` / `maxvdec` | Range for Proper Motion in Dec (km/s). |
| `ming` / `maxg` | Range for Apparent Magnitude (`phot_g_mean_mag`). |

## 🧠 DBSCAN & Weighting (Used if `optim = no`)
| Parameter | Description |
| :--- | :--- |
| `eps` | The $\epsilon$ distance parameter for DBSCAN. |
| `mnei` | `min_neighbor`: Minimum stars to define a core point. |
| `mcl` | `min_cluster`: Minimum stars to form a cluster. |
| `w3d` | Weighting for 3D spatial positions. |
| `wvel` | Weighting for velocity components. |
| `whrd` | Weighting for Photometry (Hertzsprung-Russell Diagram). |

## 📈 MCMC Optimization (Used if `optim = yes`)
| Parameter | Description |
| :--- | :--- |
| `nchain` | Number of MCMC chains (iterations) to run. |
| `nburnout` | Number of initial iterations to discard (burn-in). |
| `maxiter` | Hard limit on total MCMC trials. |
| `minQc` | Minimum required contrast quality ($Q_c$). |
| `minQn` | Minimum required member count quality ($Q_n$). |
| `forcedminstars` | Minimum stars forced if the `minQn` criteria isn't reached. |
| `niterqminq` | Trials per relaxation step during initialization (Final trials = `niterqminq * 30`). |
| `mingoodsolution` | Minimum "good" solutions needed during the pre-check. |

## 🔄 Cycle & Stop Parameters
| Parameter | Description |
| :--- | :--- |
| `cyclemax` | Maximum number of search cycles. |
| `minstarstop` | Minimum stars required to continue cycling. |
| `minchainreached` | Minimum chain length required to analyze a solution. |
| `qcminstop` | $Q_c$ threshold to stop cycling after the first pass. |
| `wratiominstop` | Minimum ratio between `w3d` and `wvel` (used to verify if it's an OC). |

## 🎯 Priors & Metrics
| Parameter | Description |
| :--- | :--- |
| `epsmean` / `epsdisp` | Prior mean and dispersion for the $\epsilon$ parameter. |
| `min_nei` / `min_cl` | Prior mean for min neighbors and min cluster size. |
| `ncoredisp` | Dispersion used for both `min_nei` and `min_cl` priors. |
| `w3dmean` / `w3ddisp` | Prior for spatial weighting. |
| `wvelmean` / `wveldisp` | Prior for velocity weighting. |
| `whrdmean` / `whrddisp` | Prior for CMD weighting. |
| `clustermax` | Maximum clusters allowed; if exceeded, `find_clusters` returns 0. |
| `labels` | Method for labeling solutions (e.g., `Qc`, `Qn`, or `QcQn`). |
| `nboot` | Number of bootstrap iterations for error estimation. |
