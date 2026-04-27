# Gaia Clustering Parameters Guide

This document explains the parameters used in the `.ext` configuration files for the Gaia clustering pipeline.

## 📂 Directories & Infrastructure
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `rootdir` | `$GAIA_ROOT` | Base directory of the project. |
| `wdir` | `.` | Working directory for products. |
| `votdir` | `.` | Directory to store/read VOTables. |
| `plotdir` | `.` | Directory for generated plots. |
| `ocdir` | `.` | Directory for Open Cluster (OC) result files. |
| `savedb` | `no` | Whether to save results to a PostgreSQL database (`yes` or `no`). |

## ⚙️ Workflow Control
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `optim` | `no` | If `yes`, performs ABC/MCMC optimization for DBSCAN parameters. |
| `tail` | `no` | If `yes`, performs a 2nd step extraction for isolated members. |
| `iso` | `no` | If `yes`, performs isochrone fitting on detected clusters. |
| `zpt` | `no` | Apply Zero Point Correction (Lindegren 2020). |
| `pca` | `no` | Output PCA components in the OC result file. |

## 🔍 Data Filtering (Star Selection)
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `maxdist` | `1e9` | Maximum distance in parsecs. |
| `mindist` | `0.0` | Minimum distance in parsecs. |
| `minvra` / `maxvra` | `+/- 250` | Range for Proper Motion in RA (km/s). |
| `minvdec` / `maxvdec` | `+/- 250` | Range for Proper Motion in Dec (km/s). |
| `ming` / `maxg` | `0.0` / `22.0` | Range for `phot_g_mean_mag`. |

## 🧠 DBSCAN Parameters (Used if `optim = no`)
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `eps` | `2.1` | The $\epsilon$ distance parameter for DBSCAN. |
| `mnei` | `7` | `min_neighbor`: Minimum stars to define a core point. |
| `mcl` | `18` | `min_cluster`: Minimum stars to form a cluster. |
| `w3d` | `7.0` | Weighting for 3D spatial positions. |
| `wvel` | `8.0` | Weighting for velocity components. |
| `whrd` | `2.5` | Weighting for Photometry (Hertzsprung-Russell Diagram). |

## 📈 MCMC Optimization (Used if `optim = yes`)
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `nchain` | `5000` | Number of MCMC chains (iterations) to run. |
| `nburnout` | `500` | Number of initial iterations to discard (burn-in). |
| `maxiter` | `50000` | Hard limit on MCMC trials. |
| `minQc` | `2.6` | Minimum required quality for the contrast metric ($Q_c$). |
| `minQn` | `40` | Minimum required quality for the member count metric ($Q_n$). |
| `niterqminq` | `500` | Iterations per "relaxation step" during initialization. |

> [!NOTE]
> **Initialization vs Optimization**: During the pre-check (initialization), the code will try up to `niterqminq * 30` trials to find a valid starting point. This is why you might see a "maxiter" of 30,000 if `niterqminq` is 1,000, even if your main `maxiter` is set lower.

## ☄️ Tail Extraction (2nd Step)
| Parameter | Default | Description |
| :--- | :--- | :--- |
| `maxRadTail` | `250.0` | Radius cut (pc) from the cluster center. |
| `maxVelTail` | `5.0` | Velocity cut (km/s) from the cluster mean velocity. |
| `maxDistCmdTail` | `0.05` | Maximum distance to the cluster CMD sequence. |
