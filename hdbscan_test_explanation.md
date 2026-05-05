# HDBSCAN Integration for Gaia Clustering
## Technical Report

### 1. Overview
This report explains the implementation and testing of **HDBSCAN** (Hierarchical Density-Based Spatial Clustering of Applications with Noise) within the `GaiaClustering` pipeline via the `extra_hdbscan.jl` script. 

While the core pipeline is built around classic DBSCAN, HDBSCAN offers significant advantages for stellar cluster extraction, particularly in its ability to find clusters of varying densities and its reduced sensitivity to the `epsilon` parameter.

### 2. Implementation Strategy: Dynamic Method Overriding
To integrate HDBSCAN without modifying the core `GaiaClustering` package (which would affect all other scripts), we utilized Julia's **dynamic method overriding**.

The script `extra_hdbscan.jl` imports the package and then redefines the `clusters` function:

```julia
import GaiaClustering: clusters

function GaiaClustering.clusters(data, epsilon, leaf, minneigh, mincluster)
    # New HDBSCAN implementation using PyCall
    ...
end
```

Because Julia uses dynamic dispatch, any subsequent calls to `clusters()`—even those inside pre-compiled functions like `find_clusters2()` or the MCMC optimization loop—are automatically routed to the new HDBSCAN implementation.

### 3. Parameter Mapping
HDBSCAN and DBSCAN share some conceptual similarities but use different parameters. The mapping used in the test script is as follows:

| DBSCAN Parameter | HDBSCAN Equivalent | Role in HDBSCAN |
| :--- | :--- | :--- |
| `min_cluster` (`mcl`) | `min_cluster_size` | The minimum number of stars to form a valid cluster. |
| `min_neighbor` (`mnei`) | `min_samples` | Controls how conservative the clustering is (higher values = more noise). |
| `epsilon` (`eps`) | `cluster_selection_epsilon` | A merge threshold. Clusters closer than this value are kept together. |

### 4. Integration with MCMC Optimization
The HDBSCAN test is fully integrated with the existing **ABC-MCMC** (Approximate Bayesian Computation Markov Chain Monte Carlo) optimizer. 

1. **Random Walk**: The optimizer generates random combinations of `w3d`, `wvel`, `whrd`, `eps`, `mcl`, and `mnei`.
2. **Evaluation**: For each step, it calls the overridden `clusters` function (now using HDBSCAN).
3. **Metric Calculation**: The quality metrics ($Q_c$, $Q_n$) are calculated on the HDBSCAN results.
4. **Convergence**: The loop finds the optimal parameters specifically for the HDBSCAN algorithm.

> [!NOTE]
> HDBSCAN requires `min_cluster_size >= 2`. The script includes a safety check to force any MCMC-generated value of `1` up to `2` to prevent Python errors.

### 5. Usage
To run the HDBSCAN version of the analysis:

```bash
julia scripts/extra_hdbscan.jl [options] votable_file.vot
```

Ensure the Python environment has the required dependencies:
```bash
pip install hdbscan scikit-learn
```
