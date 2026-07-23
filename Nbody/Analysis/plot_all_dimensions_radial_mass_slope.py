#!/usr/bin/env python3
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def compute_mass_slope(masses, m_min=0.1, m_max=50.0, num_bins=12):
    """
    Computes the Salpeter mass function slope alpha (dN/dM propto M^-alpha)
    for a given set of masses using log-spaced bins.
    Returns (alpha, alpha_err).
    """
    if len(masses) < 20:
        return np.nan, np.nan
        
    bins = np.logspace(np.log10(m_min), np.log10(m_max), num_bins)
    counts, _ = np.histogram(masses, bins=bins)
    
    dn_dm = counts / (bins[1:] - bins[:-1])
    bin_centers = np.sqrt(bins[1:] * bins[:-1])
    
    valid = dn_dm > 0
    if np.sum(valid) < 3:
        return np.nan, np.nan
        
    try:
        slope, cov = np.polyfit(np.log10(bin_centers[valid]), np.log10(dn_dm[valid]), 1, cov=True)
        slope_err = np.sqrt(cov[0, 0])
        fit_slope = slope[0]
    except Exception:
        try:
            fit_slope = np.polyfit(np.log10(bin_centers[valid]), np.log10(dn_dm[valid]), 1)[0]
            slope_err = np.nan
        except Exception:
            return np.nan, np.nan
            
    return -fit_slope, slope_err

def compute_radial_slope_profile(dists, masses, num_radial_bins=8, m_min=0.1, m_max=50.0):
    """
    Splits particles into quantile radial bins and fits the mass function slope alpha in each bin.
    """
    quantiles = np.linspace(0, 1, num_radial_bins + 1)
    bin_edges = np.quantile(dists, quantiles)
    
    bin_radii = []
    slopes = []
    slope_errs = []
    
    for i in range(num_radial_bins):
        lower, upper = bin_edges[i], bin_edges[i+1]
        mask = (dists >= lower) & (dists <= upper)
        
        bin_masses = masses[mask]
        bin_dists = dists[mask]
        
        if len(bin_masses) >= 20:
            bin_radii.append(np.median(bin_dists))
            alpha, alpha_err = compute_mass_slope(bin_masses, m_min, m_max)
            slopes.append(alpha)
            slope_errs.append(alpha_err)
            
    return np.array(bin_radii), np.array(slopes), np.array(slope_errs)

def main():
    parser = argparse.ArgumentParser(description="Plot and compare the final radial mass function slope alpha(R) across all dimensions.")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--m-min", type=float, default=0.1, help="Minimum mass in solar masses")
    parser.add_argument("--m-max", type=float, default=50.0, help="Maximum mass in solar masses")
    parser.add_argument("--num-radial-bins", type=int, default=8, help="Number of radial quantile bins")
    parser.add_argument("--output-grid", type=str, default="radial_mass_slope_all_dims.png", help="Multi-panel plot output filename")
    parser.add_argument("--output-overlay", type=str, default="radial_mass_slope_overlay_all_dims.png", help="Overlay plot output filename")
    args = parser.parse_args()

    print(f"Connecting to database '{args.dbname}' on '{args.host}'...")
    try:
        conn = psycopg2.connect(
            dbname=args.dbname,
            user=args.user,
            password=args.password,
            host=args.host
        )
        cursor = conn.cursor()
    except Exception as e:
        print(f"Failed to connect to database: {e}")
        return

    # Query distinct spatial dimensions
    cursor.execute("SELECT DISTINCT dim_space FROM star_snapshots ORDER BY dim_space ASC;")
    dims = [row[0] for row in cursor.fetchall()]

    if not dims:
        print("No simulation data found in table 'star_snapshots'.")
        cursor.close()
        conn.close()
        return

    print(f"Found {len(dims)} spatial dimensions to process: {dims}")

    dim_data = {}

    for dim in dims:
        print(f"\nProcessing dim_space = {dim}D...")
        cursor.execute("SELECT MIN(snapshot_id), MAX(snapshot_id) FROM star_snapshots WHERE dim_space = %s;", (dim,))
        min_snap, max_snap = cursor.fetchone()

        if min_snap is None or max_snap is None:
            print(f"  Skipping {dim}D: No snapshot range found.")
            continue

        # Load initial positions and masses
        cursor.execute("""
            SELECT position, mass FROM (
                SELECT DISTINCT ON (star_id) position, mass 
                FROM star_snapshots 
                WHERE dim_space = %s AND snapshot_id = %s 
                ORDER BY star_id, id DESC
            ) as sub;
        """, (dim, min_snap))
        initial_rows = cursor.fetchall()

        # Load final positions and masses
        cursor.execute("""
            SELECT position, mass FROM (
                SELECT DISTINCT ON (star_id) position, mass 
                FROM star_snapshots 
                WHERE dim_space = %s AND snapshot_id = %s 
                ORDER BY star_id, id DESC
            ) as sub;
        """, (dim, max_snap))
        final_rows = cursor.fetchall()

        if not initial_rows or not final_rows:
            print(f"  Skipping {dim}D: Insufficient data rows.")
            continue

        pos_init = np.array([r[0] for r in initial_rows])
        mass_init = np.array([r[1] for r in initial_rows])

        pos_final = np.array([r[0] for r in final_rows])
        mass_final = np.array([r[1] for r in final_rows])

        # Barycenter computation
        center_init = np.sum(pos_init * mass_init[:, np.newaxis], axis=0) / np.sum(mass_init)
        center_final = np.sum(pos_final * mass_final[:, np.newaxis], axis=0) / np.sum(mass_final)

        dist_init = np.sqrt(np.sum((pos_init - center_init)**2, axis=1))
        dist_final = np.sqrt(np.sum((pos_final - center_final)**2, axis=1))

        # Radial slope profiles
        r_i, alpha_i, err_i = compute_radial_slope_profile(dist_init, mass_init, args.num_radial_bins, args.m_min, args.m_max)
        r_f, alpha_f, err_f = compute_radial_slope_profile(dist_final, mass_final, args.num_radial_bins, args.m_min, args.m_max)

        dim_data[dim] = {
            "r_init": r_i, "alpha_init": alpha_i, "err_init": err_i,
            "r_final": r_f, "alpha_final": alpha_f, "err_final": err_f,
            "max_snap": max_snap
        }

    cursor.close()
    conn.close()

    if not dim_data:
        print("No valid dimension data processed.")
        return

    # 1. Multi-panel Subplot Grid
    n_dims = len(dim_data)
    cols = min(3, n_dims)
    rows = math.ceil(n_dims / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows), sharex=False, sharey=False)
    axes_flat = axes.flatten() if n_dims > 1 else [axes]

    for idx, (dim, data) in enumerate(dim_data.items()):
        ax = axes_flat[idx]

        ax.errorbar(data["r_init"], data["alpha_init"], yerr=data["err_init"],
                    fmt='bo:', label='Initial Slope (t=0)', linewidth=1.5, capsize=3, alpha=0.6)
        ax.errorbar(data["r_final"], data["alpha_final"], yerr=data["err_final"],
                    fmt='ro-', label=f'Final Slope (Snap {data["max_snap"]})', linewidth=2, capsize=3)

        ax.axhline(y=2.35, color='k', linestyle='--', alpha=0.7, label='Salpeter (α=2.35)')

        ax.set_xscale('log')
        ax.set_title(f'Dimension N = {dim}D', fontsize=12, fontweight='bold')
        ax.set_xlabel('Radial Distance R (pc)', fontsize=10)
        ax.set_ylabel('Mass Function Slope (α)', fontsize=10)
        ax.grid(True, which="both", ls="--", alpha=0.4)
        ax.legend(fontsize=8, loc='best')

    # Hide extra subplot panels
    for j in range(idx + 1, len(axes_flat)):
        fig.delaxes(axes_flat[j])

    plt.suptitle('Radial Mass Function Slope Profile α(R) Across Dimensions', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(args.output_grid, dpi=300, bbox_inches='tight')
    print(f"\nMulti-panel mass function slope profile saved to '{args.output_grid}'!")

    # 2. Overlay Comparison Plot for Final Slopes
    plt.figure(figsize=(10, 7))
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(dim_data)))

    plt.axhline(y=2.35, color='black', linestyle='--', linewidth=1.5, label='Salpeter Slope (α = 2.35)')

    for (dim, data), col in zip(dim_data.items(), colors):
        plt.errorbar(data["r_final"], data["alpha_final"], yerr=data["err_final"],
                     fmt='o-', color=col, label=f'N = {dim}D Final', linewidth=2, capsize=3, markersize=5)

    plt.xscale('log')
    plt.xlabel('Radial Distance R from Center (pc, log scale)', fontsize=12)
    plt.ylabel('Final Mass Function Slope α (dN/dM ∝ M^-α)', fontsize=12)
    plt.title('Comparison of Final Radial Mass Function Slope α(R) Across Dimensions', fontsize=13, fontweight='bold')
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.legend(fontsize=10, loc='best')

    plt.savefig(args.output_overlay, dpi=300, bbox_inches='tight')
    print(f"Overlay comparison plot saved to '{args.output_overlay}'!")

if __name__ == "__main__":
    main()
