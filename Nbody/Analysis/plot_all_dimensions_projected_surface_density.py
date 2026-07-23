#!/usr/bin/env python3
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def main():
    parser = argparse.ArgumentParser(description="Plot initial and final projected 2D surface density comparison for all dimensions (N >= 2).")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--proj-x", type=int, default=0, help="First dimension index for 2D projection (default: 0)")
    parser.add_argument("--proj-y", type=int, default=1, help="Second dimension index for 2D projection (default: 1)")
    parser.add_argument("--min-dim", type=int, default=2, help="Minimum spatial dimension to include (default: 2)")
    parser.add_argument("--output", type=str, default="projected_surface_density_all_dims.png", help="Output plot filename")
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

    # Get available spatial dimensions >= min_dim
    cursor.execute("""
        SELECT DISTINCT dim_space 
        FROM star_snapshots 
        WHERE dim_space >= %s
        ORDER BY dim_space ASC;
    """, (args.min_dim,))
    dims = [row[0] for row in cursor.fetchall()]

    if not dims:
        print(f"No spatial dimensions >= {args.min_dim} found in table 'star_snapshots'.")
        cursor.close()
        conn.close()
        return

    print(f"Found {len(dims)} spatial dimensions to process: {dims}")

    dim_data = {}

    for dim in dims:
        print(f"\nProcessing dim_space = {dim}D...")
        # Get min and max snapshot_id
        cursor.execute("SELECT MIN(snapshot_id), MAX(snapshot_id) FROM star_snapshots WHERE dim_space = %s;", (dim,))
        min_snap, max_snap = cursor.fetchone()
        
        if min_snap is None or max_snap is None:
            print(f"  Skipping {dim}D: No snapshot range found.")
            continue

        # Number of stars per single run snapshot
        cursor.execute("SELECT COUNT(DISTINCT star_id) FROM star_snapshots WHERE dim_space = %s AND snapshot_id = %s;", (dim, max_snap))
        num_stars = cursor.fetchone()[0]

        if num_stars == 0:
            print(f"  Skipping {dim}D: 0 stars in final snapshot.")
            continue

        # Load initial positions & masses for single run
        cursor.execute("""
            SELECT position, mass FROM (
                SELECT DISTINCT ON (star_id) position, mass 
                FROM star_snapshots 
                WHERE dim_space = %s AND snapshot_id = %s 
                ORDER BY star_id, id DESC
            ) as sub;
        """, (dim, min_snap))
        initial_rows = cursor.fetchall()

        # Load final positions & masses for single run
        cursor.execute("""
            SELECT position, mass FROM (
                SELECT DISTINCT ON (star_id) position, mass 
                FROM star_snapshots 
                WHERE dim_space = %s AND snapshot_id = %s 
                ORDER BY star_id, id DESC
            ) as sub;
        """, (dim, max_snap))
        final_rows = cursor.fetchall()

        pos_init = np.array([r[0] for r in initial_rows])
        mass_init = np.array([r[1] for r in initial_rows])

        pos_final = np.array([r[0] for r in final_rows])
        mass_final = np.array([r[1] for r in final_rows])

        # 2D Projections
        proj_x = min(args.proj_x, dim - 1)
        proj_y = min(args.proj_y, dim - 1)
        if proj_x == proj_y:
            proj_y = (proj_x + 1) % dim

        pos_p_init = pos_init[:, [proj_x, proj_y]]
        pos_p_final = pos_final[:, [proj_x, proj_y]]

        center_init = np.sum(pos_p_init * mass_init[:, np.newaxis], axis=0) / np.sum(mass_init)
        center_final = np.sum(pos_p_final * mass_final[:, np.newaxis], axis=0) / np.sum(mass_final)

        dist_init = np.sqrt(np.sum((pos_p_init - center_init)**2, axis=1))
        dist_final = np.sqrt(np.sum((pos_p_final - center_final)**2, axis=1))

        # Compute log-spaced surface density profile
        max_r = max(np.max(dist_init), np.max(dist_final))
        min_r = max(min(np.min(dist_init[dist_init > 0]) if np.any(dist_init > 0) else 1e-2,
                        np.min(dist_final[dist_final > 0]) if np.any(dist_final > 0) else 1e-2), 1e-3)

        bins = np.logspace(np.log10(min_r), np.log10(max_r), 40)
        counts_init, edges = np.histogram(dist_init, bins=bins)
        counts_final, _ = np.histogram(dist_final, bins=edges)

        areas = math.pi * (edges[1:]**2 - edges[:-1]**2)
        areas = np.maximum(areas, 1e-15)

        sd_init = counts_init / areas
        sd_final = counts_final / areas
        bin_centers = (edges[:-1] + edges[1:]) / 2.0

        err_init = np.sqrt(counts_init) / areas
        err_final = np.sqrt(counts_final) / areas

        dim_data[dim] = {
            "r": bin_centers,
            "sd_init": sd_init,
            "sd_final": sd_final,
            "err_init": err_init,
            "err_final": err_final,
            "min_snap": min_snap,
            "max_snap": max_snap,
            "num_stars": num_stars
        }

    cursor.close()
    conn.close()

    if not dim_data:
        print("No valid dimension data extracted.")
        return

    # Plot creation: Multi-panel grid + Summary comparison overlay
    n_dims = len(dim_data)
    cols = min(3, n_dims)
    rows = math.ceil(n_dims / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), sharex=False, sharey=False)
    axes_flat = axes.flatten() if n_dims > 1 else [axes]

    for idx, (dim, data) in enumerate(dim_data.items()):
        ax = axes_flat[idx]
        
        valid_i = data["sd_init"] > 0
        valid_f = data["sd_final"] > 0

        ax.errorbar(data["r"][valid_i], data["sd_init"][valid_i], yerr=data["err_init"][valid_i],
                    label='Initial (t=0)', color='#1f77b4', fmt='o-', markersize=4, capsize=2, alpha=0.85)
        ax.errorbar(data["r"][valid_f], data["sd_final"][valid_f], yerr=data["err_final"][valid_f],
                    label=f'Final (Snap {data["max_snap"]})', color='#d62728', fmt='s--', markersize=4, capsize=2, alpha=0.85)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(f'Dimension N = {dim}D', fontsize=12, fontweight='bold')
        ax.set_xlabel('Projected Radius R (pc)', fontsize=10)
        ax.set_ylabel('Surface Density Σ (Stars / pc²)', fontsize=10)
        ax.grid(True, which="both", ls="--", alpha=0.4)
        ax.legend(fontsize=9)

    # Hide unused subplots if any
    for j in range(idx + 1, len(axes_flat)):
        fig.delaxes(axes_flat[j])

    plt.suptitle('Comparison of Initial vs Final Projected 2D Surface Density Across Dimensions (N ≥ 2)', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    grid_filename = args.output
    plt.savefig(grid_filename, dpi=300, bbox_inches='tight')
    print(f"\nMulti-panel comparison plot saved to '{grid_filename}'!")

    # Summary Overlay Plot: All dimensions final vs initial
    plt.figure(figsize=(10, 7))
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(dim_data)))

    for (dim, data), col in zip(dim_data.items(), colors):
        valid_f = data["sd_final"] > 0
        valid_i = data["sd_init"] > 0
        plt.plot(data["r"][valid_i], data["sd_init"][valid_i], linestyle=':', color=col, alpha=0.4)
        plt.plot(data["r"][valid_f], data["sd_final"][valid_f], label=f'{dim}D Final', linestyle='-', color=col, linewidth=2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Projected Radius R from Center (pc)', fontsize=12)
    plt.ylabel('Surface Density Σ (Stars / pc²)', fontsize=12)
    plt.title('Overlay Comparison: Final Radial Projected Surface Density by Dimension (N ≥ 2)', fontsize=13, fontweight='bold')
    plt.grid(True, which="both", ls="--", alpha=0.4)
    plt.legend(fontsize=10, loc='best')
    overlay_filename = "projected_surface_density_overlay_all_dims.png"
    plt.savefig(overlay_filename, dpi=300, bbox_inches='tight')
    print(f"Overlay comparison plot saved to '{overlay_filename}'!")

if __name__ == "__main__":
    main()
