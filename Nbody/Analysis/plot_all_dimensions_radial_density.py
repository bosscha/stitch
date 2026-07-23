#!/usr/bin/env python3
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def volume_nd(r, dim):
    """Calculate the volume of an n-dimensional sphere of radius r."""
    return (math.pi**(dim / 2.0) / math.gamma(dim / 2.0 + 1.0)) * (r**dim)

def main():
    parser = argparse.ArgumentParser(description="Plot initial vs final radial density (N-dimensional volume) for all dimensions from hypercluster database.")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--output", type=str, default="radial_density_all_dims.png", help="Output plot filename")
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
            print(f"  Skipping {dim}D: No snapshot data.")
            continue

        # Count distinct stars
        cursor.execute("SELECT COUNT(DISTINCT star_id) FROM star_snapshots WHERE dim_space = %s AND snapshot_id = %s;", (dim, max_snap))
        num_stars = cursor.fetchone()[0]

        if num_stars == 0:
            print(f"  Skipping {dim}D: 0 stars found.")
            continue

        # Load initial positions & masses (using DISTINCT ON star_id for single run)
        cursor.execute("""
            SELECT position, mass FROM (
                SELECT DISTINCT ON (star_id) position, mass 
                FROM star_snapshots 
                WHERE dim_space = %s AND snapshot_id = %s 
                ORDER BY star_id, id DESC
            ) as sub;
        """, (dim, min_snap))
        initial_rows = cursor.fetchall()

        # Load final positions & masses
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

        # Calculate N-dimensional barycenter (Center of Mass)
        center_init = np.sum(pos_init * mass_init[:, np.newaxis], axis=0) / np.sum(mass_init)
        center_final = np.sum(pos_final * mass_final[:, np.newaxis], axis=0) / np.sum(mass_final)

        # N-dimensional radial distances
        dist_init = np.sqrt(np.sum((pos_init - center_init)**2, axis=1))
        dist_final = np.sqrt(np.sum((pos_final - center_final)**2, axis=1))

        # Log-spaced radial bins
        max_r = max(np.max(dist_init), np.max(dist_final))
        min_r = max(min(np.min(dist_init[dist_init > 0]) if np.any(dist_init > 0) else 1e-2,
                        np.min(dist_final[dist_final > 0]) if np.any(dist_final > 0) else 1e-2), 1e-3)

        bins = np.logspace(np.log10(min_r), np.log10(max_r), 40)
        counts_init, edges = np.histogram(dist_init, bins=bins)
        counts_final, _ = np.histogram(dist_final, bins=edges)

        # Volume of N-dimensional spherical shell: V_N(R_outer) - V_N(R_inner)
        v_outer = volume_nd(edges[1:], dim)
        v_inner = volume_nd(edges[:-1], dim)
        volumes = np.maximum(v_outer - v_inner, 1e-20)

        rho_init = counts_init / volumes
        rho_final = counts_final / volumes
        bin_centers = (edges[:-1] + edges[1:]) / 2.0

        err_init = np.sqrt(counts_init) / volumes
        err_final = np.sqrt(counts_final) / volumes

        dim_data[dim] = {
            "r": bin_centers,
            "rho_init": rho_init,
            "rho_final": rho_final,
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

    # Multi-panel grid layout
    n_dims = len(dim_data)
    cols = min(3, n_dims)
    rows = math.ceil(n_dims / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows), sharex=False, sharey=False)
    axes_flat = axes.flatten() if n_dims > 1 else [axes]

    for idx, (dim, data) in enumerate(dim_data.items()):
        ax = axes_flat[idx]

        valid_i = data["rho_init"] > 0
        valid_f = data["rho_final"] > 0

        ax.errorbar(data["r"][valid_i], data["rho_init"][valid_i], yerr=data["err_init"][valid_i],
                    label='Initial (t=0)', color='#1f77b4', fmt='o-', markersize=4, capsize=2, alpha=0.85)
        ax.errorbar(data["r"][valid_f], data["rho_final"][valid_f], yerr=data["err_final"][valid_f],
                    label=f'Final (Snap {data["max_snap"]})', color='#d62728', fmt='s--', markersize=4, capsize=2, alpha=0.85)

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_title(f'Dimension N = {dim}D', fontsize=12, fontweight='bold')
        ax.set_xlabel('Radial Distance R (pc)', fontsize=10)
        ax.set_ylabel(f'3D/ND Density $\\rho$ (Stars / pc$^{{{dim}}}$)', fontsize=10)
        ax.grid(True, which="both", ls="--", alpha=0.4)
        ax.legend(fontsize=9)

    # Hide extra axes if any
    for j in range(idx + 1, len(axes_flat)):
        fig.delaxes(axes_flat[j])

    plt.suptitle('Initial vs Final N-Dimensional Radial Density Comparison Across Dimensions', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches='tight')
    print(f"\nMulti-panel radial density comparison saved to '{args.output}'!")

if __name__ == "__main__":
    main()
