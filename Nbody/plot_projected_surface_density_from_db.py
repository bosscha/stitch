#!/usr/bin/env python3
import sys
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def main():
    parser = argparse.ArgumentParser(description="Plot initial and final projected 2D surface density from the hypercluster database.")
    parser.add_argument("--dim", type=int, help="Spatial dimension of the simulation (e.g., 2, 3, 4, 5, 6, 100)")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--proj-x", type=int, default=0, help="First dimension index for 2D projection (default: 0)")
    parser.add_argument("--proj-y", type=int, default=1, help="Second dimension index for 2D projection (default: 1)")
    parser.add_argument("--linear", action="store_true", help="Plot in linear scale with linearly-spaced bins (default: log-log scale)")
    args = parser.parse_args()

    # Prompt if dimension not provided
    dim_space = args.dim
    if dim_space is None:
        try:
            dim_space = int(input("Enter the spatial dimension (dim_space) to plot (e.g., 2, 3, 4, 5, 6, 100): "))
        except ValueError:
            print("Invalid input. Please enter an integer.")
            sys.exit(1)

    if dim_space < 2:
        print("Error: Spatial dimension must be at least 2 to compute a 2D projected surface density.")
        sys.exit(1)

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
        sys.exit(1)

    print(f"Connected successfully. Querying latest run for dim_space = {dim_space}...")

    # Find the maximum snapshot_id for this dimension (belongs to the latest run)
    cursor.execute("""
        SELECT MAX(snapshot_id) 
        FROM star_snapshots 
        WHERE dim_space = %s;
    """, (dim_space,))
    
    max_snapshot_id = cursor.fetchone()[0]
    if max_snapshot_id is None:
        print(f"No simulation data found in the database for dim_space = {dim_space}.")
        cursor.close()
        conn.close()
        sys.exit(1)

    # Determine the number of stars in this run by counting stars in the final snapshot
    cursor.execute("""
        SELECT COUNT(*) 
        FROM star_snapshots 
        WHERE dim_space = %s AND snapshot_id = %s;
    """, (dim_space, max_snapshot_id))
    num_stars = cursor.fetchone()[0]

    if num_stars == 0:
        print(f"Error: Found 0 stars for the final snapshot {max_snapshot_id}.")
        cursor.close()
        conn.close()
        sys.exit(1)

    print(f"Found simulation snapshots up to {max_snapshot_id} with {num_stars} stars per snapshot.")

    # Load initial positions (snapshot_id = 0) from the latest run
    print("Loading initial snapshot positions...")
    cursor.execute("""
        SELECT position FROM (
            SELECT id, position 
            FROM star_snapshots 
            WHERE dim_space = %s AND snapshot_id = 0 
            ORDER BY id DESC 
            LIMIT %s
        ) as sub 
        ORDER BY id ASC;
    """, (dim_space, num_stars))
    initial_rows = cursor.fetchall()
    
    # Load final positions (snapshot_id = max_snapshot_id)
    print("Loading final snapshot positions...")
    cursor.execute("""
        SELECT position FROM (
            SELECT id, position 
            FROM star_snapshots 
            WHERE dim_space = %s AND snapshot_id = %s 
            ORDER BY id DESC 
            LIMIT %s
        ) as sub 
        ORDER BY id ASC;
    """, (dim_space, max_snapshot_id, num_stars))
    final_rows = cursor.fetchall()

    cursor.close()
    conn.close()

    if len(initial_rows) < num_stars or len(final_rows) < num_stars:
        print(f"Warning: Expected {num_stars} stars, but loaded {len(initial_rows)} initial and {len(final_rows)} final stars.")
        if not initial_rows or not final_rows:
            print("Error: Could not retrieve snapshot data.")
            sys.exit(1)

    # Convert positions to numpy arrays
    pos_init = np.array([row[0] for row in initial_rows])
    pos_final = np.array([row[0] for row in final_rows])
    
    # Project positions onto the selected 2D plane (default: dimensions 0 and 1)
    pos_proj_init = pos_init[:, [args.proj_x, args.proj_y]]
    pos_proj_final = pos_final[:, [args.proj_x, args.proj_y]]
    
    # Compute center of projected positions (mean)
    center_proj_init = np.mean(pos_proj_init, axis=0)
    center_proj_final = np.mean(pos_proj_final, axis=0)
    
    # Calculate projected radial distances from the projected center
    dist_proj_init = np.sqrt(np.sum((pos_proj_init - center_proj_init)**2, axis=1))
    dist_proj_final = np.sqrt(np.sum((pos_proj_final - center_proj_final)**2, axis=1))

    print(f"Projected center (dimensions {args.proj_x} & {args.proj_y}) - Initial: {center_proj_init}, Final: {center_proj_final}")
    print(f"Projected distance range - Initial: {np.min(dist_proj_init):.4f} to {np.max(dist_proj_init):.4f}")
    print(f"Projected distance range - Final: {np.min(dist_proj_final):.4f} to {np.max(dist_proj_final):.4f}")

    # Plot Setup
    plt.figure(figsize=(10, 6))
    
    max_dist = max(np.max(dist_proj_init), np.max(dist_proj_final))
    
    if args.linear:
        # Use linearly-spaced bins for linear scale
        bins = np.linspace(0.0, max_dist, 50)
    else:
        # Use log-spaced bins for log-log scale
        min_dist_val = max(np.min(dist_proj_init[dist_proj_init > 0]) if np.any(dist_proj_init > 0) else 1e-2, 1e-3)
        bins = np.logspace(np.log10(min_dist_val), np.log10(max_dist), 50)
        
    counts_init, bins_edges = np.histogram(dist_proj_init, bins=bins)
    counts_final, _ = np.histogram(dist_proj_final, bins=bins_edges)
    
    # Calculate areas of 2D concentric circular rings (Area = pi * (R_outer^2 - R_inner^2))
    areas = math.pi * (bins_edges[1:]**2 - bins_edges[:-1]**2)
    areas = np.maximum(areas, 1e-15)
    
    surf_density_init = counts_init / areas
    surf_density_final = counts_final / areas
    bin_centers = (bins_edges[:-1] + bins_edges[1:]) / 2.0
    
    # Plot density
    if args.linear:
        plt.plot(bin_centers, surf_density_init, label='Initial Surface Density (DB)', color='b', marker='.', linestyle='-')
        plt.plot(bin_centers, surf_density_final, label='Final Surface Density (DB)', color='r', marker='.', linestyle='-')
    else:
        # Only plot where density > 0 to avoid log(0) issues in log-log plot
        valid_init = surf_density_init > 0
        valid_final = surf_density_final > 0
        plt.plot(bin_centers[valid_init], surf_density_init[valid_init], label='Initial Surface Density (DB)', color='b', marker='.', linestyle='-')
        plt.plot(bin_centers[valid_final], surf_density_final[valid_final], label='Final Surface Density (DB)', color='r', marker='.', linestyle='-')
        plt.xscale('log')
        plt.yscale('log')
        
    plt.xlabel('Projected Radial Distance R from Center')
    plt.ylabel('Surface Density Σ (Stars / Area)')
    
    scale_type = "Linear Scale" if args.linear else "Log-Log Scale"
    plt.title(f'Initial vs Final Projected Surface Density ({scale_type}, Projected from {dim_space}D to 2D)')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    
    suffix = "_linear" if args.linear else ""
    plot_filename = f'projected_surface_density_{dim_space}D{suffix}_from_db.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Projected surface density plot ({scale_type}) successfully saved to '{plot_filename}'!")

if __name__ == "__main__":
    main()
