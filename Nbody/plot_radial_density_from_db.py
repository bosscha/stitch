#!/usr/bin/env python3
import sys
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def volume_nd(r, dim):
    """Calculate the volume of an n-dimensional sphere of radius r."""
    return (math.pi**(dim/2.0) / math.gamma(dim/2.0 + 1.0)) * (r**dim)

def main():
    parser = argparse.ArgumentParser(description="Plot initial and final radial density from the hypercluster database.")
    parser.add_argument("--dim", type=int, help="Spatial dimension of the simulation (e.g., 2, 3, 4, 5, 6, 100)")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    args = parser.parse_args()

    # Prompt if dimension not provided
    dim_space = args.dim
    if dim_space is None:
        try:
            dim_space = int(input("Enter the spatial dimension (dim_space) to plot (e.g., 2, 3, 4, 5, 6, 100): "))
        except ValueError:
            print("Invalid input. Please enter an integer.")
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
    
    # Compute center of positions (mean of positions)
    center_init = np.mean(pos_init, axis=0)
    center_final = np.mean(pos_final, axis=0)
    
    # Calculate radial distances from center of mass
    dist_init = np.sqrt(np.sum((pos_init - center_init)**2, axis=1))
    dist_final = np.sqrt(np.sum((pos_final - center_final)**2, axis=1))

    print(f"Initial positions center: {center_init}")
    print(f"Final positions center: {center_final}")
    print(f"Min/Max initial distance: {np.min(dist_init):.4f} / {np.max(dist_init):.4f}")
    print(f"Min/Max final distance: {np.min(dist_final):.4f} / {np.max(dist_final):.4f}")

    # Plot Setup
    plt.figure(figsize=(10, 6))
    
    # Use log-spaced bins for a better representation in log-log scale
    max_dist = max(np.max(dist_init), np.max(dist_final))
    min_dist_val = max(np.min(dist_init[dist_init > 0]) if np.any(dist_init > 0) else 1e-2, 1e-3)
    bins_log = np.logspace(np.log10(min_dist_val), np.log10(max_dist), 50)
    
    counts_init, bins_edges = np.histogram(dist_init, bins=bins_log)
    counts_final, _ = np.histogram(dist_final, bins=bins_edges)
    
    # Calculate volumes of n-dimensional shells
    volumes = volume_nd(bins_edges[1:], dim_space) - volume_nd(bins_edges[:-1], dim_space)
    volumes = np.maximum(volumes, 1e-15)
    
    density_init = counts_init / volumes
    density_final = counts_final / volumes
    bin_centers = (bins_edges[:-1] + bins_edges[1:]) / 2.0
    
    # Only plot where density > 0 to avoid log(0) issues in log-log plot
    valid_init = density_init > 0
    valid_final = density_final > 0
    
    plt.plot(bin_centers[valid_init], density_init[valid_init], label='Initial Density (DB)', color='b', marker='.', linestyle='-')
    plt.plot(bin_centers[valid_final], density_final[valid_final], label='Final Density (DB)', color='r', marker='.', linestyle='-')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Radial Distance from Center')
    plt.ylabel('Density (Stars / Volume)')
    plt.title(f'Initial vs Final Radial Density (N={dim_space} dimensions) - from DB')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    
    plot_filename = f'radial_density_{dim_space}D_from_db.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Radial density plot successfully saved to '{plot_filename}'!")

if __name__ == "__main__":
    main()
