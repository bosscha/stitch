#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def main():
    parser = argparse.ArgumentParser(description="Plot initial and final radial velocity dispersion from the hypercluster database.")
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

    # Load initial positions and velocities (snapshot_id = 0) from the latest run
    print("Loading initial snapshot positions and velocities...")
    cursor.execute("""
        SELECT position, velocity FROM (
            SELECT id, position, velocity 
            FROM star_snapshots 
            WHERE dim_space = %s AND snapshot_id = 0 
            ORDER BY id DESC 
            LIMIT %s
        ) as sub 
        ORDER BY id ASC;
    """, (dim_space, num_stars))
    initial_rows = cursor.fetchall()
    
    # Load final positions and velocities (snapshot_id = max_snapshot_id)
    print("Loading final snapshot positions and velocities...")
    cursor.execute("""
        SELECT position, velocity FROM (
            SELECT id, position, velocity 
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

    # Convert to numpy arrays
    pos_init = np.array([row[0] for row in initial_rows])
    vel_init = np.array([row[1] for row in initial_rows])
    
    pos_final = np.array([row[0] for row in final_rows])
    vel_final = np.array([row[1] for row in final_rows])
    
    def calculate_radial_velocities(pos, vel):
        # 1. Center of position (mean of positions)
        center_pos = np.mean(pos, axis=0)
        # 2. Bulk velocity of the cluster (mean of velocities)
        bulk_vel = np.mean(vel, axis=0)
        
        # 3. Relative positions and velocities
        rel_pos = pos - center_pos
        rel_vel = vel - bulk_vel
        
        # 4. Radial distances
        dists = np.sqrt(np.sum(rel_pos**2, axis=1))
        
        # Avoid division by zero for any star at the exact center
        eps = 1e-10
        dists_safe = np.where(dists == 0, eps, dists)
        
        # 5. Unit radial vector
        unit_r = rel_pos / dists_safe[:, np.newaxis]
        
        # 6. Radial velocity: dot product of relative velocity and unit radial vector
        v_radial = np.sum(rel_vel * unit_r, axis=1)
        
        return dists, v_radial

    print("Calculating initial radial velocities...")
    dist_init, v_rad_init = calculate_radial_velocities(pos_init, vel_init)
    
    print("Calculating final radial velocities...")
    dist_final, v_rad_final = calculate_radial_velocities(pos_final, vel_final)

    # Calculate velocity dispersion in radial bins
    max_dist = max(np.max(dist_init), np.max(dist_final))
    min_dist_val = max(np.min(dist_init[dist_init > 0]) if np.any(dist_init > 0) else 1e-2, 1e-3)
    bins_log = np.logspace(np.log10(min_dist_val), np.log10(max_dist), 40)
    
    def compute_dispersion(dists, v_rad, bins):
        bin_centers = []
        dispersion = []
        for i in range(len(bins)-1):
            lower, upper = bins[i], bins[i+1]
            mask = (dists >= lower) & (dists < upper)
            subset = v_rad[mask]
            if len(subset) >= 5: # Require at least 5 stars in a bin
                bin_centers.append((lower + upper) / 2.0)
                dispersion.append(np.std(subset))
        return np.array(bin_centers), np.array(dispersion)

    print("Computing dispersion profiles...")
    bin_centers_init, sigma_r_init = compute_dispersion(dist_init, v_rad_init, bins_log)
    bin_centers_final, sigma_r_final = compute_dispersion(dist_final, v_rad_final, bins_log)

    # Plot Setup
    plt.figure(figsize=(10, 6))
    
    plt.plot(bin_centers_init, sigma_r_init, label='Initial Radial σ_r (DB)', color='b', marker='.', linestyle='-')
    plt.plot(bin_centers_final, sigma_r_final, label='Final Radial σ_r (DB)', color='r', marker='.', linestyle='-')
    
    plt.xscale('log')
    plt.xlabel('Radial Distance from Center')
    plt.ylabel('Radial Velocity Dispersion σ_r')
    plt.title(f'Initial vs Final Radial Velocity Dispersion (N={dim_space} dimensions) - from DB')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    
    plot_filename = f'radial_velocity_dispersion_{dim_space}D_from_db.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Radial velocity dispersion plot successfully saved to '{plot_filename}'!")

if __name__ == "__main__":
    main()
