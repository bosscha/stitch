#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def compute_lagrangian_radii(dists, masses, percentiles):
    """
    Compute Lagrangian radii (mass percentiles) for a set of stellar distances and masses.
    """
    sort_idx = np.argsort(dists)
    sorted_dists = dists[sort_idx]
    sorted_masses = masses[sort_idx]
    cum_mass = np.cumsum(sorted_masses)
    cum_mass_frac = cum_mass / cum_mass[-1]
    
    # Interpolate to find the distance at each cumulative mass fraction
    radii = np.interp(np.array(percentiles) / 100.0, cum_mass_frac, sorted_dists)
    return radii

def main():
    parser = argparse.ArgumentParser(description="Plot the time series of radial distances of stars from the hypercluster database.")
    parser.add_argument("--dim", type=int, help="Spatial dimension of the simulation (e.g., 25, 100)")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--num-stars", type=int, default=50, help="Number of individual star trajectories to plot (default: 50)")
    parser.add_argument("--max-points", type=int, default=200, help="Maximum number of snapshots to analyze (default: 200)")
    parser.add_argument("--percentiles", type=str, default="10,30,50,75,90", help="Comma-separated list of Lagrangian radii percentiles to plot")
    parser.add_argument("--log-y", action="store_true", help="Plot radial distance on a logarithmic scale")
    args = parser.parse_args()

    # Prompt if dimension not provided
    dim_space = args.dim
    if dim_space is None:
        try:
            dim_space = int(input("Enter the spatial dimension (dim_space) to plot (e.g., 25, 100): "))
        except ValueError:
            print("Invalid input. Please enter an integer.")
            sys.exit(1)

    # Parse percentiles
    try:
        pct_list = [float(p.strip()) for p in args.percentiles.split(",")]
    except ValueError:
        print("Error: Invalid percentiles list. Please provide comma-separated numbers.")
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

    print(f"Connected successfully. Querying run history for dim_space = {dim_space}...")

    # 1. Find the maximum snapshot_id for this dimension (belongs to the latest run)
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

    # 2. Get the number of stars from the latest snapshot
    cursor.execute("""
        SELECT COUNT(*) 
        FROM star_snapshots 
        WHERE dim_space = %s AND snapshot_id = %s;
    """, (dim_space, max_snapshot_id))
    num_stars = cursor.fetchone()[0]

    # Find the minimum snapshot_id for this dimension
    cursor.execute("""
        SELECT MIN(snapshot_id)
        FROM star_snapshots
        WHERE dim_space = %s;
    """, (dim_space,))
    min_snapshot_id = cursor.fetchone()[0]

    # 3. Find the starting database ID for this run (minimum id of the latest min_snapshot_id)
    cursor.execute("""
        SELECT MIN(id) FROM (
            SELECT id FROM star_snapshots 
            WHERE dim_space = %s AND snapshot_id = %s 
            ORDER BY id DESC 
            LIMIT %s
        ) as sub;
    """, (dim_space, min_snapshot_id, num_stars))
    start_db_id = cursor.fetchone()[0]

    # 4. Find the step size dynamically
    cursor.execute("""
        SELECT MIN(snapshot_id) 
        FROM star_snapshots 
        WHERE dim_space = %s AND snapshot_id > %s AND id >= %s;
    """, (dim_space, min_snapshot_id, start_db_id))
    second_snap = cursor.fetchone()[0]
    step_size = (second_snap - min_snapshot_id) if second_snap is not None else 100

    # Generate full snapshot IDs list
    snap_ids = list(range(min_snapshot_id, max_snapshot_id + 1, step_size))
    total_snapshots = len(snap_ids)

    print(f"Detected snapshots up to ID {max_snapshot_id} (step size: {step_size}). Total potential snapshots: {total_snapshots}")

    # Downsample snapshots
    sample_interval = max(1, total_snapshots // args.max_points)
    sampled_snap_ids = snap_ids[::sample_interval]
    print(f"Downsampled to {len(sampled_snap_ids)} snapshots for analysis (interval: every {sample_interval} steps).")

    from collections import defaultdict

    times = []
    # To store radial distances for all stars at each snapshot
    # Shape: (num_sampled_snapshots, num_stars)
    rad_dists_all = []
    masses_all = []

    print("Loading all snapshot data in one query...")
    # Fetch all data in a single query
    cursor.execute("""
        SELECT snapshot_id, time_myr, star_id, position, mass
        FROM star_snapshots
        WHERE dim_space = %s AND id >= %s AND snapshot_id IN %s
        ORDER BY snapshot_id ASC, star_id ASC;
    """, (dim_space, start_db_id, tuple(sampled_snap_ids)))
    
    print("Fetching results...")
    all_rows = cursor.fetchall()
    print(f"Loaded {len(all_rows)} rows. Processing...")

    # Group by snapshot_id
    snap_data = defaultdict(list)
    for row in all_rows:
        snap_id = row[0]
        snap_data[snap_id].append(row)

    # Process each snapshot in chronological order
    for snap_id in sorted(snap_data.keys()):
        rows = snap_data[snap_id]
        
        # Sort by star_id to ensure consistent ordering of stars across snapshots
        rows.sort(key=lambda x: x[2])
        
        pos = np.array([row[3] for row in rows])
        mass = np.array([row[4] for row in rows])
        time_val = float(rows[0][1])
        
        # Compute center of mass (barycenter)
        cm = np.sum(pos * mass[:, np.newaxis], axis=0) / np.sum(mass)
        
        # Compute radial distance of each star from the center of mass
        dists = np.sqrt(np.sum((pos - cm)**2, axis=1))
        
        times.append(time_val)
        rad_dists_all.append(dists)
        masses_all.append(mass)

    print("Analysis complete. Processing data for plotting...")
    cursor.close()
    conn.close()

    rad_dists_all = np.array(rad_dists_all) # Shape: (num_snapshots, num_stars)
    masses_all = np.array(masses_all)       # Shape: (num_snapshots, num_stars)
    times = np.array(times)

    plt.figure(figsize=(12, 8))

    # 1. Plot individual star trajectories
    if args.num_stars > 0:
        # Choose a subset of star indices randomly (reproducible with seed)
        np.random.seed(42)
        sampled_star_indices = np.random.choice(num_stars, size=min(args.num_stars, num_stars), replace=False)
        
        print(f"Plotting trajectories of {len(sampled_star_indices)} sample stars...")
        for star_idx in sampled_star_indices:
            plt.plot(times, rad_dists_all[:, star_idx], color='lightgray', alpha=0.3, linewidth=0.6)

    # 2. Plot Lagrangian radii (percentiles)
    if pct_list:
        print(f"Plotting Lagrangian radii for percentiles: {pct_list}...")
        
        # Calculate Lagrangian radii at each snapshot
        lagrangian_radii = []
        for i in range(len(times)):
            snap_radii = compute_lagrangian_radii(rad_dists_all[i], masses_all[i], pct_list)
            lagrangian_radii.append(snap_radii)
        
        lagrangian_radii = np.array(lagrangian_radii) # Shape: (num_snapshots, len(pct_list))
        
        # Plot each Lagrangian radius line with a distinct color from a colormap
        cmap = plt.get_cmap('plasma')
        colors = [cmap(i) for i in np.linspace(0.1, 0.9, len(pct_list))]
        
        for p_idx, pct in enumerate(pct_list):
            plt.plot(times, lagrangian_radii[:, p_idx], color=colors[p_idx], 
                     label=f'{pct}% Lagrangian Radius', linewidth=2.5)

    if args.log_y:
        plt.yscale('log')
        plt.ylabel('Radial Distance R from Center of Mass (pc, log scale)')
    else:
        plt.ylabel('Radial Distance R from Center of Mass (pc)')

    plt.xlabel('Time (Myr)')
    plt.title(f'Time Series of Stellar Radial Distances ({dim_space}D)')
    plt.grid(True, which="both", linestyle='--', alpha=0.5)
    plt.legend(loc='upper right', framealpha=0.9)
    
    suffix = "_log" if args.log_y else ""
    plot_filename = f'radial_distance_time_series_{dim_space}D{suffix}.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Plot successfully saved to '{plot_filename}'!")

if __name__ == "__main__":
    main()
