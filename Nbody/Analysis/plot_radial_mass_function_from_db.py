#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def compute_mass_slope(masses, m_min, m_max, num_bins=15):
    """
    Computes the Salpeter mass function slope alpha (dN/dM propto M^-alpha)
    for a given set of masses. Returns the slope and the bin data for plotting.
    """
    if len(masses) < 15:
        return np.nan, None, None, None
        
    bins = np.logspace(np.log10(m_min), np.log10(m_max), num_bins)
    counts, _ = np.histogram(masses, bins=bins)
    
    # dN/dM
    dn_dm = counts / (bins[1:] - bins[:-1])
    bin_centers = np.sqrt(bins[1:] * bins[:-1])
    
    valid = dn_dm > 0
    if np.sum(valid) < 3:
        return np.nan, None, None, None
        
    slope, _ = np.polyfit(np.log10(bin_centers[valid]), np.log10(dn_dm[valid]), 1)
    return -slope, bin_centers, dn_dm, valid

def main():
    parser = argparse.ArgumentParser(description="Plot the initial and final radial mass function from the hypercluster database.")
    parser.add_argument("--dim", type=int, help="Spatial dimension of the simulation (e.g., 2, 3, 4, 5, 6, 100)")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--m-min", type=float, default=0.1, help="Minimum stellar mass in solar masses (default: 0.1)")
    parser.add_argument("--m-max", type=float, default=50.0, help="Maximum stellar mass in solar masses (default: 50.0)")
    parser.add_argument("--num-bins", type=int, default=12, help="Number of mass bins for the histogram (default: 12)")
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

    # Find the minimum snapshot_id for this dimension
    cursor.execute("""
        SELECT MIN(snapshot_id)
        FROM star_snapshots
        WHERE dim_space = %s;
    """, (dim_space,))
    min_snapshot_id = cursor.fetchone()[0]

    # Find the starting database ID for this run (minimum id of the latest min_snapshot_id)
    cursor.execute("""
        SELECT MIN(id) FROM (
            SELECT id FROM star_snapshots 
            WHERE dim_space = %s AND snapshot_id = %s 
            ORDER BY id DESC 
            LIMIT %s
        ) as sub;
    """, (dim_space, min_snapshot_id, num_stars))
    start_db_id = cursor.fetchone()[0]

    print(f"Found simulation snapshots up to {max_snapshot_id} with {num_stars} stars per snapshot.")

    # Load initial positions and masses (snapshot_id = min_snapshot_id) from the latest run
    print("Loading initial snapshot positions and masses...")
    cursor.execute("""
        SELECT position, mass FROM (
            SELECT id, position, mass 
            FROM star_snapshots 
            WHERE dim_space = %s AND snapshot_id = %s AND id >= %s
            ORDER BY id DESC 
            LIMIT %s
        ) as sub 
        ORDER BY id ASC;
    """, (dim_space, min_snapshot_id, start_db_id, num_stars))
    initial_rows = cursor.fetchall()
    
    # Load final positions and masses (snapshot_id = max_snapshot_id)
    print("Loading final snapshot positions and masses...")
    cursor.execute("""
        SELECT position, mass FROM (
            SELECT id, position, mass 
            FROM star_snapshots 
            WHERE dim_space = %s AND snapshot_id = %s AND id >= %s
            ORDER BY id DESC 
            LIMIT %s
        ) as sub 
        ORDER BY id ASC;
    """, (dim_space, max_snapshot_id, start_db_id, num_stars))
    final_rows = cursor.fetchall()

    cursor.close()
    conn.close()

    # Convert to numpy arrays
    pos_init = np.array([row[0] for row in initial_rows])
    mass_init = np.array([row[1] for row in initial_rows])
    
    pos_final = np.array([row[0] for row in final_rows])
    mass_final = np.array([row[1] for row in final_rows])
    
    # Compute center of mass (barycenter)
    center_init = np.sum(pos_init * mass_init[:, np.newaxis], axis=0) / np.sum(mass_init)
    center_final = np.sum(pos_final * mass_final[:, np.newaxis], axis=0) / np.sum(mass_final)
    
    # Calculate radial distances
    dist_init = np.sqrt(np.sum((pos_init - center_init)**2, axis=1))
    dist_final = np.sqrt(np.sum((pos_final - center_final)**2, axis=1))

    # Define radial zones based on quantiles or physical boundaries (e.g. 50% radius)
    # R50 divides the cluster into equal numbers of stars (inner vs outer half)
    r50_init = np.quantile(dist_init, 0.5)
    r50_final = np.quantile(dist_final, 0.5)
    
    print(f"Initial R50 (50% mass radius): {r50_init:.4f} pc")
    print(f"Final R50 (50% mass radius): {r50_final:.4f} pc")

    # Masks for inner vs outer stars
    mask_inner_init = dist_init <= r50_init
    mask_outer_init = dist_init > r50_init
    
    mask_inner_final = dist_final <= r50_final
    mask_outer_final = dist_final > r50_final

    # Calculate mass functions and slopes
    slope_in_init, bins_in_init, dn_in_init, val_in_init = compute_mass_slope(mass_init[mask_inner_init], args.m_min, args.m_max, args.num_bins)
    slope_out_init, bins_out_init, dn_out_init, val_out_init = compute_mass_slope(mass_init[mask_outer_init], args.m_min, args.m_max, args.num_bins)
    
    slope_in_final, bins_in_final, dn_in_final, val_in_final = compute_mass_slope(mass_final[mask_inner_final], args.m_min, args.m_max, args.num_bins)
    slope_out_final, bins_out_final, dn_out_final, val_out_final = compute_mass_slope(mass_final[mask_outer_final], args.m_min, args.m_max, args.num_bins)

    print(f"Initial slopes - Inner (R <= R50): {slope_in_init:.3f}, Outer (R > R50): {slope_out_init:.3f}")
    print(f"Final slopes   - Inner (R <= R50): {slope_in_final:.3f}, Outer (R > R50): {slope_out_final:.3f}")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

    # Initial snapshot plotting
    ax1.set_title("Initial Mass Function (by Radial Zone)")
    if bins_in_init is not None:
        ax1.loglog(bins_in_init[val_in_init], dn_in_init[val_in_init], 'bo-', label=f'Inner Core (α={slope_in_init:.2f})')
    if bins_out_init is not None:
        ax1.loglog(bins_out_init[val_out_init], dn_out_init[val_out_init], 'co--', label=f'Outer Envelope (α={slope_out_init:.2f})')
    ax1.set_xlabel('Stellar Mass (M_sun)')
    ax1.set_ylabel('dN/dM')
    ax1.legend()
    ax1.grid(True, which="both", ls="--", alpha=0.5)

    # Final snapshot plotting
    ax2.set_title("Final Mass Function (by Radial Zone)")
    if bins_in_final is not None:
        ax2.loglog(bins_in_final[val_in_final], dn_in_final[val_in_final], 'ro-', label=f'Inner Core (α={slope_in_final:.2f})')
    if bins_out_final is not None:
        ax2.loglog(bins_out_final[val_out_final], dn_out_final[val_out_final], 'yo--', label=f'Outer Envelope (α={slope_out_final:.2f})')
    ax2.set_xlabel('Stellar Mass (M_sun)')
    ax2.legend()
    ax2.grid(True, which="both", ls="--", alpha=0.5)

    plt.tight_layout()
    plot_filename = f'radial_mass_function_{dim_space}D_from_db.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Radial mass function plot successfully saved to '{plot_filename}'!")

if __name__ == "__main__":
    main()
