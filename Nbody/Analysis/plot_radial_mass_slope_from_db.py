#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def compute_mass_slope(masses, m_min, m_max, num_bins=12):
    """
    Computes the Salpeter mass function slope alpha (dN/dM propto M^-alpha)
    for a given set of masses. Returns the fitted slope and its standard error.
    """
    if len(masses) < 20:
        return np.nan, np.nan
        
    bins = np.logspace(np.log10(m_min), np.log10(m_max), num_bins)
    counts, _ = np.histogram(masses, bins=bins)
    
    # dN/dM
    dn_dm = counts / (bins[1:] - bins[:-1])
    bin_centers = np.sqrt(bins[1:] * bins[:-1])
    
    valid = dn_dm > 0
    if np.sum(valid) < 3:
        return np.nan, np.nan
        
    try:
        # Fit with covariance matrix to get standard error
        slope, cov = np.polyfit(np.log10(bin_centers[valid]), np.log10(dn_dm[valid]), 1, cov=True)
        slope_err = np.sqrt(cov[0, 0])
        fit_slope = slope[0]
    except Exception:
        # Fallback if fit fails or cov cannot be computed
        try:
            fit_slope = np.polyfit(np.log10(bin_centers[valid]), np.log10(dn_dm[valid]), 1)[0]
            slope_err = np.nan
        except Exception:
            return np.nan, np.nan
            
    return -fit_slope, slope_err

def main():
    parser = argparse.ArgumentParser(description="Plot the initial and final mass function slope as a function of radial distance.")
    parser.add_argument("--dim", type=int, help="Spatial dimension of the simulation (e.g., 2, 3, 4, 5, 6, 100)")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--m-min", type=float, default=0.1, help="Minimum stellar mass in solar masses (default: 0.1)")
    parser.add_argument("--m-max", type=float, default=50.0, help="Maximum stellar mass in solar masses (default: 50.0)")
    parser.add_argument("--num-radial-bins", type=int, default=8, help="Number of radial bins for the slope profile (default: 8)")
    parser.add_argument("--num-mass-bins", type=int, default=12, help="Number of mass bins for fitting the slope (default: 12)")
    parser.add_argument("--linear", action="store_true", help="Plot in linear scale (default: log-X scale)")
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

    # Helper function to compute the radial profile of slopes using equal-number bins (quantiles)
    def compute_radial_slope_profile(dists, masses, num_radial_bins):
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
                # Use median distance as the radius of the bin
                bin_radii.append(np.median(bin_dists))
                alpha, alpha_err = compute_mass_slope(bin_masses, args.m_min, args.m_max, args.num_mass_bins)
                slopes.append(alpha)
                slope_errs.append(alpha_err)
                
        return np.array(bin_radii), np.array(slopes), np.array(slope_errs)

    print("Computing initial slope profile...")
    radii_init, slopes_init, slope_errs_init = compute_radial_slope_profile(dist_init, mass_init, args.num_radial_bins)
    
    print("Computing final slope profile...")
    radii_final, slopes_final, slope_errs_final = compute_radial_slope_profile(dist_final, mass_final, args.num_radial_bins)
 
    # Plot
    plt.figure(figsize=(10, 6))
    
    plt.errorbar(radii_init, slopes_init, yerr=slope_errs_init, fmt='bo-', label='Initial Slope Profile', linewidth=2, capsize=4)
    plt.errorbar(radii_final, slopes_final, yerr=slope_errs_final, fmt='ro-', label='Final Slope Profile', linewidth=2, capsize=4)
    
    # Horizontal line representing Salpeter slope (2.35)
    plt.axhline(y=2.35, color='k', linestyle='--', label='Initial Salpeter Slope (α=2.35)')
    
    if not args.linear:
        plt.xscale('log')
        plt.xlabel('Radial Distance from Center (pc, log scale)')
    else:
        plt.xlabel('Radial Distance from Center (pc)')
        
    plt.ylabel('Fitted Mass Function Slope (α)')
    scale_type = "Linear Scale" if args.linear else "Log-X Scale"
    plt.title(f'Radial Mass Function Slope Profile ({scale_type}, N={dim_space} dimensions)')
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    
    suffix = "_linear" if args.linear else ""
    plot_filename = f'radial_mass_slope_profile_{dim_space}D{suffix}_from_db.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Radial mass function slope profile ({scale_type}) successfully saved to '{plot_filename}'!")

if __name__ == "__main__":
    main()
