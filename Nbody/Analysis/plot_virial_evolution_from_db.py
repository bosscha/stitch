#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def get_bound_subsystem(pos, vel, mass, G, softening, max_iter=10):
    """
    Iteratively removes unbound stars (E_i = T_i + V_i >= 0) to find the bound subsystem,
    and returns its indices, total kinetic energy, total potential energy, and bound count.
    """
    idx_bound = np.arange(len(mass))
    
    T_tot = 0.0
    V_tot = 0.0
    
    for iteration in range(max_iter):
        N_curr = len(idx_bound)
        if N_curr < 3:
            return idx_bound, 0.0, 0.0
            
        pos_curr = pos[idx_bound]
        vel_curr = vel[idx_bound]
        mass_curr = mass[idx_bound]
        
        # Center of mass position & velocity of current bound subset
        total_m = np.sum(mass_curr)
        cm_pos = np.sum(pos_curr * mass_curr[:, np.newaxis], axis=0) / total_m
        cm_vel = np.sum(vel_curr * mass_curr[:, np.newaxis], axis=0) / total_m
        
        # Relative coordinates
        rel_pos = pos_curr - cm_pos
        rel_vel = vel_curr - cm_vel
        
        # Kinetic energy of each star relative to bound CM
        T_i = 0.5 * mass_curr * np.sum(rel_vel**2, axis=1)
        
        # Potential of each star due to other bound stars
        pos_sq = np.sum(rel_pos**2, axis=1)
        dist_sq = pos_sq[:, np.newaxis] + pos_sq[np.newaxis, :] - 2 * np.dot(rel_pos, rel_pos.T)
        dist_sq = np.maximum(dist_sq, 0.0)  # Avoid negative values due to floating-point precision
        inv_dist = 1.0 / np.sqrt(dist_sq + softening**2)
        np.fill_diagonal(inv_dist, 0)
        
        Phi_i = -G * np.sum(mass_curr * inv_dist, axis=1)
        V_i = mass_curr * Phi_i
        
        # Total energy of each star
        E_i = T_i + V_i
        
        # Identify bound stars
        bound_mask = E_i < 0
        N_next = np.sum(bound_mask)
        
        T_tot = np.sum(T_i)
        V_tot = 0.5 * np.sum(V_i) # 0.5 to avoid double counting pairs
        
        if N_next == N_curr:
            break
            
        idx_bound = idx_bound[bound_mask]
        
    return idx_bound, T_tot, V_tot

def main():
    parser = argparse.ArgumentParser(description="Plot the evolution of the virial ratio of bound stars from the hypercluster database.")
    parser.add_argument("--dim", type=int, help="Spatial dimension of the simulation (e.g., 2, 3, 4, 5, 6, 100)")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--softening", type=float, default=0.01, help="Softening parameter used in simulation (default: 0.01)")
    parser.add_argument("--g-constant", type=float, default=1.0, help="Gravitational constant G (default: 1.0)")
    parser.add_argument("--max-points", type=int, default=200, help="Maximum number of snapshots to analyze to keep runtimes fast (default: 200)")
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

    times = []
    virial_ratios = []
    virial_ratios_all = []
    bound_fractions = []

    # Get N_STARS
    cursor.execute("""
        SELECT COUNT(*) 
        FROM star_snapshots 
        WHERE dim_space = %s AND snapshot_id = %s AND id >= %s;
    """, (dim_space, sampled_snap_ids[0], start_db_id))
    n_stars = cursor.fetchone()[0]
    print(f"Number of stars per snapshot: {n_stars}")

    # Process each snapshot
    for idx, snap_id in enumerate(sampled_snap_ids):
        sys.stdout.write(f"\rAnalyzing snapshot {idx+1}/{len(sampled_snap_ids)} (ID: {snap_id})...")
        sys.stdout.flush()
        
        # Load positions, velocities, masses, and time for this snapshot
        cursor.execute("""
            SELECT position, velocity, mass, time_myr FROM (
                SELECT id, position, velocity, mass, time_myr 
                FROM star_snapshots 
                WHERE dim_space = %s AND snapshot_id = %s AND id >= %s
                ORDER BY id DESC 
                LIMIT %s
            ) as sub 
            ORDER BY id ASC;
        """, (dim_space, snap_id, start_db_id, n_stars))
        
        rows = cursor.fetchall()
        if len(rows) == 0:
            continue
            
        pos = np.array([row[0] for row in rows])
        vel = np.array([row[1] for row in rows])
        mass = np.array([row[2] for row in rows])
        # Normalize masses to N-body units (sum to 1.0) to match N-body positions and velocities
        mass = mass / np.sum(mass)
        time_val = float(rows[0][3])
        
        # 1. Compute bound subsystem virial
        bound_indices, T_bound, V_bound = get_bound_subsystem(
            pos, vel, mass, args.g_constant, args.softening
        )
        if V_bound != 0:
            Q_bound = T_bound / abs(V_bound)
        else:
            Q_bound = np.nan
            
        # 2. Compute entire cluster virial (all stars)
        total_m_all = np.sum(mass)
        cm_pos_all = np.sum(pos * mass[:, np.newaxis], axis=0) / total_m_all
        cm_vel_all = np.sum(vel * mass[:, np.newaxis], axis=0) / total_m_all
        rel_pos_all = pos - cm_pos_all
        rel_vel_all = vel - cm_vel_all
        
        T_all = 0.5 * np.sum(mass * np.sum(rel_vel_all**2, axis=1))
        
        pos_sq_all = np.sum(rel_pos_all**2, axis=1)
        dist_sq_all = pos_sq_all[:, np.newaxis] + pos_sq_all[np.newaxis, :] - 2 * np.dot(rel_pos_all, rel_pos_all.T)
        dist_sq_all = np.maximum(dist_sq_all, 0.0)
        inv_dist_all = 1.0 / np.sqrt(dist_sq_all + args.softening**2)
        np.fill_diagonal(inv_dist_all, 0)
        
        Phi_all = -args.g_constant * np.sum(mass * inv_dist_all, axis=1)
        V_all = 0.5 * np.sum(mass * Phi_all)
        
        if V_all != 0:
            Q_all = T_all / abs(V_all)
        else:
            Q_all = np.nan
            
        bound_frac = len(bound_indices) / len(mass)
        
        times.append(time_val)
        virial_ratios.append(Q_bound)
        virial_ratios_all.append(Q_all)
        bound_fractions.append(bound_frac)

    print("\nAnalysis complete. Generating plots...")
    cursor.close()
    conn.close()

    # Create figure with 2 subplots (Virial Ratio and Bound Fraction)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Plot Virial Ratio
    ax1.plot(times, virial_ratios, color='crimson', label='Bound Core Virial Q = T/|V|', linewidth=2)
    ax1.plot(times, virial_ratios_all, color='darkorange', linestyle='--', label='Entire Cluster Virial Q = T/|V|', linewidth=2)
    ax1.axhline(y=0.5, color='gray', linestyle=':', label='Virial Equilibrium (Q=0.5)')
    ax1.set_ylabel('Virial Ratio Q')
    ax1.set_title(f'Evolution of the Virial Ratio (N={dim_space} dimensions)')
    ax1.legend()
    ax1.grid(True, linestyle='--', alpha=0.5)

    # Plot Bound Fraction
    ax2.plot(times, np.array(bound_fractions) * 100, color='royalblue', label='Bound Star %', linewidth=2)
    ax2.set_xlabel('Time (Myr)')
    ax2.set_ylabel('Bound Stars (%)')
    ax2.legend()
    ax2.grid(True, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plot_filename = f'virial_evolution_{dim_space}D.png'
    plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
    print(f"Plot successfully saved to '{plot_filename}'!")

if __name__ == "__main__":
    main()
