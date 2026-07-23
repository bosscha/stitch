#!/usr/bin/env python3
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def compute_energies(pos, vel, mass, G=1.0, softening=0.001, dim=3, max_particles=2500):
    """
    Computes Kinetic Energy T, Potential Energy V, and Total Energy E = T + V
    for a given particle subset using vectorized float32 matrix math.
    Subsamples particles to max_particles for performance on larger systems.
    Uses the mathematically correct D-dimensional potential law.
    """
    N = len(mass)
    if N < 3:
        return 0.0, 0.0, 0.0

    if N > max_particles:
        idx = np.random.choice(N, max_particles, replace=False)
        pos = pos[idx]
        vel = vel[idx]
        mass = mass[idx]
        N = max_particles

    pos = pos.astype(np.float32)
    vel = vel.astype(np.float32)
    # Normalize masses to sum to 1.0 for consistent N-body units
    mass = (mass / np.sum(mass)).astype(np.float32)

    total_m = np.sum(mass)
    cm_pos = np.sum(pos * mass[:, np.newaxis], axis=0) / total_m
    cm_vel = np.sum(vel * mass[:, np.newaxis], axis=0) / total_m

    rel_pos = pos - cm_pos
    rel_vel = vel - cm_vel

    # Kinetic Energy T
    T = float(0.5 * np.sum(mass * np.sum(rel_vel**2, axis=1)))

    # Pairwise squared distances
    pos_sq = np.sum(rel_pos**2, axis=1)
    dist_sq = pos_sq[:, np.newaxis] + pos_sq[np.newaxis, :] - 2.0 * np.dot(rel_pos, rel_pos.T)
    np.maximum(dist_sq, 0.0, out=dist_sq)

    # Softened distance squared
    r2_soft = dist_sq + np.float32(softening**2)

    # Compute Potential Energy V depending on dimension
    if dim == 1:
        # V_ij = -G * m_i * m_j * sqrt(r_ij^2 + eps^2)
        r_soft = np.sqrt(r2_soft)
        v_mat = -np.float32(G) * mass[:, np.newaxis] * mass[np.newaxis, :] * r_soft
        np.fill_diagonal(v_mat, 0.0)
        V = float(0.5 * np.sum(v_mat))
    elif dim == 2:
        # V_ij = -0.5 * G * m_i * m_j * ln(r_ij^2 + eps^2)
        ln_r2 = np.log(r2_soft)
        v_mat = -0.5 * np.float32(G) * mass[:, np.newaxis] * mass[np.newaxis, :] * ln_r2
        np.fill_diagonal(v_mat, 0.0)
        V = float(0.5 * np.sum(v_mat))
    else:
        # V_ij = -G * m_i * m_j / ((d - 2) * (r_ij^2 + eps^2)^((d - 2) / 2))
        exp = (dim - 2) / 2.0
        inv_r_power = 1.0 / np.power(r2_soft, exp)
        v_mat = -np.float32(G) * mass[:, np.newaxis] * mass[np.newaxis, :] * inv_r_power / np.float32(dim - 2)
        np.fill_diagonal(v_mat, 0.0)
        V = float(0.5 * np.sum(v_mat))

    E = T + V
    return T, V, E

def main():
    parser = argparse.ArgumentParser(description="Plot and compare the evolution of kinetic, potential, and total energy across all dimensions.")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--softening", type=float, default=0.001, help="Softening parameter")
    parser.add_argument("--g-constant", type=float, default=1.0, help="Gravitational constant G")
    parser.add_argument("--max-points", type=int, default=100, help="Max snapshots to analyze per dimension for performance")
    parser.add_argument("--output-grid", type=str, default="energy_evolution_all_dims.png", help="Multi-panel plot output filename")
    parser.add_argument("--output-conservation", type=str, default="energy_conservation_all_dims.png", help="Relative energy error plot output filename")
    args = parser.parse_args()

    np.random.seed(42)  # Reproducible subsampling

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

    # Use LIMIT 1 dimension check technique to quickly identify active dimensions
    candidate_dims = [1, 2, 3, 4, 5, 6, 25, 50, 100]
    dims = []
    for dim in candidate_dims:
        cursor.execute("SELECT 1 FROM star_snapshots WHERE dim_space = %s LIMIT 1;", (dim,))
        if cursor.fetchone():
            dims.append(dim)

    if not dims:
        print("No simulation data found in table 'star_snapshots'.")
        cursor.close()
        conn.close()
        return

    print(f"Found {len(dims)} spatial dimensions to process: {dims}")

    dim_data = {}

    for dim in dims:
        print(f"\nProcessing dim_space = {dim}D...")

        # 1. Get the min and max snapshot ID
        cursor.execute("SELECT MIN(snapshot_id), MAX(snapshot_id) FROM star_snapshots WHERE dim_space = %s;", (dim,))
        min_snap, max_snap = cursor.fetchone()

        if min_snap is None or max_snap is None:
            print(f"  Skipping {dim}D: No snapshot data.")
            continue

        # 2. Get the number of stars from the latest snapshot
        cursor.execute("SELECT COUNT(*) FROM star_snapshots WHERE dim_space = %s AND snapshot_id = %s;", (dim, max_snap))
        num_stars = cursor.fetchone()[0]

        if num_stars == 0:
            print(f"  Skipping {dim}D: 0 stars found.")
            continue

        # 3. Find the starting database ID for this run (minimum id of the latest min_snapshot_id)
        cursor.execute("""
            SELECT MIN(id) FROM (
                SELECT id FROM star_snapshots 
                WHERE dim_space = %s AND snapshot_id = %s 
                ORDER BY id DESC 
                LIMIT %s
            ) as sub;
        """, (dim, min_snap, num_stars))
        start_db_id = cursor.fetchone()[0]

        # 4. Find the step size dynamically
        cursor.execute("""
            SELECT MIN(snapshot_id) 
            FROM star_snapshots 
            WHERE dim_space = %s AND snapshot_id > %s AND id >= %s;
        """, (dim, min_snap, start_db_id))
        second_snap = cursor.fetchone()[0]
        step_size = (second_snap - min_snap) if second_snap is not None else 100

        # Generate full snapshot IDs list
        snap_ids = list(range(min_snap, max_snap + 1, step_size))
        total_snapshots = len(snap_ids)

        # Downsample snapshots
        sample_interval = max(1, total_snapshots // args.max_points)
        sampled_snaps = snap_ids[::sample_interval]

        print(f"  Analyzing {len(sampled_snaps)} sampled snapshots out of {total_snapshots} total...")

        times = []
        t_list = []
        v_list = []
        e_list = []

        for snap_id in sampled_snaps:
            cursor.execute("""
                SELECT position, velocity, mass, time_myr FROM (
                    SELECT id, position, velocity, mass, time_myr 
                    FROM star_snapshots 
                    WHERE dim_space = %s AND snapshot_id = %s AND id >= %s
                    ORDER BY id DESC 
                    LIMIT %s
                ) as sub ORDER BY id ASC;
            """, (dim, snap_id, start_db_id, num_stars))

            rows = cursor.fetchall()
            if not rows:
                continue

            pos = np.array([r[0] for r in rows], dtype=np.float32)
            vel = np.array([r[1] for r in rows], dtype=np.float32)
            mass = np.array([r[2] for r in rows], dtype=np.float32)
            time_myr = float(rows[0][3])

            T, V, E = compute_energies(pos, vel, mass, args.g_constant, args.softening, dim)

            times.append(time_myr)
            t_list.append(T)
            v_list.append(V)
            e_list.append(E)

        dim_data[dim] = {
            "times": np.array(times),
            "T": np.array(t_list),
            "V": np.array(v_list),
            "E": np.array(e_list)
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

    # Use a professional, clean style for plotting
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.size': 10,
        'axes.labelsize': 11,
        'axes.titlesize': 12,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'figure.titlesize': 14
    })

    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 5 * rows), sharex=False, sharey=False)
    axes_flat = axes.flatten() if n_dims > 1 else [axes]

    # Elegant, premium color palette
    c_kinetic = '#3b82f6'    # Vibrant Blue
    c_potential = '#ef4444'  # Vibrant Red
    c_total = '#10b981'      # Vibrant Emerald Green

    for idx, (dim, data) in enumerate(dim_data.items()):
        ax = axes_flat[idx]

        ax.plot(data["times"], data["T"], color=c_kinetic, label='Kinetic Energy ($T$)', linewidth=2.0)
        ax.plot(data["times"], data["V"], color=c_potential, label='Potential Energy ($V$)', linewidth=2.0)
        ax.plot(data["times"], data["E"], color=c_total, label='Total Energy ($E = T+V$)', linewidth=2.2)

        # Plot a reference dashed line at the initial total energy
        e_init = data["E"][0]
        ax.axhline(y=e_init, color='#047857', linestyle='--', alpha=0.6, label='Initial Energy ($E_0$)', linewidth=1.2)

        ax.set_title(f'Dimension N = {dim}D', fontsize=12, fontweight='bold', pad=10)
        ax.set_xlabel('Time (Myr)', fontsize=10)
        ax.set_ylabel('Energy (N-body units)', fontsize=10)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.legend(fontsize=9, loc='best', framealpha=0.9)

    # Hide extra axes if the grid has empty panels
    for j in range(idx + 1, len(axes_flat)):
        fig.delaxes(axes_flat[j])

    plt.suptitle('Evolution of Kinetic, Potential, and Total Energy Across Dimensions', fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(args.output_grid, dpi=300, bbox_inches='tight')
    print(f"\nMulti-panel energy evolution plot saved to '{args.output_grid}'!")

    # 2. Relative Energy Conservation Plot (Overlay)
    plt.figure(figsize=(10, 6.5))
    
    # We will use plasma colormap to distinguish different dimensions
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(dim_data)))

    for (dim, data), col in zip(dim_data.items(), colors):
        e_init = data["E"][0]
        if abs(e_init) > 1e-10:
            rel_err = (data["E"] - e_init) / abs(e_init)
        else:
            rel_err = data["E"] - e_init
            
        plt.plot(data["times"], rel_err, color=col, label=f'{dim}D Simulation', linewidth=2.0)

    plt.axhline(y=0.0, color='black', linestyle=':', linewidth=1.5, alpha=0.7)
    plt.xlabel('Time (Myr)', fontsize=12)
    plt.ylabel('Relative Energy Error $\\Delta E / |E_0|$', fontsize=12)
    plt.title('Relative Energy Conservation Error $\\Delta E(t) / |E_0|$ Across Dimensions', fontsize=13, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend(fontsize=10, loc='best')

    plt.savefig(args.output_conservation, dpi=300, bbox_inches='tight')
    print(f"Energy conservation comparison plot saved to '{args.output_conservation}'!")

if __name__ == "__main__":
    main()
