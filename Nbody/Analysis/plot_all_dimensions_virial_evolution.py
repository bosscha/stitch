#!/usr/bin/env python3
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt
import psycopg2

def compute_virial_for_subset(pos, vel, mass, G=1.0, softening=0.001, dim=3, max_particles=2500):
    """
    Computes kinetic energy T, potential energy V, and virial ratio Q = T/|V|
    for a given particle subset using vectorized float32 matrix math.
    Subsamples particles to max_particles for ultrafast potential computation.
    """
    N = len(mass)
    if N < 3:
        return np.nan, 0.0, 0.0

    if N > max_particles:
        idx = np.random.choice(N, max_particles, replace=False)
        pos = pos[idx]
        vel = vel[idx]
        mass = mass[idx]
        N = max_particles

    pos = pos.astype(np.float32)
    vel = vel.astype(np.float32)
    mass = (mass / np.sum(mass)).astype(np.float32)

    total_m = np.sum(mass)
    cm_pos = np.sum(pos * mass[:, np.newaxis], axis=0) / total_m
    cm_vel = np.sum(vel * mass[:, np.newaxis], axis=0) / total_m

    rel_pos = pos - cm_pos
    rel_vel = vel - cm_vel

    # Kinetic Energy T
    T = float(0.5 * np.sum(mass * np.sum(rel_vel**2, axis=1)))

    # Fast pairwise potential calculation
    pos_sq = np.sum(rel_pos**2, axis=1)
    dist_sq = pos_sq[:, np.newaxis] + pos_sq[np.newaxis, :] - 2.0 * np.dot(rel_pos, rel_pos.T)
    np.maximum(dist_sq, 0.0, out=dist_sq)

    exp = max(0.5, (dim - 1) / 2.0)
    inv_dist_pow = 1.0 / np.power(dist_sq + np.float32(softening**2), exp)
    np.fill_diagonal(inv_dist_pow, 0.0)

    Phi_i = -np.float32(G) * np.sum(mass * inv_dist_pow, axis=1)
    V = float(0.5 * np.sum(mass * Phi_i))

    if abs(V) > 1e-15:
        Q = T / abs(V)
    else:
        Q = np.nan

    return Q, T, V

def main():
    parser = argparse.ArgumentParser(description="Plot and compare the evolution of inner, outer, and total virial ratios across all dimensions.")
    parser.add_argument("--host", type=str, default="localhost", help="Database host")
    parser.add_argument("--user", type=str, default="stephane", help="Database user")
    parser.add_argument("--password", type=str, default="tallis", help="Database password")
    parser.add_argument("--dbname", type=str, default="hypercluster", help="Database name")
    parser.add_argument("--softening", type=float, default=0.001, help="Softening parameter")
    parser.add_argument("--g-constant", type=float, default=1.0, help="Gravitational constant G")
    parser.add_argument("--max-points", type=int, default=40, help="Max snapshots to analyze per dimension for performance")
    parser.add_argument("--output-grid", type=str, default="virial_evolution_all_dims.png", help="Multi-panel plot output filename")
    parser.add_argument("--output-overlay", type=str, default="virial_evolution_overlay_all_dims.png", help="Overlay plot output filename")
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

        cursor.execute("SELECT DISTINCT snapshot_id FROM star_snapshots WHERE dim_space = %s ORDER BY snapshot_id ASC;", (dim,))
        all_snapshots = [r[0] for r in cursor.fetchall()]

        if not all_snapshots:
            print(f"  Skipping {dim}D: No snapshot IDs found.")
            continue

        sample_interval = max(1, len(all_snapshots) // args.max_points)
        sampled_snaps = [int(s) for s in all_snapshots[::sample_interval]]

        print(f"  Analyzing {len(sampled_snaps)} sampled snapshots out of {len(all_snapshots)} total...")

        times = []
        q_total_list = []
        q_inner_list = []
        q_outer_list = []

        for snap_id in sampled_snaps:
            cursor.execute("""
                SELECT position, velocity, mass, time_myr FROM (
                    SELECT id, position, velocity, mass, time_myr 
                    FROM star_snapshots 
                    WHERE dim_space = %s AND snapshot_id = %s 
                    ORDER BY id DESC 
                    LIMIT 10000
                ) as sub ORDER BY id ASC;
            """, (dim, snap_id))

            rows = cursor.fetchall()
            if not rows:
                continue

            pos = np.array([r[0] for r in rows], dtype=np.float32)
            vel = np.array([r[1] for r in rows], dtype=np.float32)
            mass = np.array([r[2] for r in rows], dtype=np.float32)
            mass_norm = mass / np.sum(mass)
            time_myr = float(rows[0][3])

            # Total cluster virial ratio
            q_tot, _, _ = compute_virial_for_subset(pos, vel, mass_norm, args.g_constant, args.softening, dim)

            # Radial distances from barycenter
            cm_pos = np.sum(pos * mass_norm[:, np.newaxis], axis=0) / np.sum(mass_norm)
            dist = np.sqrt(np.sum((pos - cm_pos)**2, axis=1))

            # Split into inner (R <= median R50) and outer (R > R50)
            r50 = np.median(dist)
            inner_mask = dist <= r50
            outer_mask = dist > r50

            q_in, _, _ = compute_virial_for_subset(pos[inner_mask], vel[inner_mask], mass_norm[inner_mask], args.g_constant, args.softening, dim)
            q_out, _, _ = compute_virial_for_subset(pos[outer_mask], vel[outer_mask], mass_norm[outer_mask], args.g_constant, args.softening, dim)

            times.append(time_myr)
            q_total_list.append(q_tot)
            q_inner_list.append(q_in)
            q_outer_list.append(q_out)

        dim_data[dim] = {
            "times": np.array(times),
            "q_total": np.array(q_total_list),
            "q_inner": np.array(q_inner_list),
            "q_outer": np.array(q_outer_list)
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

    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows), sharex=False, sharey=False)
    axes_flat = axes.flatten() if n_dims > 1 else [axes]

    for idx, (dim, data) in enumerate(dim_data.items()):
        ax = axes_flat[idx]

        ax.plot(data["times"], data["q_inner"], color='#1f77b4', label='Inner Virial ($R \\le R_{50}$)', linewidth=2)
        ax.plot(data["times"], data["q_outer"], color='#ff7f0e', label='Outer Virial ($R > R_{50}$)', linewidth=2)
        ax.plot(data["times"], data["q_total"], color='#d62728', linestyle='--', label='Total Cluster Virial', linewidth=2)

        ax.axhline(y=0.5, color='black', linestyle=':', label='Equilibrium ($Q = 0.5$)', alpha=0.7)

        ax.set_title(f'Dimension N = {dim}D', fontsize=12, fontweight='bold')
        ax.set_xlabel('Time (Myr)', fontsize=10)
        ax.set_ylabel('Virial Ratio $Q = T / |V|$', fontsize=10)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.legend(fontsize=8, loc='best')

    for j in range(idx + 1, len(axes_flat)):
        fig.delaxes(axes_flat[j])

    plt.suptitle('Evolution of Inner, Outer, and Total Virial Ratio Q(t) Across Dimensions', fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.savefig(args.output_grid, dpi=300, bbox_inches='tight')
    print(f"\nMulti-panel virial evolution plot saved to '{args.output_grid}'!")

    # 2. Overlay Plot for Total & Inner Virial Ratios
    plt.figure(figsize=(10, 7))
    colors = plt.cm.plasma(np.linspace(0.1, 0.9, len(dim_data)))

    plt.axhline(y=0.5, color='black', linestyle=':', label='Virial Equilibrium ($Q = 0.5$)', linewidth=1.5)

    for (dim, data), col in zip(dim_data.items(), colors):
        plt.plot(data["times"], data["q_inner"], linestyle='-', color=col, label=f'N = {dim}D Inner', linewidth=2)
        plt.plot(data["times"], data["q_total"], linestyle='--', color=col, label=f'N = {dim}D Total', linewidth=1.5, alpha=0.7)

    plt.xlabel('Time (Myr)', fontsize=12)
    plt.ylabel('Virial Ratio $Q = T / |V|$', fontsize=12)
    plt.title('Comparison of Inner and Total Virial Ratio Evolution Q(t) Across Dimensions', fontsize=13, fontweight='bold')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend(fontsize=9, loc='best')

    plt.savefig(args.output_overlay, dpi=300, bbox_inches='tight')
    print(f"Overlay virial evolution plot saved to '{args.output_overlay}'!")

if __name__ == "__main__":
    main()
