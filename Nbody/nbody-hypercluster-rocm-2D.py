import torch
import numpy as np
import time
import math
import matplotlib.pyplot as plt
import psycopg2
from psycopg2 import extras

# NBody simulation for GPU AMD 8060S 
# python env: astro_env

# 1. Simulation Constants
N_STARS = 5000
DIM = 2          # Number of spatial dimensions
G = 1.0          # Astrophysical units
DT = 0.0001        # Time step
STEPS = 5000000      # Number of integration steps
SOFTENING = 0.01 # Prevent numerical infinities during close encounters

# Salpeter IMF Constants
M_MIN = 0.1      # Minimum stellar mass (solar masses)
M_MAX = 50.0     # Maximum stellar mass (solar masses)
ALPHA = 2.35     # Salpeter exponent

# 2. Initialize Positions & Velocities on CPU
np.random.seed(42)

def sample_king_nd(n_stars, dim, r_c=1.0, r_t=10.0):
    """
    Sample an N-dimensional King-like distribution.
    Uses rejection sampling for the radial coordinate based on the generalized King profile density.
    """
    r_grid = np.linspace(0, r_t, 1000)
    p_grid = r_grid**(dim-1) * (1/np.sqrt(1 + (r_grid/r_c)**2) - 1/np.sqrt(1 + (r_t/r_c)**2))**2
    p_max = np.max(p_grid) * 1.1
    
    r_samples = []
    while len(r_samples) < n_stars:
        r_cand = np.random.uniform(0, r_t, n_stars)
        p_cand = np.random.uniform(0, p_max, n_stars)
        p_eval = r_cand**(dim-1) * (1/np.sqrt(1 + (r_cand/r_c)**2) - 1/np.sqrt(1 + (r_t/r_c)**2))**2
        accepted = r_cand[p_cand < p_eval]
        r_samples.extend(accepted.tolist())
        
    r_samples = np.array(r_samples[:n_stars])
    
    v = np.random.randn(n_stars, dim)
    v_norm = np.linalg.norm(v, axis=1, keepdims=True)
    unit_v = v / v_norm
    
    pos = unit_v * r_samples[:, np.newaxis]
    return pos.astype(np.float32)

pos = sample_king_nd(N_STARS, DIM, r_c=1.0, r_t=10.0)
vel = np.random.randn(N_STARS, DIM).astype(np.float32) * 0.5  # Initial orbital speed

# 3. Push Base Arrays to the AMD 8060S
pos_gpu = torch.tensor(pos, device="cuda")
vel_gpu = torch.tensor(vel, device="cuda")

# 4. Generate Salpeter IMF directly on the GPU
with torch.no_grad():
    U = torch.rand((N_STARS, 1), device="cuda", dtype=torch.float32)
    pow_idx = 1.0 - ALPHA
    m_diff = (M_MAX ** pow_idx) - (M_MIN ** pow_idx)
    mass_gpu = (U * m_diff + (M_MIN ** pow_idx)).pow(1.0 / pow_idx)
    total_mass = torch.sum(mass_gpu)
    mass_gpu = mass_gpu / total_mass

    # Shift velocities to stop CoM drift
    mass_sum = mass_gpu.sum()
    vel_cm = (vel_gpu * mass_gpu).sum(dim=0) / mass_sum
    vel_gpu -= vel_cm

    # Compute energies and scale velocities to Virial Equilibrium (Q = 0.5)
    T = 0.5 * torch.sum(mass_gpu * torch.sum(vel_gpu**2, dim=1, keepdim=True)).item()
    dist_mat = torch.cdist(pos_gpu, pos_gpu, p=2)
    inv_dist = 1.0 / torch.sqrt(dist_mat**2 + SOFTENING**2)
    inv_dist.fill_diagonal_(0)
    V = -0.5 * G * torch.sum(mass_gpu * torch.matmul(inv_dist, mass_gpu)).item()
    
    f_scale = math.sqrt(0.5 * abs(V) / T)
    vel_gpu *= f_scale
    print(f"➔ Initial Virial Ratio before scaling: {T/abs(V):.4f}")
    print(f"➔ Scaled velocities by factor {f_scale:.4f} to achieve virial equilibrium (Q = 0.5)")

G_PHYS = 0.00449
L_UNIT = 1.0
time_to_myr = np.sqrt((L_UNIT**3) / (G_PHYS * total_mass.item()))

print(f"IMF Generated on GPU.")
print(f"➔ Total initial mass: {total_mass.item():.2f} M_sun")
print(f"➔ Heaviest Star relative mass: {torch.max(mass_gpu).item():.5f}")
print(f"➔ Lightest Star relative mass: {torch.min(mass_gpu).item():.5f}")
print(f"➔ Time Conversion: 1 N-body time unit = {time_to_myr:.4f} Myr")

# 5. Core N-Body Gravitational Kernel
def compute_gravitational_accelerations(pos, mass):
    acc = torch.zeros_like(pos)
    chunk_size = 1000  # Process in batches to avoid GPU Hang / TDR timeout on AMD ROCm
    
    # Process target stars in chunks
    for i in range(0, pos.size(0), chunk_size):
        end = min(i + chunk_size, pos.size(0))
        
        # diff shape: (chunk_size, N_STARS, DIM)
        diff = pos[i:end].unsqueeze(1) - pos.unsqueeze(0)
        
        dist_sq = torch.sum(diff ** 2, dim=-1)
        
        # Softened inverse distance calculation
        inv_dist_power = (dist_sq + SOFTENING**2).pow(-DIM / 2.0)
        
        # Calculate forces. Self-interaction evaluates to zero since diff is zero.
        forces = diff * inv_dist_power.unsqueeze(-1) * mass.view(1, -1, 1)
        acc[i:end] = G * torch.sum(forces, dim=1)
        
    return acc

# 6. Main Integration Loop (Velocity Verlet)
print(f"\nStarting Salpeter N-body simulation on {torch.cuda.get_device_name(0)}...")
start_time = time.time()

# Calculate initial 80% radius
center_initial = pos_gpu.mean(dim=0)
dist_initial = torch.sqrt(torch.sum((pos_gpu - center_initial)**2, dim=1))
r80_initial = torch.quantile(dist_initial, 0.8).item()
print(f"➔ Initial 80% radius (R80): {r80_initial:.4f}")

TRACK_INTERVAL = 100
escaped_counts = []
time_steps_list = []
slopes_inside_50 = []
slopes_outside_50 = []

mass_phys_cpu = (mass_gpu * total_mass).cpu().numpy().flatten()

# Database setup
try:
    db_conn = psycopg2.connect(dbname='hypercluster', user='stephane', password='tallis', host='localhost')
    db_cursor = db_conn.cursor()
    db_cursor.execute("""
        CREATE TABLE IF NOT EXISTS star_snapshots (
            id SERIAL PRIMARY KEY,
            snapshot_id INTEGER,
            time_myr DOUBLE PRECISION,
            star_id INTEGER,
            mass DOUBLE PRECISION,
            dim_space INTEGER,
            position DOUBLE PRECISION[],
            velocity DOUBLE PRECISION[]
        );
    """)
    db_conn.commit()
except Exception as e:
    print(f"Failed to connect to database or create table: {e}")
    db_conn = None

def compute_mass_slope(masses, m_min, m_max, num_bins=15):
    if len(masses) < 10: return np.nan
    bins = np.logspace(np.log10(m_min), np.log10(m_max), num_bins)
    counts, _ = np.histogram(masses, bins=bins)
    dn_dm = counts / (bins[1:] - bins[:-1])
    bin_centers = np.sqrt(bins[1:] * bins[:-1])
    valid = dn_dm > 0
    if np.sum(valid) < 3: return np.nan
    slope, _ = np.polyfit(np.log10(bin_centers[valid]), np.log10(dn_dm[valid]), 1)
    return -slope

acc_gpu = compute_gravitational_accelerations(pos_gpu, mass_gpu)

for step in range(STEPS):
    pos_gpu += vel_gpu * DT + 0.5 * acc_gpu * (DT ** 2)
    vel_half = vel_gpu + 0.5 * acc_gpu * DT
    acc_gpu = compute_gravitational_accelerations(pos_gpu, mass_gpu)
    vel_gpu = vel_half + 0.5 * acc_gpu * DT
    
    if step % TRACK_INTERVAL == 0:
        center = pos_gpu.mean(dim=0)
        dists = torch.sqrt(torch.sum((pos_gpu - center)**2, dim=1))
        escaped = torch.sum(dists > r80_initial).item()
        escaped_counts.append(escaped)
        time_steps_list.append(step * DT * time_to_myr)
        
        dists_cpu = dists.cpu().numpy()
        inside_mask = dists_cpu <= 50.0
        outside_mask = dists_cpu > 50.0
        slopes_inside_50.append(compute_mass_slope(mass_phys_cpu[inside_mask], M_MIN, M_MAX))
        slopes_outside_50.append(compute_mass_slope(mass_phys_cpu[outside_mask], M_MIN, M_MAX))
        
        if db_conn is not None:
            pos_cpu = pos_gpu.cpu().numpy()
            vel_cpu = vel_gpu.cpu().numpy()
            time_val = step * DT * time_to_myr
            
            # Prepare data for bulk insert
            records = []
            for i in range(N_STARS):
                records.append((
                    step,
                    float(time_val),
                    i,
                    float(mass_phys_cpu[i]),
                    DIM,
                    pos_cpu[i].tolist(),
                    vel_cpu[i].tolist()
                ))
            
            insert_query = """
                INSERT INTO star_snapshots 
                (snapshot_id, time_myr, star_id, mass, dim_space, position, velocity) 
                VALUES %s
            """
            try:
                extras.execute_values(db_cursor, insert_query, records)
                db_conn.commit()
            except Exception as e:
                print(f"Database insertion failed at step {step}: {e}")
                db_conn.rollback()

    if step % 1000 == 0:
        time_myr = step * DT * time_to_myr
        print(f"Step {step:04d}/{STEPS} completed. Physical Time: {time_myr:.2f} Myr")

if db_conn is not None:
    db_cursor.close()
    db_conn.close()

torch.cuda.synchronize()
print(f"\nSimulation completed in: {time.time() - start_time:.4f} seconds!")
print(f"Total physical time simulated: {STEPS * DT * time_to_myr:.2f} Myr")

# Final positions and masses ready for plotting or analysis
final_positions = pos_gpu.cpu().numpy()
final_masses = mass_gpu.cpu().numpy()

# Generate summary plot
plt.figure(figsize=(10, 6))
plt.plot(time_steps_list, escaped_counts, label='Escaped Stars (r > initial R80)', color='b', linewidth=2)
plt.axhline(y=N_STARS * 0.2, color='r', linestyle='--', label='Initial 20% outside R80')
plt.xlabel('Time (Myr)')
plt.ylabel('Number of Escaped Stars')
plt.title(f'Evolution of Escaped Stars (King Profile, N={DIM} dimensions)')
plt.legend()
plt.grid(True)
plot_filename = f'escaped_stars_{DIM}D.png'
plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
print(f"Summary plot of escaped stars saved to '{plot_filename}'")

# Generate radial distribution plot
dist_initial_cpu = dist_initial.cpu().numpy()
center_final = pos_gpu.mean(dim=0)
dist_final_cpu = torch.sqrt(torch.sum((pos_gpu - center_final)**2, dim=1)).cpu().numpy()

plt.figure(figsize=(10, 6))
max_dist = max(np.max(dist_initial_cpu), np.max(dist_final_cpu))
bins = np.linspace(0, max_dist, 100)
plt.hist(dist_initial_cpu, bins=bins, alpha=0.5, label='Initial Distribution', color='b')
plt.hist(dist_final_cpu, bins=bins, alpha=0.5, label='Final Distribution', color='r')
plt.yscale('log')
plt.xlabel('Radial Distance from Center')
plt.ylabel('Number of Stars (log scale)')
plt.title(f'Initial vs Final Radial Distribution (N={DIM} dimensions)')
plt.legend()
plt.grid(True)
plot_rad_filename = f'radial_distribution_{DIM}D.png'
plt.savefig(plot_rad_filename, dpi=300, bbox_inches='tight')
print(f"Radial distribution plot saved to '{plot_rad_filename}'")

# Generate radial density plot
plt.figure(figsize=(10, 6))

# Use log-spaced bins for a better representation in log-log scale
min_dist = max(np.min(dist_initial_cpu[dist_initial_cpu > 0]) if np.any(dist_initial_cpu > 0) else 1e-2, 1e-3)
bins_log = np.logspace(np.log10(min_dist), np.log10(max_dist), 50)

counts_init, bins_edges = np.histogram(dist_initial_cpu, bins=bins_log)
counts_final, _ = np.histogram(dist_final_cpu, bins=bins_edges)

def volume_nd(r, dim):
    return (math.pi**(dim/2.0) / math.gamma(dim/2.0 + 1.0)) * (r**dim)

volumes = volume_nd(bins_edges[1:], DIM) - volume_nd(bins_edges[:-1], DIM)
volumes = np.maximum(volumes, 1e-15)

density_init = counts_init / volumes
density_final = counts_final / volumes
bin_centers = (bins_edges[:-1] + bins_edges[1:]) / 2.0

# Plot density; only plot where density > 0 to avoid log(0) issues in plot lines
valid_init = density_init > 0
valid_final = density_final > 0

plt.plot(bin_centers[valid_init], density_init[valid_init], label='Initial Density', color='b', marker='.', linestyle='-')
plt.plot(bin_centers[valid_final], density_final[valid_final], label='Final Density', color='r', marker='.', linestyle='-')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Radial Distance from Center')
plt.ylabel('Density (Stars / Volume)')
plt.title(f'Initial vs Final Radial Density (N={DIM} dimensions)')
plt.legend()
plt.grid(True)
plot_dens_filename = f'radial_density_{DIM}D.png'
plt.savefig(plot_dens_filename, dpi=300, bbox_inches='tight')
print(f"Radial density plot saved to '{plot_dens_filename}'")

# Generate mass slope evolution plot
plt.figure(figsize=(10, 6))
plt.plot(time_steps_list, slopes_inside_50, label='Inside 50 pc', color='purple', linewidth=2)
plt.plot(time_steps_list, slopes_outside_50, label='Outside 50 pc', color='orange', linewidth=2)
plt.axhline(y=ALPHA, color='k', linestyle='--', label=f'Initial Salpeter (α={ALPHA})')
plt.xlabel('Time (Myr)')
plt.ylabel('Mass Slope α (from dN/dM ∝ M^(-α))')
plt.title('Evolution of Mass Function Slope Inside/Outside 50 pc')
plt.legend()
plt.grid(True)
plot_slope_filename = f'mass_slope_evolution_{DIM}D.png'
plt.savefig(plot_slope_filename, dpi=300, bbox_inches='tight')
print(f"Mass slope evolution plot saved to '{plot_slope_filename}'")
