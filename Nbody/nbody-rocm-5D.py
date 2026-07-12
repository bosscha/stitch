import torch
import numpy as np
import time

# NBody simulation for GPU AMD 8060S 
# python env: astro_env

# 1. Simulation Constants
N_STARS = 1000
G = 1.0          # Astrophysical units
DT = 0.001        # Time step
STEPS = 50000      # Number of integration steps
SOFTENING = 0.05 # Prevent numerical infinities during close encounters

# Salpeter IMF Constants
M_MIN = 0.1      # Minimum stellar mass (solar masses)
M_MAX = 50.0     # Maximum stellar mass (solar masses)
ALPHA = 2.35     # Salpeter exponent

# 2. Initialize Positions & Velocities on CPU
np.random.seed(42)
pos = np.random.randn(N_STARS, 3).astype(np.float32) * 5.0  # Distributed in a sphere
vel = np.random.randn(N_STARS, 3).astype(np.float32) * 0.5  # Initial orbital speed

# 3. Push Base Arrays to the AMD 8060S
pos_gpu = torch.tensor(pos, device="cuda")
vel_gpu = torch.tensor(vel, device="cuda")

# 4. Generate Salpeter IMF directly on the GPU
# Using Inverse Transform Sampling for a power law: 
# M = [U * (M_max^(1-α) - M_min^(1-α)) + M_min^(1-α)] ^ (1 / (1-α))
with torch.no_grad():
    # Uniform random numbers between 0 and 1 on the GPU
    U = torch.rand((N_STARS, 1), device="cuda", dtype=torch.float32)
    
    pow_idx = 1.0 - ALPHA
    m_diff = (M_MAX ** pow_idx) - (M_MIN ** pow_idx)
    
    # Sample individual physical masses
    mass_gpu = (U * m_diff + (M_MIN ** pow_idx)).pow(1.0 / pow_idx)
    
    # Normalize the cluster's total mass to 1.0 for numerical stability
    total_mass = torch.sum(mass_gpu)
    mass_gpu = mass_gpu / total_mass

# Calculate Time conversion factor to Myr
# G = 0.00449 pc^3 / (M_sun * Myr^2)
# Time_unit = sqrt( L_unit^3 / (G_phys * M_unit) )
G_PHYS = 0.00449
L_UNIT = 1.0 # Assuming positions are in parsecs
time_to_myr = np.sqrt((L_UNIT**3) / (G_PHYS * total_mass.item()))

print(f"IMF Generated on GPU.")
print(f"➔ Total initial mass: {total_mass.item():.2f} M_sun")
print(f"➔ Heaviest Star relative mass: {torch.max(mass_gpu).item():.5f}")
print(f"➔ Lightest Star relative mass: {torch.min(mass_gpu).item():.5f}")
print(f"➔ Time Conversion: 1 N-body time unit = {time_to_myr:.4f} Myr")

# 5. Core N-Body Gravitational Kernel
def compute_gravitational_accelerations(pos, mass):
    # dx[i, j] = pos[j] - pos[i] -> Shape: (N, N, 3)
    diff = pos.unsqueeze(0) - pos.unsqueeze(1) 
    
    # Pairwise squared distances -> Shape: (N, N)
    dist_sq = torch.sum(diff ** 2, dim=-1)
    
    # Softened inverse cube distance term -> Shape: (N, N)
    inv_dist_cubed = (dist_sq + SOFTENING**2).pow(-1.5)
    inv_dist_cubed.fill_diagonal_(0.0) # Self-gravity suppression
    
    # Physics adjustment: Force depends on the mass of the pulling star (mass of star j)
    # mass.unsqueeze(-1) matches the (N, N, 3) broadcast requirements
    # mass behaves as shape (1, N, 1) to act as the source mass m_j pulling on target i
    forces = diff * inv_dist_cubed.unsqueeze(-1) * mass.view(1, -1, 1)
    
    # Sum up all forces acting on each star i
    acc = G * torch.sum(forces, dim=1)
    return acc

# 6. Main Integration Loop (Velocity Verlet)
print(f"\nStarting Salpeter N-body simulation on {torch.cuda.get_device_name(0)}...")
start_time = time.time()

acc_gpu = compute_gravitational_accelerations(pos_gpu, mass_gpu)

for step in range(STEPS):
    # Update positions
    pos_gpu += vel_gpu * DT + 0.5 * acc_gpu * (DT ** 2)
    
    # Save half-step velocity
    vel_half = vel_gpu + 0.5 * acc_gpu * DT
    
    # Compute new acceleration with updated positions
    acc_gpu = compute_gravitational_accelerations(pos_gpu, mass_gpu)
    
    # Update velocities
    vel_gpu = vel_half + 0.5 * acc_gpu * DT
    
    if step % 1000 == 0:
        time_myr = step * DT * time_to_myr
        print(f"Step {step:04d}/{STEPS} completed. Physical Time: {time_myr:.2f} Myr")

torch.cuda.synchronize()
print(f"\nSimulation completed in: {time.time() - start_time:.4f} seconds!")
print(f"Total physical time simulated: {STEPS * DT * time_to_myr:.2f} Myr")

# Final positions and masses ready for plotting or analysis
final_positions = pos_gpu.cpu().numpy()
final_masses = mass_gpu.cpu().numpy()
