import sys

with open("rust_nbody_3D_incoming/src/main.rs", "r") as f:
    code = f.read()

# 1. Update Constants
code = code.replace(
    "const N_STARS: usize = 10000;\nconst DIM: usize = 3;",
    "const N_STARS: usize = 10001;\nconst CLUSTER_STARS: usize = 10000;\nconst DIM: usize = 3;"
)

# 2. Init calls
code = code.replace(
    "let mut pos = sample_king_nd(N_STARS, 1.0, 10.0);",
    "let mut pos = sample_king_nd(CLUSTER_STARS, 1.0, 10.0);"
)
code = code.replace(
    "let mut vel: Vec<Vector3> = (0..N_STARS).map(|_| {",
    "let mut vel: Vec<Vector3> = (0..CLUSTER_STARS).map(|_| {"
)

# 3. Mass generation
code = code.replace(
    "let mut mass = vec![0.0f32; N_STARS];\n    let mut total_mass = 0.0f32;\n    \n    for i in 0..N_STARS {",
    "let mut mass = vec![0.0f32; CLUSTER_STARS];\n    let mut total_mass = 0.0f32;\n    \n    for i in 0..CLUSTER_STARS {"
)
code = code.replace(
    "for i in 0..N_STARS {\n        mass[i] /= total_mass;\n    }",
    "for i in 0..CLUSTER_STARS {\n        mass[i] /= total_mass;\n    }"
)

# 4. CoM drift and T_energy
code = code.replace(
    "for i in 0..N_STARS {\n        vel_cm[0] += vel[i][0] * mass[i];",
    "for i in 0..CLUSTER_STARS {\n        vel_cm[0] += vel[i][0] * mass[i];"
)

energy_str = """    let mut t_energy = 0.0f32;
    for i in 0..N_STARS {
        vel[i][0] -= vel_cm[0];
        vel[i][1] -= vel_cm[1];
        vel[i][2] -= vel_cm[2];
        t_energy += 0.5 * mass[i] * (vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2]);
    }

    // Compute potential energy
    let mut v_energy = 0.0f32;
    let mut initial_acc = vec![[0.0f32; 3]; N_STARS];
    for i in 0..N_STARS {
        for j in 0..N_STARS {
            if i == j { continue; }
            let dx = pos[j][0] - pos[i][0];
            let dy = pos[j][1] - pos[i][1];
            let dz = pos[j][2] - pos[i][2];
            let dist_sq = dx*dx + dy*dy + dz*dz;
            let inv_dist_sq = dist_sq + SOFTENING * SOFTENING;
            let inv_dist = 1.0 / inv_dist_sq.sqrt();
            let inv_dist_cube = inv_dist * inv_dist * inv_dist;
            
            if j > i {
                v_energy += -G * mass[i] * mass[j] * inv_dist;
            }
            
            let force_mag = G * mass[j] * inv_dist_cube;
            initial_acc[i][0] += force_mag * dx;
            initial_acc[i][1] += force_mag * dy;
            initial_acc[i][2] += force_mag * dz;
        }
    }
    
    let f_scale = (0.5 * v_energy.abs() / t_energy).sqrt();
    for i in 0..N_STARS {"""

energy_rep = """    let mut t_energy = 0.0f32;
    for i in 0..CLUSTER_STARS {
        vel[i][0] -= vel_cm[0];
        vel[i][1] -= vel_cm[1];
        vel[i][2] -= vel_cm[2];
        t_energy += 0.5 * mass[i] * (vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2]);
    }

    // Compute potential energy
    let mut v_energy = 0.0f32;
    for i in 0..CLUSTER_STARS {
        for j in (i + 1)..CLUSTER_STARS {
            let dx = pos[j][0] - pos[i][0];
            let dy = pos[j][1] - pos[i][1];
            let dz = pos[j][2] - pos[i][2];
            let dist_sq = dx*dx + dy*dy + dz*dz;
            let inv_dist = 1.0 / (dist_sq + SOFTENING * SOFTENING).sqrt();
            v_energy += -G * mass[i] * mass[j] * inv_dist;
        }
    }
    
    let f_scale = (0.5 * v_energy.abs() / t_energy).sqrt();
    for i in 0..CLUSTER_STARS {"""
code = code.replace(energy_str, energy_rep)

# 5. Inject object and initial_acc
inject_str = """    println!("➔ Time Conversion: 1 N-body time unit = {:.4} Myr", time_to_myr);

    let mut db_client = Client::connect("host=localhost user=stephane password=tallis dbname=hypercluster", NoTls).ok();"""
    
inject_rep = """    println!("➔ Time Conversion: 1 N-body time unit = {:.4} Myr", time_to_myr);

    // Inject incoming object
    let incoming_mass = 300.0 / total_mass;
    let v_sim = 204.54 * time_to_myr;
    mass.push(incoming_mass);
    pos.push([200.0, 0.0, 0.0]);
    vel.push([-v_sim, 0.0, 0.0]);

    // Compute initial_acc for all N_STARS
    let mut initial_acc = vec![[0.0f32; 3]; N_STARS];
    for i in 0..N_STARS {
        for j in 0..N_STARS {
            if i == j { continue; }
            let dx = pos[j][0] - pos[i][0];
            let dy = pos[j][1] - pos[i][1];
            let dz = pos[j][2] - pos[i][2];
            let dist_sq = dx*dx + dy*dy + dz*dz;
            let inv_dist_cube = 1.0 / (dist_sq + SOFTENING * SOFTENING).powf(1.5);
            let force_mag = G * mass[j] * inv_dist_cube;
            initial_acc[i][0] += force_mag * dx;
            initial_acc[i][1] += force_mag * dy;
            initial_acc[i][2] += force_mag * dz;
        }
    }

    let mut db_client = Client::connect("host=localhost user=stephane password=tallis dbname=hypercluster", NoTls).ok();"""
code = code.replace(inject_str, inject_rep)

# 6. R80 and Center
r80_str = """    // Calculate initial 80% radius
    let mut center_initial = [0.0f32; 3];
    for p in &pos {
        center_initial[0] += p[0] / N_STARS as f32;
        center_initial[1] += p[1] / N_STARS as f32;
        center_initial[2] += p[2] / N_STARS as f32;
    }
    
    let mut dist_initial: Vec<f32> = pos.iter().map(|p| {
        let dx = p[0] - center_initial[0];
        let dy = p[1] - center_initial[1];
        let dz = p[2] - center_initial[2];
        (dx*dx + dy*dy + dz*dz).sqrt()
    }).collect();
    
    dist_initial.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let r80_initial = dist_initial[(0.8 * N_STARS as f32) as usize];"""

r80_rep = """    // Calculate initial 80% radius (for cluster only)
    let mut center_initial = [0.0f32; 3];
    for i in 0..CLUSTER_STARS {
        center_initial[0] += pos[i][0] / CLUSTER_STARS as f32;
        center_initial[1] += pos[i][1] / CLUSTER_STARS as f32;
        center_initial[2] += pos[i][2] / CLUSTER_STARS as f32;
    }
    
    let mut dist_initial: Vec<f32> = pos[0..CLUSTER_STARS].iter().map(|p| {
        let dx = p[0] - center_initial[0];
        let dy = p[1] - center_initial[1];
        let dz = p[2] - center_initial[2];
        (dx*dx + dy*dy + dz*dz).sqrt()
    }).collect();
    
    dist_initial.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let r80_initial = dist_initial[(0.8 * CLUSTER_STARS as f32) as usize];"""
code = code.replace(r80_str, r80_rep)

# Escaped check center
esc_center_str = """            let mut center = [0.0f32; 3];
            for p in result {
                center[0] += p.pos_mass[0] / N_STARS as f32;
                center[1] += p.pos_mass[1] / N_STARS as f32;
                center[2] += p.pos_mass[2] / N_STARS as f32;
            }"""
esc_center_rep = """            let mut center = [0.0f32; 3];
            for i in 0..CLUSTER_STARS {
                center[0] += result[i].pos_mass[0] / CLUSTER_STARS as f32;
                center[1] += result[i].pos_mass[1] / CLUSTER_STARS as f32;
                center[2] += result[i].pos_mass[2] / CLUSTER_STARS as f32;
            }"""
code = code.replace(esc_center_str, esc_center_rep)

with open("rust_nbody_3D_incoming/src/main.rs", "w") as f:
    f.write(code)

print("Modification complete.")
