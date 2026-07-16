struct Particle {
    pos: array<f32, 8>,
    vel: array<f32, 8>,
    acc: array<f32, 8>,
    mass_pad: array<f32, 8>,
}

struct SimParams {
    num_particles: u32,
    g_const: f32,
    dt: f32,
    softening_sq: f32,
    dim: u32,
    pad1: u32,
    pad2: u32,
    pad3: u32,
}

@group(0) @binding(0) var<storage, read_write> particles: array<Particle>;
@group(0) @binding(1) var<uniform> params: SimParams;

@compute @workgroup_size(256)
fn update_pos(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= params.num_particles) {
        return;
    }
    
    var p = particles[index];
    let dt = params.dt;
    let dim = params.dim;
    
    for (var k: u32 = 0u; k < dim; k = k + 1u) {
        p.pos[k] += p.vel[k] * dt + 0.5 * p.acc[k] * dt * dt;
        p.vel[k] += 0.5 * p.acc[k] * dt;
    }
    
    particles[index] = p;
}

@compute @workgroup_size(256)
fn update_vel(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= params.num_particles) {
        return;
    }
    
    let dim = params.dim;
    var p_target_pos = array<f32, 8>();
    for (var k: u32 = 0u; k < dim; k = k + 1u) {
        p_target_pos[k] = particles[index].pos[k];
    }
    
    var new_acc = array<f32, 8>();
    for (var k: u32 = 0u; k < 8u; k = k + 1u) {
        new_acc[k] = 0.0;
    }
    
    for (var j: u32 = 0u; j < params.num_particles; j = j + 1u) {
        let mass = particles[j].mass_pad[0];
        
        var dist_sq = 0.0;
        var diffs = array<f32, 8>();
        
        for (var k: u32 = 0u; k < dim; k = k + 1u) {
            let diff = particles[j].pos[k] - p_target_pos[k];
            diffs[k] = diff;
            dist_sq += diff * diff;
        }
        
        let inv_dist_sq = 1.0 / (dist_sq + params.softening_sq);
        let inv_dist = sqrt(inv_dist_sq);
        
        // General formula for N-dimensional gravity: force_mag = G * m * inv_dist^(N)
        let force_mag = params.g_const * mass * pow(inv_dist, f32(dim));
        
        for (var k: u32 = 0u; k < dim; k = k + 1u) {
            new_acc[k] += diffs[k] * force_mag;
        }
    }
    
    let dt = params.dt;
    for (var k: u32 = 0u; k < dim; k = k + 1u) {
        particles[index].acc[k] = new_acc[k];
        particles[index].vel[k] += 0.5 * new_acc[k] * dt;
    }
}
