struct Particle {
    pos_mass: vec4<f32>,
    vel: vec4<f32>,
    acc: vec4<f32>,
}

struct SimParams {
    num_particles: u32,
    g_const: f32,
    dt: f32,
    softening_sq: f32,
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
    
    p.pos_mass.x += p.vel.x * dt + 0.5 * p.acc.x * dt * dt;
    p.pos_mass.y += p.vel.y * dt + 0.5 * p.acc.y * dt * dt;
    p.pos_mass.z += p.vel.z * dt + 0.5 * p.acc.z * dt * dt;
    
    p.vel.x += 0.5 * p.acc.x * dt;
    p.vel.y += 0.5 * p.acc.y * dt;
    p.vel.z += 0.5 * p.acc.z * dt;
    
    particles[index] = p;
}

@compute @workgroup_size(256)
fn update_vel(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= params.num_particles) {
        return;
    }
    
    let p_target = particles[index].pos_mass;
    var new_acc = vec3<f32>(0.0);
    
    for (var j: u32 = 0u; j < params.num_particles; j = j + 1u) {
        let p_source = particles[j].pos_mass;
        
        let dx = p_source.x - p_target.x;
        let dy = p_source.y - p_target.y;
        let dz = p_source.z - p_target.z;
        
        let dist_sq = dx*dx + dy*dy + dz*dz;
        
        let inv_dist_sq = dist_sq + params.softening_sq;
        let inv_dist = 1.0 / sqrt(inv_dist_sq);
        let inv_dist_cube = inv_dist * inv_dist * inv_dist;
        
        let force_mag = params.g_const * p_source.w * inv_dist_cube;
        new_acc += vec3<f32>(dx, dy, dz) * force_mag;
    }
    
    let dt = params.dt;
    particles[index].acc = vec4<f32>(new_acc, 0.0);
    particles[index].vel.x += 0.5 * new_acc.x * dt;
    particles[index].vel.y += 0.5 * new_acc.y * dt;
    particles[index].vel.z += 0.5 * new_acc.z * dt;
}
