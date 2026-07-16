struct Particle {
    pos: vec4<f32>,
    vel: vec4<f32>,
    acc: vec4<f32>,
    mass_pad: vec4<f32>,
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
    
    p.pos += p.vel * dt + 0.5 * p.acc * dt * dt;
    p.vel += 0.5 * p.acc * dt;
    
    particles[index] = p;
}

@compute @workgroup_size(256)
fn update_vel(@builtin(global_invocation_id) global_id: vec3<u32>) {
    let index = global_id.x;
    if (index >= params.num_particles) {
        return;
    }
    
    let p_target = particles[index].pos;
    var new_acc = vec4<f32>(0.0);
    
    for (var j: u32 = 0u; j < params.num_particles; j = j + 1u) {
        let p_source = particles[j].pos;
        let mass = particles[j].mass_pad.x;
        
        let diff = p_source - p_target;
        let dist_sq = dot(diff, diff);
        
        let inv_dist_sq = 1.0 / (dist_sq + params.softening_sq);
        let force_mag = params.g_const * mass * inv_dist_sq * inv_dist_sq; // power is -2 in 4D
        new_acc += diff * force_mag;
    }
    
    let dt = params.dt;
    particles[index].acc = new_acc;
    particles[index].vel += 0.5 * new_acc * dt;
}
