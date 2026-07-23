use rand::Rng;
use rand_distr::StandardNormal;
use std::time::Instant;
use postgres::{Client, NoTls};
use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;
use std::borrow::Cow;

const N_STARS: usize = 5000;
const DIM: usize = 4;
const G: f32 = 1.0;
const DT: f32 = 0.0001;
const STEPS: usize = 1000000;
const SOFTENING: f32 = 0.001;

const M_MIN: f32 = 0.1;
const M_MAX: f32 = 50.0;
const ALPHA: f32 = 2.35;

type VectorD = [f32; 4];

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct Particle {
    pos: [f32; 4],
    vel: [f32; 4],
    acc: [f32; 4],
    mass_pad: [f32; 4],
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct SimParams {
    num_particles: u32,
    g_const: f32,
    dt: f32,
    softening_sq: f32,
}

fn sample_king_nd(n_stars: usize, r_c: f32, r_t: f32) -> Vec<VectorD> {
    let mut rng = rand::thread_rng();
    let mut pos = Vec::with_capacity(n_stars);
    
    let p_max = {
        let mut max_p = 0.0f32;
        let num_points = 1000;
        for i in 0..num_points {
            let r = r_t * (i as f32) / (num_points as f32);
            let p = (r / r_t).powi((DIM - 1) as i32) * (1.0 / (1.0 + (r/r_c).powi(2)).sqrt() - 1.0 / (1.0 + (r_t/r_c).powi(2)).sqrt()).powi(2);
            if p > max_p { max_p = p; }
        }
        max_p * 1.1
    };

    while pos.len() < n_stars {
        let r_cand: f32 = rng.gen_range(0.0..r_t);
        let p_cand: f32 = rng.gen_range(0.0..p_max);
        let p_eval = (r_cand / r_t).powi((DIM - 1) as i32) * (1.0 / (1.0 + (r_cand/r_c).powi(2)).sqrt() - 1.0 / (1.0 + (r_t/r_c).powi(2)).sqrt()).powi(2);
        
        if p_cand < p_eval {
            let v: [f32; 4] = [
                rng.sample(StandardNormal),
                rng.sample(StandardNormal),
                rng.sample(StandardNormal),
                rng.sample(StandardNormal)
            ];
            let v_norm = (v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]).sqrt();
            pos.push([
                v[0] / v_norm * r_cand,
                v[1] / v_norm * r_cand,
                v[2] / v_norm * r_cand,
                v[3] / v_norm * r_cand,
            ]);
        }
    }
    pos
}

fn compute_mass_slope(masses: &[f32], m_min: f32, m_max: f32, num_bins: usize) -> f32 {
    if masses.len() < 10 { return f32::NAN; }
    
    let log_min = m_min.log10();
    let log_max = m_max.log10();
    let bin_width = (log_max - log_min) / (num_bins as f32 - 1.0);
    
    let mut bins = vec![0.0f32; num_bins];
    for i in 0..num_bins {
        bins[i] = 10.0f32.powf(log_min + (i as f32) * bin_width);
    }
    
    let mut counts = vec![0; num_bins - 1];
    for &m in masses {
        for i in 0..(num_bins - 1) {
            if m >= bins[i] && m < bins[i+1] {
                counts[i] += 1;
                break;
            }
        }
    }
    
    let mut dn_dm = vec![0.0f32; num_bins - 1];
    let mut bin_centers = vec![0.0f32; num_bins - 1];
    
    let mut valid_x = Vec::new();
    let mut valid_y = Vec::new();
    
    for i in 0..(num_bins - 1) {
        dn_dm[i] = counts[i] as f32 / (bins[i+1] - bins[i]);
        bin_centers[i] = (bins[i+1] * bins[i]).sqrt();
        
        if dn_dm[i] > 0.0 {
            valid_x.push(bin_centers[i].log10());
            valid_y.push(dn_dm[i].log10());
        }
    }
    
    if valid_x.len() < 3 { return f32::NAN; }
    
    let n = valid_x.len() as f32;
    let sum_x: f32 = valid_x.iter().sum();
    let sum_y: f32 = valid_y.iter().sum();
    let sum_xx: f32 = valid_x.iter().map(|x| x * x).sum();
    let sum_xy: f32 = valid_x.iter().zip(valid_y.iter()).map(|(x, y)| x * y).sum();
    
    let slope = (n * sum_xy - sum_x * sum_y) / (n * sum_xx - sum_x * sum_x);
    -slope
}

async fn run() -> Result<(), Box<dyn std::error::Error>> {
    println!("Starting Salpeter N-body simulation with wgpu GPU Acceleration (4D)...");
    
    let pos = sample_king_nd(N_STARS, 1.0, 10.0);
    
    let mut rng = rand::thread_rng();
    let mut vel: Vec<VectorD> = (0..N_STARS).map(|_| {
        let v: [f32; 4] = [
            rng.sample(StandardNormal),
            rng.sample(StandardNormal),
            rng.sample(StandardNormal),
            rng.sample(StandardNormal)
        ];
        [v[0] * 0.5, v[1] * 0.5, v[2] * 0.5, v[3] * 0.5]
    }).collect();

    // Generate Salpeter IMF
    let pow_idx = 1.0 - ALPHA;
    let m_diff = M_MAX.powf(pow_idx) - M_MIN.powf(pow_idx);
    
    let mut mass = vec![0.0f32; N_STARS];
    let mut total_mass = 0.0f32;
    
    for i in 0..N_STARS {
        let u: f32 = rng.gen();
        mass[i] = (u * m_diff + M_MIN.powf(pow_idx)).powf(1.0 / pow_idx);
        total_mass += mass[i];
    }
    
    for i in 0..N_STARS {
        mass[i] /= total_mass;
    }
    
    // Shift velocities to stop CoM drift
    let mut vel_cm = [0.0f32; 4];
    for i in 0..N_STARS {
        vel_cm[0] += vel[i][0] * mass[i];
        vel_cm[1] += vel[i][1] * mass[i];
        vel_cm[2] += vel[i][2] * mass[i];
        vel_cm[3] += vel[i][3] * mass[i];
    }
    
    let mut t_energy = 0.0f32;
    for i in 0..N_STARS {
        vel[i][0] -= vel_cm[0];
        vel[i][1] -= vel_cm[1];
        vel[i][2] -= vel_cm[2];
        vel[i][3] -= vel_cm[3];
        t_energy += 0.5 * mass[i] * (vel[i][0]*vel[i][0] + vel[i][1]*vel[i][1] + vel[i][2]*vel[i][2] + vel[i][3]*vel[i][3]);
    }

    // Compute potential energy
    let mut v_energy = 0.0f32;
    let mut initial_acc = vec![[0.0f32; 4]; N_STARS];
    for i in 0..N_STARS {
        for j in 0..N_STARS {
            if i == j { continue; }
            let dx = pos[j][0] - pos[i][0];
            let dy = pos[j][1] - pos[i][1];
            let dz = pos[j][2] - pos[i][2];
            let dw = pos[j][3] - pos[i][3];
            let dist_sq = dx*dx + dy*dy + dz*dz + dw*dw;
            let inv_dist_sq = 1.0 / (dist_sq + SOFTENING * SOFTENING);
            let inv_dist = inv_dist_sq.sqrt();
            
            if j > i {
                v_energy += -G * mass[i] * mass[j] * inv_dist;
            }
            
            // In 4D, force scales with 1/r^3. So force_mag = G * m * inv_dist_sq^2
            let force_mag = G * mass[j] * inv_dist_sq * inv_dist_sq;
            initial_acc[i][0] += force_mag * dx;
            initial_acc[i][1] += force_mag * dy;
            initial_acc[i][2] += force_mag * dz;
            initial_acc[i][3] += force_mag * dw;
        }
    }
    
    let f_scale = (0.5 * v_energy.abs() / t_energy).sqrt();
    for i in 0..N_STARS {
        vel[i][0] *= f_scale;
        vel[i][1] *= f_scale;
        vel[i][2] *= f_scale;
        vel[i][3] *= f_scale;
    }
    
    println!("➔ Initial Virial Ratio before scaling: {:.4}", t_energy / v_energy.abs());
    println!("➔ Scaled velocities by factor {:.4} to achieve virial equilibrium (Q = 0.5)", f_scale);

    let g_phys = 0.00449;
    let l_unit = 1.0;
    let time_to_myr = ((l_unit as f32).powi(3) / (g_phys * total_mass)).sqrt();

    let mut max_mass = 0.0f32;
    let mut min_mass = f32::MAX;
    for &m in &mass {
        if m > max_mass { max_mass = m; }
        if m < min_mass { min_mass = m; }
    }
    
    println!("IMF Generated.");
    println!("➔ Total initial mass: {:.2} M_sun", total_mass);
    println!("➔ Heaviest Star relative mass: {:.5}", max_mass);
    println!("➔ Lightest Star relative mass: {:.5}", min_mass);
    println!("➔ Time Conversion: 1 N-body time unit = {:.4} Myr", time_to_myr);

    let mut db_client = Client::connect("host=localhost user=stephane password=tallis dbname=hypercluster", NoTls).ok();
    
    if let Some(ref mut client) = db_client {
        let _ = client.execute("
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
        ", &[]);
    }

    // Calculate initial 80% radius
    let mut center_initial = [0.0f32; 4];
    for p in &pos {
        center_initial[0] += p[0] / N_STARS as f32;
        center_initial[1] += p[1] / N_STARS as f32;
        center_initial[2] += p[2] / N_STARS as f32;
        center_initial[3] += p[3] / N_STARS as f32;
    }
    
    let mut dist_initial: Vec<f32> = pos.iter().map(|p| {
        let dx = p[0] - center_initial[0];
        let dy = p[1] - center_initial[1];
        let dz = p[2] - center_initial[2];
        let dw = p[3] - center_initial[3];
        (dx*dx + dy*dy + dz*dz + dw*dw).sqrt()
    }).collect();
    
    dist_initial.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
    let r80_initial = dist_initial[(0.8 * N_STARS as f32) as usize];
    println!("➔ Initial 80% radius (R80): {:.4}", r80_initial);

    // Initialize wgpu
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends: wgpu::Backends::VULKAN,
        ..Default::default()
    });
    
    let adapter = instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        force_fallback_adapter: false,
        compatible_surface: None,
    }).await.ok_or("Failed to find suitable GPU adapter")?;
    
    println!("➔ Selected GPU Adapter: {:?}", adapter.get_info().name);
    
    let (device, queue) = adapter.request_device(
        &wgpu::DeviceDescriptor {
            label: Some("N-Body GPU"),
            required_features: wgpu::Features::empty(),
            required_limits: wgpu::Limits::default(),
            memory_hints: Default::default(),
        },
        None,
    ).await?;

    let mut particles = vec![Particle { pos: [0.0; 4], vel: [0.0; 4], acc: [0.0; 4], mass_pad: [0.0; 4] }; N_STARS];
    for i in 0..N_STARS {
        particles[i].pos = [pos[i][0], pos[i][1], pos[i][2], pos[i][3]];
        particles[i].vel = [vel[i][0], vel[i][1], vel[i][2], vel[i][3]];
        particles[i].acc = [initial_acc[i][0], initial_acc[i][1], initial_acc[i][2], initial_acc[i][3]];
        particles[i].mass_pad = [mass[i], 0.0, 0.0, 0.0];
    }
    
    let particle_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("Particle Buffer"),
        contents: bytemuck::cast_slice(&particles),
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC | wgpu::BufferUsages::COPY_DST,
    });
    
    let staging_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("Staging Buffer"),
        size: (particles.len() * std::mem::size_of::<Particle>()) as u64,
        usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    
    let sim_params = SimParams {
        num_particles: N_STARS as u32,
        g_const: G,
        dt: DT,
        softening_sq: SOFTENING * SOFTENING,
    };
    
    let param_buffer = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
        label: Some("Params Buffer"),
        contents: bytemuck::bytes_of(&sim_params),
        usage: wgpu::BufferUsages::UNIFORM,
    });

    let shader_src = include_str!("compute.wgsl");
    let cs_module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("N-Body Shader"),
        source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(shader_src)),
    });
    
    let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        entries: &[
            wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Storage { read_only: false },
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            },
            wgpu::BindGroupLayoutEntry {
                binding: 1,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }
        ],
        label: None,
    });
    
    let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
        layout: &bind_group_layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: particle_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: param_buffer.as_entire_binding(),
            }
        ],
        label: None,
    });
    
    let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: None,
        bind_group_layouts: &[&bind_group_layout],
        push_constant_ranges: &[],
    });
    
    let update_pos_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("Update Pos Pipeline"),
        layout: Some(&pipeline_layout),
        module: &cs_module,
        entry_point: "update_pos",
        compilation_options: Default::default(),
        cache: None,
    });
    
    let update_vel_pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
        label: Some("Update Vel Pipeline"),
        layout: Some(&pipeline_layout),
        module: &cs_module,
        entry_point: "update_vel",
        compilation_options: Default::default(),
        cache: None,
    });

    let track_interval = 500;
    let mut escaped_counts = Vec::new();
    let mut time_steps_list = Vec::new();
    let mass_phys: Vec<f32> = mass.iter().map(|&m| m * total_mass).collect();
    
    let start_time = Instant::now();
    let num_batches = STEPS / track_interval;
    let workgroup_count = (N_STARS as u32 + 255) / 256;

    for batch in 0..num_batches {
        for _ in 0..track_interval {
            let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            {
                let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor { label: None, timestamp_writes: None });
                cpass.set_pipeline(&update_pos_pipeline);
                cpass.set_bind_group(0, &bind_group, &[]);
                cpass.dispatch_workgroups(workgroup_count, 1, 1);
            }
            {
                let mut cpass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor { label: None, timestamp_writes: None });
                cpass.set_pipeline(&update_vel_pipeline);
                cpass.set_bind_group(0, &bind_group, &[]);
                cpass.dispatch_workgroups(workgroup_count, 1, 1);
            }
            queue.submit(Some(encoder.finish()));
        }
        
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        
        encoder.copy_buffer_to_buffer(&particle_buffer, 0, &staging_buffer, 0, (particles.len() * std::mem::size_of::<Particle>()) as u64);
        queue.submit(Some(encoder.finish()));
        
        let buffer_slice = staging_buffer.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        buffer_slice.map_async(wgpu::MapMode::Read, move |v| tx.send(v).unwrap());
        device.poll(wgpu::Maintain::Wait);
        
        if let Ok(Ok(())) = rx.recv() {
            let data = buffer_slice.get_mapped_range();
            let result: &[Particle] = bytemuck::cast_slice(&data);
            
            let step = (batch + 1) * track_interval;
            let mut center = [0.0f32; 4];
            for p in result {
                center[0] += p.pos[0] / N_STARS as f32;
                center[1] += p.pos[1] / N_STARS as f32;
                center[2] += p.pos[2] / N_STARS as f32;
                center[3] += p.pos[3] / N_STARS as f32;
            }
            
            let mut escaped = 0;
            for p in result {
                let dx = p.pos[0] - center[0];
                let dy = p.pos[1] - center[1];
                let dz = p.pos[2] - center[2];
                let dw = p.pos[3] - center[3];
                if (dx*dx + dy*dy + dz*dz + dw*dw).sqrt() > r80_initial {
                    escaped += 1;
                }
            }
            escaped_counts.push(escaped);
            
            let time_myr = (step as f32) * DT * time_to_myr;
            time_steps_list.push(time_myr);
            
            if let Some(ref mut client) = db_client {
                let mut trans = client.transaction()?;
                for i in 0..N_STARS {
                    let pos_vec = vec![result[i].pos[0] as f64, result[i].pos[1] as f64, result[i].pos[2] as f64, result[i].pos[3] as f64];
                    let vel_vec = vec![result[i].vel[0] as f64, result[i].vel[1] as f64, result[i].vel[2] as f64, result[i].vel[3] as f64];
                    let _ = trans.execute("
                        INSERT INTO star_snapshots 
                        (snapshot_id, time_myr, star_id, mass, dim_space, position, velocity) 
                        VALUES ($1, $2, $3, $4, $5, $6, $7)
                    ", &[&(step as i32), &(time_myr as f64), &(i as i32), &(mass_phys[i] as f64), &(DIM as i32), &pos_vec, &vel_vec]);
                }
                trans.commit()?;
            }
            
            if step % 1000 == 0 {
                println!("Step {:04}/{} completed. Physical Time: {:.2} Myr", step, STEPS, time_myr);
            }
            
            drop(data);
            staging_buffer.unmap();
        } else {
            println!("Failed to read buffer from GPU!");
            break;
        }
    }
    
    let duration = start_time.elapsed();
    println!("\nSimulation completed in: {:.4} seconds!", duration.as_secs_f64());
    
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    pollster::block_on(run())
}
