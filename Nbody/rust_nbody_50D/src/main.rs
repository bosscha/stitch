use rand::Rng;
use rand_distr::StandardNormal;
use std::time::Instant;
use postgres::{Client, NoTls};
use bytemuck::{Pod, Zeroable};
use wgpu::util::DeviceExt;
use std::borrow::Cow;

const N_STARS: usize = 10000;
const DIM: usize = 50;
const MAX_DIM: usize = 64;
const G: f32 = 1.0;
const DT: f32 = 0.0001;
const STEPS: usize = 5000000;
const SOFTENING: f32 = 0.01;

const M_MIN: f32 = 0.1;
const M_MAX: f32 = 50.0;
const ALPHA: f32 = 2.35;

type VectorD = [f32; MAX_DIM];

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
struct Particle {
    pos: [f32; MAX_DIM],
    vel: [f32; MAX_DIM],
    acc: [f32; MAX_DIM],
    mass_pad: [f32; MAX_DIM],
}

#[repr(C)]
#[derive(Copy, Clone, Debug, Pod, Zeroable)]
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

fn sample_king_nd(n_stars: usize, r_c: f32, r_t: f32) -> Vec<VectorD> {
    let mut rng = rand::thread_rng();
    let mut pos = Vec::with_capacity(n_stars);
    
    let p_max = {
        let mut max_p = 0.0f32;
        let num_points = 1000;
        for i in 0..num_points {
            let r = r_t * (i as f32) / (num_points as f32);
            let p = r.powi((DIM - 1) as i32) * (1.0 / (1.0 + (r/r_c).powi(2)).sqrt() - 1.0 / (1.0 + (r_t/r_c).powi(2)).sqrt()).powi(2);
            if p > max_p { max_p = p; }
        }
        max_p * 1.1
    };

    while pos.len() < n_stars {
        let r_cand: f32 = rng.gen_range(0.0..r_t);
        let p_cand: f32 = rng.gen_range(0.0..p_max);
        let p_eval = r_cand.powi((DIM - 1) as i32) * (1.0 / (1.0 + (r_cand/r_c).powi(2)).sqrt() - 1.0 / (1.0 + (r_t/r_c).powi(2)).sqrt()).powi(2);
        
        if p_cand < p_eval {
            let mut v = [0.0f32; MAX_DIM];
            let mut v_norm_sq = 0.0;
            for i in 0..DIM {
                v[i] = rng.sample(StandardNormal);
                v_norm_sq += v[i] * v[i];
            }
            let v_norm = v_norm_sq.sqrt();
            
            let mut p = [0.0f32; MAX_DIM];
            for i in 0..DIM {
                p[i] = v[i] / v_norm * r_cand;
            }
            pos.push(p);
        }
    }
    pos
}

async fn run() -> Result<(), Box<dyn std::error::Error>> {
    println!("Starting Salpeter N-body simulation with wgpu GPU Acceleration ({}D)...", DIM);
    
    let pos = sample_king_nd(N_STARS, 1.0, 10.0);
    
    let mut rng = rand::thread_rng();
    let mut vel: Vec<VectorD> = (0..N_STARS).map(|_| {
        let mut v = [0.0f32; MAX_DIM];
        for i in 0..DIM {
            let val: f32 = rng.sample(StandardNormal);
            v[i] = val * 0.5;
        }
        v
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
    let mut vel_cm = [0.0f32; MAX_DIM];
    for i in 0..N_STARS {
        for k in 0..DIM {
            vel_cm[k] += vel[i][k] * mass[i];
        }
    }
    
    let mut t_energy = 0.0f32;
    for i in 0..N_STARS {
        for k in 0..DIM {
            vel[i][k] -= vel_cm[k];
            t_energy += 0.5 * mass[i] * vel[i][k] * vel[i][k];
        }
    }

    // Compute potential energy
    let mut v_energy = 0.0f32;
    let mut initial_acc = vec![[0.0f32; MAX_DIM]; N_STARS];
    for i in 0..N_STARS {
        for j in 0..N_STARS {
            if i == j { continue; }
            let mut dist_sq = 0.0;
            let mut diffs = [0.0f32; MAX_DIM];
            for k in 0..DIM {
                diffs[k] = pos[j][k] - pos[i][k];
                dist_sq += diffs[k] * diffs[k];
            }
            
            let inv_dist_sq = 1.0 / (dist_sq + SOFTENING * SOFTENING);
            let inv_dist = inv_dist_sq.sqrt();
            
            if j > i {
                v_energy += -G * mass[i] * mass[j] * inv_dist;
            }
            
            let force_mag = G * mass[j] * inv_dist.powf(DIM as f32);
            for k in 0..DIM {
                initial_acc[i][k] += force_mag * diffs[k];
            }
        }
    }
    
    let f_scale = (0.5 * v_energy.abs() / t_energy).sqrt();
    for i in 0..N_STARS {
        for k in 0..DIM {
            vel[i][k] *= f_scale;
        }
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
    let mut center_initial = [0.0f32; MAX_DIM];
    for p in &pos {
        for k in 0..DIM {
            center_initial[k] += p[k] / N_STARS as f32;
        }
    }
    
    let mut dist_initial: Vec<f32> = pos.iter().map(|p| {
        let mut dist_sq = 0.0;
        for k in 0..DIM {
            let dx = p[k] - center_initial[k];
            dist_sq += dx * dx;
        }
        dist_sq.sqrt()
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

    let mut particles = vec![Particle { pos: [0.0; MAX_DIM], vel: [0.0; MAX_DIM], acc: [0.0; MAX_DIM], mass_pad: [0.0; MAX_DIM] }; N_STARS];
    for i in 0..N_STARS {
        for k in 0..DIM {
            particles[i].pos[k] = pos[i][k];
            particles[i].vel[k] = vel[i][k];
            particles[i].acc[k] = initial_acc[i][k];
        }
        particles[i].mass_pad[0] = mass[i];
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
        dim: DIM as u32,
        pad1: 0,
        pad2: 0,
        pad3: 0,
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
            let mut center = [0.0f32; MAX_DIM];
            for p in result {
                for k in 0..DIM {
                    center[k] += p.pos[k] / N_STARS as f32;
                }
            }
            
            let mut escaped = 0;
            for p in result {
                let mut dist_sq = 0.0;
                for k in 0..DIM {
                    let dx = p.pos[k] - center[k];
                    dist_sq += dx * dx;
                }
                if dist_sq.sqrt() > r80_initial {
                    escaped += 1;
                }
            }
            escaped_counts.push(escaped);
            
            let time_myr = (step as f32) * DT * time_to_myr;
            time_steps_list.push(time_myr);
            
            if let Some(ref mut client) = db_client {
                let mut trans = client.transaction()?;
                for i in 0..N_STARS {
                    let mut pos_vec = Vec::new();
                    let mut vel_vec = Vec::new();
                    for k in 0..DIM {
                        pos_vec.push(result[i].pos[k] as f64);
                        vel_vec.push(result[i].vel[k] as f64);
                    }
                    
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
