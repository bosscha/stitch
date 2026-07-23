#!/usr/bin/env julia

import Pkg
try
    using LibPQ, DataFrames, ArgParse, Plots
catch
    println("Installing missing dependencies...")
    Pkg.add(["LibPQ", "DataFrames", "ArgParse", "Plots"])
    using LibPQ, DataFrames, ArgParse, Plots
end

function get_bound_subsystem(pos::Matrix{Float64}, vel::Matrix{Float64}, mass::Vector{Float64}, G::Float64, softening::Float64; max_iter::Int=10)
    idx_bound = collect(1:length(mass))
    
    T_tot = 0.0
    V_tot = 0.0
    
    for iteration in 1:max_iter
        N_curr = length(idx_bound)
        if N_curr < 3
            return idx_bound, 0.0, 0.0
        end
        
        pos_curr = pos[idx_bound, :]
        vel_curr = vel[idx_bound, :]
        mass_curr = mass[idx_bound]
        
        total_m = sum(mass_curr)
        cm_pos = sum(pos_curr .* mass_curr, dims=1) ./ total_m
        cm_vel = sum(vel_curr .* mass_curr, dims=1) ./ total_m
        
        rel_pos = pos_curr .- cm_pos
        rel_vel = vel_curr .- cm_vel
        
        # Kinetic energy of each star
        T_i = 0.5 .* mass_curr .* vec(sum(rel_vel.^2, dims=2))
        
        # Potential energy of each star
        V_i = zeros(N_curr)
        
        Threads.@threads for i in 1:N_curr
            phi_i = 0.0
            for j in 1:N_curr
                if i != j
                    dist_sq = 0.0
                    for d in 1:size(pos_curr, 2)
                        dist_sq += (rel_pos[i, d] - rel_pos[j, d])^2
                    end
                    phi_i += mass_curr[j] / sqrt(dist_sq + softening^2)
                end
            end
            V_i[i] = -G * mass_curr[i] * phi_i
        end
        
        E_i = T_i .+ V_i
        
        bound_mask = E_i .< 0
        N_next = sum(bound_mask)
        
        T_tot = sum(T_i)
        V_tot = 0.5 * sum(V_i)
        
        if N_next == N_curr
            break
        end
        
        idx_bound = idx_bound[bound_mask]
    end
    
    return idx_bound, T_tot, V_tot
end

function parse_commandline()
    s = ArgParseSettings(description="Plot the evolution of the virial ratio of bound stars from the hypercluster database (Julia version).")

    @add_arg_table! s begin
        "--dim"
            help = "Spatial dimension of the simulation (e.g., 2, 3, 4, 5, 6, 100)"
            arg_type = Int
        "--host"
            help = "Database host"
            default = "localhost"
        "--user"
            help = "Database user"
            default = "stephane"
        "--password"
            help = "Database password"
            default = "tallis"
        "--dbname"
            help = "Database name"
            default = "hypercluster"
        "--softening"
            help = "Softening parameter used in simulation (default: 0.001)"
            arg_type = Float64
            default = 0.001
        "--g-constant"
            help = "Gravitational constant G (default: 1.0)"
            arg_type = Float64
            default = 1.0
        "--max-points"
            help = "Maximum number of snapshots to analyze (default: 200)"
            arg_type = Int
            default = 200
    end

    return parse_args(s)
end

function main()
    args = parse_commandline()
    
    dim_space = args["dim"]
    if dim_space === nothing
        print("Enter the spatial dimension (dim_space) to plot (e.g., 2, 3, 4, 5, 6, 100): ")
        dim_space_str = readline()
        dim_space = parse(Int, dim_space_str)
    end
    
    println("Connecting to database '$(args["dbname"])' on '$(args["host"])'...")
    
    conn_str = "dbname=$(args["dbname"]) user=$(args["user"]) password=$(args["password"]) host=$(args["host"])"
    conn = try
        LibPQ.Connection(conn_str)
    catch e
        println("Failed to connect to database: ", e)
        exit(1)
    end
    
    println("Connected successfully. Querying run history for dim_space = $dim_space...")
    
    # 1. max snapshot
    result = execute(conn, "SELECT MAX(snapshot_id) FROM star_snapshots WHERE dim_space = $dim_space;")
    max_snapshot_id = first(result).max
    if ismissing(max_snapshot_id)
        println("No simulation data found in the database for dim_space = $dim_space.")
        close(conn)
        exit(1)
    end
    
    # 2. min snapshot
    result = execute(conn, "SELECT MIN(snapshot_id) FROM star_snapshots WHERE dim_space = $dim_space;")
    min_snapshot_id = first(result).min
    
    # 3. get number of stars from the final snapshot to verify
    result = execute(conn, "SELECT COUNT(*) FROM star_snapshots WHERE dim_space = $dim_space AND snapshot_id = $max_snapshot_id;")
    num_stars = first(result).count
    
    # 4. start db id
    query_start_id = """
        SELECT MIN(id) FROM (
            SELECT id FROM star_snapshots 
            WHERE dim_space = $dim_space AND snapshot_id = $min_snapshot_id 
            ORDER BY id DESC 
            LIMIT $num_stars
        ) as sub;
    """
    result = execute(conn, query_start_id)
    start_db_id = first(result).min
    
    # 5. step size
    query_step = """
        SELECT MIN(snapshot_id) 
        FROM star_snapshots 
        WHERE dim_space = $dim_space AND snapshot_id > $min_snapshot_id AND id >= $start_db_id;
    """
    result = execute(conn, query_step)
    second_snap = first(result).min
    step_size = !ismissing(second_snap) ? (second_snap - min_snapshot_id) : 100
    
    snap_ids = collect(min_snapshot_id:step_size:max_snapshot_id)
    total_snapshots = length(snap_ids)
    
    println("Detected snapshots up to ID $max_snapshot_id (step size: $step_size). Total potential snapshots: $total_snapshots")
    
    sample_interval = max(1, total_snapshots ÷ args["max-points"])
    sampled_snap_ids = snap_ids[1:sample_interval:end]
    println("Downsampled to $(length(sampled_snap_ids)) snapshots for analysis (interval: every $sample_interval steps).")
    
    # get n_stars for actual query limit
    query_nstars = """
        SELECT COUNT(*) 
        FROM star_snapshots 
        WHERE dim_space = $dim_space AND snapshot_id = $(sampled_snap_ids[1]) AND id >= $start_db_id;
    """
    result = execute(conn, query_nstars)
    n_stars = first(result).count
    println("Number of stars per snapshot: $n_stars")
    
    times = Float64[]
    virial_ratios = Float64[]
    virial_ratios_all = Float64[]
    bound_fractions = Float64[]
    
    for (idx, snap_id) in enumerate(sampled_snap_ids)
        print("\rAnalyzing snapshot $idx/$(length(sampled_snap_ids)) (ID: $snap_id)...")
        
        query_snap = """
            SELECT position, velocity, mass, time_myr FROM (
                SELECT id, position, velocity, mass, time_myr 
                FROM star_snapshots 
                WHERE dim_space = $dim_space AND snapshot_id = $snap_id AND id >= $start_db_id
                ORDER BY id DESC 
                LIMIT $n_stars
            ) as sub 
            ORDER BY id ASC;
        """
        result = execute(conn, query_snap)
        df = DataFrame(result)
        
        if nrow(df) == 0
            continue
        end
        
        pos = Array{Float64, 2}(undef, nrow(df), dim_space)
        vel = Array{Float64, 2}(undef, nrow(df), dim_space)
        mass = zeros(Float64, nrow(df))
        
        for i in 1:nrow(df)
            for d in 1:dim_space
                # Handle cases where LibPQ arrays might be nested or direct
                pos[i, d] = df.position[i][d]
                vel[i, d] = df.velocity[i][d]
            end
            mass[i] = df.mass[i]
        end
        
        mass ./= sum(mass)
        time_val = Float64(df.time_myr[1])
        
        bound_indices, T_bound, V_bound = get_bound_subsystem(pos, vel, mass, args["g-constant"], args["softening"])
        
        Q_bound = V_bound != 0 ? T_bound / abs(V_bound) : NaN
        
        # Entire cluster
        cm_pos_all = sum(pos .* mass, dims=1) ./ sum(mass)
        cm_vel_all = sum(vel .* mass, dims=1) ./ sum(mass)
        
        rel_pos_all = pos .- cm_pos_all
        rel_vel_all = vel .- cm_vel_all
        
        T_all = 0.5 * sum(mass .* vec(sum(rel_vel_all.^2, dims=2)))
        
        V_i_all = zeros(nrow(df))
        Threads.@threads for i in 1:nrow(df)
            phi_i = 0.0
            for j in 1:nrow(df)
                if i != j
                    dist_sq = 0.0
                    for d in 1:dim_space
                        dist_sq += (rel_pos_all[i, d] - rel_pos_all[j, d])^2
                    end
                    phi_i += mass[j] / sqrt(dist_sq + args["softening"]^2)
                end
            end
            V_i_all[i] = -args["g-constant"] * mass[i] * phi_i
        end
        
        V_all = 0.5 * sum(V_i_all)
        Q_all = V_all != 0 ? T_all / abs(V_all) : NaN
        
        bound_frac = length(bound_indices) / nrow(df)
        
        push!(times, time_val)
        push!(virial_ratios, Q_bound)
        push!(virial_ratios_all, Q_all)
        push!(bound_fractions, bound_frac)
    end
    
    println("\nAnalysis complete. Generating plots...")
    close(conn)
    
    p1 = plot(times, virial_ratios, color=:crimson, label="Bound Core Virial Q = T/|V|", linewidth=2)
    plot!(p1, times, virial_ratios_all, color=:darkorange, linestyle=:dash, label="Entire Cluster Virial Q = T/|V|", linewidth=2)
    hline!(p1, [0.5], color=:gray, linestyle=:dot, label="Virial Equilibrium (Q=0.5)")
    ylabel!(p1, "Virial Ratio Q")
    title!(p1, "Evolution of the Virial Ratio (N=$dim_space dimensions)")
    
    p2 = plot(times, bound_fractions .* 100, color=:royalblue, label="Bound Star %", linewidth=2)
    xlabel!(p2, "Time (Myr)")
    ylabel!(p2, "Bound Stars (%)")
    
    p = plot(p1, p2, layout=(2, 1), size=(1000, 800), link=:x)
    
    plot_filename = "virial_evolution_$(dim_space)D_julia.png"
    savefig(p, plot_filename)
    println("Plot successfully saved to '$plot_filename'!")
end

main()
