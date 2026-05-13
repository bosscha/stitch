push!(LOAD_PATH, joinpath(@__DIR__, "..", "src"))
using GaiaClustering
using NearestNeighbors
using Statistics
using LinearAlgebra
using ArgParse
using CSV
using DataFrames

# Note: volume_unit_sphere, digamma_int, kozachenko_leonenko_entropy,
# and compute_virial_ratio are now available in the GaiaClustering package.

function test_entropy()
    # Simulate a stellar cluster in 6D (3D Space + 3D Velocity)
    N_core = 500
    N_noise = 500
    D = 6
    
    println("Generating 6D phase-space data...\n")
    
    # 1. Structured Core (Low Entropy)
    # Simulates a tightly bound, relaxed cluster core
    core_data = randn(D, N_core) .* 0.5
    
    # 2. Uniform Background Noise (High Entropy)
    # Simulates random field stars across a large volume
    noise_data = (rand(D, N_noise) .- 0.5) .* 20.0 
    
    # 3. Mixed Cluster (Core + Halo/Noise)
    mixed_data = hcat(core_data, noise_data)
    
    println("--- k-NN Phase-Space Entropy (6D) ---")
    
    H_core = kozachenko_leonenko_entropy(core_data, k=3)
    println("Core only (tightly packed) : ", round(H_core, digits=3))
    
    H_noise = kozachenko_leonenko_entropy(noise_data, k=3)
    println("Noise only (uniform field) : ", round(H_noise, digits=3))
    
    H_mixed = kozachenko_leonenko_entropy(mixed_data, k=3)
    println("Mixed Cluster              : ", round(H_mixed, digits=3))
    
    println("\nInterpretation:")
    println("- Lower entropy = More clustered, dynamically relaxed, highly structured.")
    println("- Higher entropy = Dispersed, random, mixed.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    function parse_commandline()
        s = ArgParseSettings()
        @add_arg_table! s begin
            "--file", "-f"
                help = "Path to the .oc.csv file"
                arg_type = String
                default = ""
            "--k"
                help = "Number of neighbors for k-NN"
                arg_type = Int
                default = 3
            "--mass", "-m"
                help = "Mean mass per star in solar masses"
                arg_type = Float64
                default = 1.0
        end
        return parse_args(s)
    end

    function process_real_data(file_path::String, k::Int, mass::Float64)
        println("Reading real cluster data from: $file_path")
        df_raw = CSV.read(file_path, DataFrame)
        
        # Define columns we need for analysis (ignoring vrad as requested)
        analysis_cols = ["X", "Y", "Z", "vl", "vb"]
        
        # 1. Drop missing values in the columns we care about
        df = dropmissing(df_raw, analysis_cols)
        
        # 2. Ensure all values are finite (no NaN or Inf)
        for col in analysis_cols
            df = df[isfinite.(df[!, col]), :]
        end
        
        required_cols = ["X", "Y", "Z", "vl", "vb"]
        missing_cols = filter(c -> !(c in names(df)), required_cols)
        if !isempty(missing_cols)
            error("Missing required columns for 5D entropy/virial: ", join(missing_cols, ", "))
        end
        
        N = nrow(df)
        println("Found $N total stars (after cleaning NaN/Missing).")
        println("Assumed mean mass per star: $mass M_sun")
        
        # 1. Entropy and Virial for ALL stars
        if N <= k
            println("Not enough stars to compute $k-NN entropy for all stars.")
        else
            data_5d = Matrix{Float64}(df[:, required_cols])'
            H_5d = kozachenko_leonenko_entropy(data_5d, k=k)
            
            # Virial Ratio (Using 2D velocities + isotropy assumption)
            pos_all = Matrix{Float64}(df[:, ["X", "Y", "Z"]])'
            vel_all = Matrix{Float64}(df[:, ["vl", "vb"]])'
            Q_all = compute_virial_ratio(pos_all, vel_all, m=mass)

            println("\n--- ALL STARS ---")
            println("Entropy (H_5D) : ", round(H_5d, digits=3))
            println("Virial Ratio Q : ", round(Q_all, digits=3))
            println("(Using 2D velocities + isotropy assumption)")
        end

        # 2. Entropy and Virial for CORE stars (Type 1)
        if "type" in names(df)
            df_core = df[df.type .== 1, :]
            N_core = nrow(df_core)
            println("\nFound $N_core core stars (type 1).")
            
            if N_core <= k
                println("Not enough core stars to compute metrics.")
            else
                # 5D Phase-Space
                data_core_5d = Matrix{Float64}(df_core[:, required_cols])'
                H_core_5d = kozachenko_leonenko_entropy(data_core_5d, k=k)
                
                # Virial Ratio
                pos_core = Matrix{Float64}(df_core[:, ["X", "Y", "Z"]])'
                vel_core = Matrix{Float64}(df_core[:, ["vl", "vb"]])'
                Q_core = compute_virial_ratio(pos_core, vel_core, m=mass)

                println("--- CORE STARS ONLY ---")
                println("Entropy (H_5D) : ", round(H_core_5d, digits=3))
                println("Virial Ratio Q : ", round(Q_core, digits=3))
                
                # 3D Spatial
                spatial_core = Matrix{Float64}(df_core[:, ["X", "Y", "Z"]])'
                H_core_3d = kozachenko_leonenko_entropy(spatial_core, k=k)
                println("Spatial (H_3D) : ", round(H_core_3d, digits=3))
            end

            # 3. Entropy and Virial for TAIL stars (Type 2)
            df_tail = df[df.type .== 2, :]
            N_tail = nrow(df_tail)
            if N_tail > k
                println("\nFound $N_tail tail stars (type 2).")
                # 5D Phase-Space
                data_tail_5d = Matrix{Float64}(df_tail[:, required_cols])'
                H_tail_5d = kozachenko_leonenko_entropy(data_tail_5d, k=k)
                
                # Virial Ratio
                pos_tail = Matrix{Float64}(df_tail[:, ["X", "Y", "Z"]])'
                vel_tail = Matrix{Float64}(df_tail[:, ["vl", "vb"]])'
                Q_tail = compute_virial_ratio(pos_tail, vel_tail, m=mass)

                println("--- TAIL STARS ONLY ---")
                println("Entropy (H_5D) : ", round(H_tail_5d, digits=3))
                println("Virial Ratio Q : ", round(Q_tail, digits=3))
                
                # 3D Spatial
                spatial_tail = Matrix{Float64}(df_tail[:, ["X", "Y", "Z"]])'
                H_tail_3d = kozachenko_leonenko_entropy(spatial_tail, k=k)
                println("Spatial (H_3D) : ", round(H_tail_3d, digits=3))
            end
        else
            println("\nColumn 'type' not found. Skipping core/tail analysis.")
        end
        
        # 3. Spatial entropy for comparison (All stars)
        if nrow(df) > k
            spatial_data = Matrix{Float64}(df[:, ["X", "Y", "Z"]])'
            H_3d_spatial = kozachenko_leonenko_entropy(spatial_data, k=k)
            println("\n--- ALL STARS Spatial Entropy (3D) ---")
            println("Entropy (H_3D) : ", round(H_3d_spatial, digits=3))
        end
    end

    parsed_args = parse_commandline()
    file_path = parsed_args["file"]
    k_val = parsed_args["k"]
    mass_val = parsed_args["mass"]
    
    if file_path != ""
        process_real_data(file_path, k_val, mass_val)
    else
        test_entropy()
        println("\n(To test on real data, use: julia --project=. scripts/prototype_entropy.jl -f <path_to.oc.csv>)")
    end
end
