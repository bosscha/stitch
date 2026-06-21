## Test script to download from local DB the Gaia pixel of NGC 2682 (M67),
## run cluster extraction using custom.ext via extra.jl,
## and upload results to local PostgreSQL gaiadb.

using LibPQ, DataFrames, PyCall

println("==========================================================")
println("## Starting NGC 2682 (M67) Gaia Clustering & DB test script")
println("==========================================================")

# Ensure target directory exists
target_dir = "./test_run_ngc2682"
if !isdir(target_dir)
    println("## Creating directory: $target_dir")
    mkpath(target_dir)
end

# 1. Connect to local PostgreSQL database
db_host = "192.168.50.228"
db_user = "stephane"
db_pass = "tallis"
db_name = "gaiadb"
conn_str = "host=$db_host user=$db_user dbname=$db_name password=$db_pass"

println("## Connecting to PostgreSQL database on $db_host...")
conn = LibPQ.Connection(conn_str)

# 2. Query gaia_source for pixel 1074 (NGC 2682's level 5 HEALPix pixel)
# Level 5 HEALPix is stored at bit 49 of the source_id
pixel_id = 1074
query = "SELECT * FROM gaia_source WHERE (source_id >> 49) = $pixel_id;"
println("## Downloading stars for pixel $pixel_id from gaia_source...")
res = execute(conn, query)
df = DataFrame(res)
close(conn)

num_stars = size(df, 1)
println("## Downloaded $num_stars stars.")
if num_stars == 0
    error("## ERROR: No stars retrieved for pixel $pixel_id! Ensure gaia_source is populated.")
end

# Add columns expected by filter_data that are absent in the local PostgreSQL schema
df[!, :pmra_error] .= 0.0
df[!, :pmdec_error] .= 0.0

# 3. Clean missing values to prepare for PyCall conversion
println("## Cleaning missing values in retrieved data...")
for col in names(df)
    col_type = eltype(df[!, col])
    if Union{Missing, Float64} <: col_type || Union{Missing, Float32} <: col_type
        df[!, col] = [ismissing(x) ? NaN : Float64(x) for x in df[!, col]]
    elseif Union{Missing, Int64} <: col_type || Union{Missing, Int32} <: col_type
        df[!, col] = [ismissing(x) ? 0 : Int64(x) for x in df[!, col]]
    elseif Union{Missing, String} <: col_type
        df[!, col] = [ismissing(x) ? "" : String(x) for x in df[!, col]]
    end
end

# 4. Save to VOTable format using Astropy via PyCall
votable_path = joinpath(target_dir, "NGC2682.vot")
println("## Converting to Astropy Table and writing VOTable to $votable_path...")
py"""
from astropy.table import Table
import pandas as pd

def write_votable(data_dict, filepath):
    df = pd.DataFrame(data_dict)
    # clean object/none types
    for col in df.columns:
        df[col] = df[col].apply(lambda x: None if (isinstance(x, float) and pd.isna(x)) or str(x) == "missing" else x)
    t = Table.from_pandas(df)
    t.write(filepath, format="votable", overwrite=True)
"""

data_dict = Dict(string(col) => df[!, col] for col in names(df))
py"write_votable"(data_dict, votable_path)
println("## VOTable saved successfully.")

# 5. Run the Julia implementation in scripts/extra.jl
extra_script = "./scripts/extra.jl"
config_file = "./test/custom.ext"
println("## Running cluster extraction: julia --project=. $extra_script -m $config_file...")
run(`julia --project=. $extra_script -m $config_file`)

# 6. Verify that the cluster was successfully uploaded to PostgreSQL
println("## Verifying database upload of the extracted clusters...")
conn = LibPQ.Connection(conn_str)
res_verify = execute(conn, "SELECT c.cluster_id, count(*) FROM clusters c JOIN clusters_metadata m ON c.cluster_id = m.cluster_id WHERE m.votname LIKE '%NGC2682%' GROUP BY c.cluster_id;")
df_verify = DataFrame(res_verify)
close(conn)

println("## Verification query result:")
show(df_verify)
println()

if size(df_verify, 1) > 0
    println("## SUCCESS: Extracted clusters and their stars were successfully uploaded to postgresql gaiadb!")
else
    error("## FAILURE: No cluster entries for NGC2682 found in the clusters table after run.")
end

println("==========================================================")
println("## Test script finished successfully!")
println("==========================================================")
