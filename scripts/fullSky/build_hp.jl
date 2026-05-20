## standalone script to (re)build sets of candidate cluster (cycle) from a postgresql DB using HEALPix index
## Can be resumed  if it fails.
## The open cluster candidates are saved in the postgresql DB only if they are in the central pixel. All 8 neighbours 
## are processed with the central one for edge cases, HEALPix level 5 is used for all-sky

using DataFrames, CSV, TOML, ArgParse, LibPQ
using Dates, Printf

using PyCall

# This block executes native Python code inside PyCall's environment
# and permanently disables SSL verification for this run.
py"""
import ssl
try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context
"""



rootdir = ENV["GAIA_ROOT"]

push!(LOAD_PATH, "$rootdir/run/src")
using GaiaClustering



########################  HP 
### All-sky processing
function fullSky(meta)
    tstart = now()
    println(blue("## Processing Full Sky with HEALPix level 5 ..."))
    println(blue("## Starting at $tstart"))

    if !haskey(meta, "fullsky")
        error("The meta TOML file is missing the [fullsky] section!")
    end

    mfull = meta["fullsky"]
    mgene = meta["general"]
    mextra = read_params(mfull["extrafile"], false)

    # ensure directories exist and override result paths
    wdir = abspath(mgene["wdir"])
    mkpath(wdir)
    
    # Force results and plots into dedicated subdirectories
    mextra.ocdir = joinpath(wdir, "results")
    mextra.plotdir = joinpath(wdir, "plots")
    
    mkpath(mextra.ocdir)
    mkpath(mextra.plotdir)
    
    cd(wdir)

    progress_table = "$(mextra.dbtable)_processed_pixels"
    
    # Connect to DB to load progress
    pwd_arg = (mextra.dbpass != "") ? "password=$(mextra.dbpass)" : ""
    conn_str = "host=$(mextra.dbhost) user=$(mextra.dbuser) dbname=$(mextra.dbname) $pwd_arg"
    
    conn = LibPQ.Connection(conn_str)
    execute(conn, "CREATE TABLE IF NOT EXISTS $progress_table (pix BIGINT PRIMARY KEY, datetime VARCHAR(255));")
    # Ensure column exists for existing tables
    execute(conn, "ALTER TABLE $progress_table ADD COLUMN IF NOT EXISTS datetime VARCHAR(255);")
    res = execute(conn, "SELECT pix FROM $progress_table;")
    dfp = DataFrame(res)
    close(conn)
    
    if !("pix" in names(dfp))
        dfp = DataFrame(pix=Int[])
    end

    if mextra.optim == "yes"
        optim = true
    elseif mextra.optim == "no"
        optim = false
    end

    # Initialize Healpix
    healpy = pyimport("healpy")
    nside = 32
    npix = 12288

    println("## Total HEALPix Level 5 pixels: $npix")

    nbatch = mextra.nbatch
    printstyled("## Running with $nbatch simultaneous batches\n", color=:blue)

    file_lock = ReentrantLock()
    active_lock_file = joinpath(wdir, "active_pixels.lock")

    function process_pixel_inner(P)
        if P in dfp.pix
            println("## Pixel $P already processed. Skipping...")
            return
        end

        lock(file_lock)
        active_pixels = isfile(active_lock_file) ? readlines(active_lock_file) : String[]
        my_pid = getpid()
        
        valid_locks = String[]
        for line in active_pixels
            parts = split(line, ",")
            if length(parts) == 2
                locked_P, locked_pid = parts
                if isdir("/proc/" * locked_pid)
                    push!(valid_locks, line)
                else
                    println("## Removing stale lock for pixel $locked_P (PID $locked_pid is dead)")
                end
            end
        end
        
        if any(startswith(l, string(P) * ",") for l in valid_locks)
            println("## Pixel $P is currently being processed (locked). Skipping...")
            open(active_lock_file, "w") do f
                for l in valid_locks
                    println(f, l)
                end
            end
            unlock(file_lock)
            return
        end
        
        push!(valid_locks, "$(P),$(my_pid)")
        open(active_lock_file, "w") do f
            for l in valid_locks
                println(f, l)
            end
        end
        unlock(file_lock)

        println("=====================================================")
        println("## Processing central pixel: $P")

        neighbors = healpy.get_all_neighbours(nside, P, nest=true)
        # Filter out -1 which indicates missing neighbor (e.g. at corners, though nside=32 has neighbors)
        valid_neighbors = filter(x -> x != -1, neighbors)
        target_pixels = vcat([P], valid_neighbors)

        mextra_local = deepcopy(mextra)
        mextra_local.votname = "HP5_$P"

        pwd_arg = (mextra_local.dbpass != "") ? "password=$(mextra_local.dbpass)" : ""
        conn_str = "host=$(mextra_local.dbhost) user=$(mextra_local.dbuser) dbname=$(mextra_local.dbname) $pwd_arg"

        try
            df, dfcart, dfcartnorm = get_data_pg(mextra_local, target_pixels)
            extra_db(mextra_local, optim, df, dfcart, dfcartnorm)
            
            # --- Cleanup out-of-bounds clusters ---
            # Any cluster saved to the DB must have its center within the central pixel P.
            # If not, we remove it from the DB.
            if mextra_local.savedb == "yes"
                conn_local = LibPQ.Connection(conn_str)
                
                # Fetch clusters from this votname
                query = "SELECT cluster_id, ra, dec FROM $(mextra_local.dbtable)_metadata WHERE votname = \$1;"
                res = execute(conn_local, query, [mextra_local.votname])
                df_clusters = DataFrame(res)
                
                for row in eachrow(df_clusters)
                    # compute cluster's healpix
                    cpix = healpy.ang2pix(nside, row.ra, row.dec, nest=true, lonlat=true)
                    if cpix != P
                        println("## Cleanup: Cluster $(row.cluster_id) is centered in pixel $cpix instead of $P. Deleting from DB...")
                        execute(conn_local, "DELETE FROM $(mextra_local.dbtable) WHERE cluster_id = \$1;", [row.cluster_id])
                        execute(conn_local, "DELETE FROM $(mextra_local.dbtable)_metadata WHERE cluster_id = \$1;", [row.cluster_id])
                    else
                        println("## Cluster $(row.cluster_id) confirmed in central pixel $P.")
                    end
                end
                close(conn_local)
            end
            
            # Save progress
            conn_local2 = LibPQ.Connection(conn_str)
            execute(conn_local2, "INSERT INTO $progress_table (pix, datetime) VALUES (\$1, \$2) ON CONFLICT (pix) DO UPDATE SET datetime = EXCLUDED.datetime;", [P, Dates.format(now(), "yyyy-mm-dd HH:MM:SS")])
            close(conn_local2)
            
            # Update progress plot after saving
            lock(file_lock)
            push!(dfp, [P])
            plot_hp_sky(dfp.pix, nside, figname="allsky_progress.png")
            unlock(file_lock)
            
        catch e
            println("## Error processing pixel $P : $e")
            if isa(e, ErrorException) && occursin("No data found", e.msg)
                println("## No stars fetched for $P, recording as done.")
                conn_local2 = LibPQ.Connection(conn_str)
                execute(conn_local2, "INSERT INTO $progress_table (pix, datetime) VALUES (\$1, \$2) ON CONFLICT (pix) DO UPDATE SET datetime = EXCLUDED.datetime;", [P, Dates.format(now(), "yyyy-mm-dd HH:MM:SS")])
                close(conn_local2)
                
                # Update progress plot after saving
                lock(file_lock)
                push!(dfp, [P])
                plot_hp_sky(dfp.pix, nside, figname="allsky_progress.png")
                unlock(file_lock)
            else
                println("## Unexpected error. Resuming on next run might retry this pixel.")
                # Depending on how the user wants to handle db drops, we may want to break or continue.
                # Continuing allows robustness against random failure of a specific pixel
                # break
            end
        finally
            lock(file_lock)
            if isfile(active_lock_file)
                active_pixels = readlines(active_lock_file)
                open(active_lock_file, "w") do f
                    for l in active_pixels
                        if !startswith(l, string(P) * ",") && l != string(P)
                            println(f, l)
                        end
                    end
                end
            end
            unlock(file_lock)
        end
    end

    function process_pixel(P)
        task_local_storage(:pix, P) do
            process_pixel_inner(P)
        end
    end

    if nbatch > 1
        asyncmap(process_pixel, 0:(npix-1); ntasks=nbatch)
    else
        for P in 0:(npix-1)
            process_pixel(P)
        end
    end
    
    println("## Full Sky Processing completed...")
end


#########################  reprocess function
function reprocess(meta)

    tstart = now()
    println(blue("## Reprocessing ..."))
    println(blue("## Starting at $tstart"))
    mrepro = meta["reprocess"]
    mextra = read_params(mrepro["extrafile"], false)
    println(blue("## Clustering algorithm: $(mextra.algo)"))
    # Force results and plots into dedicated subdirectories
    wdir = abspath(mgene["wdir"])
    mextra.ocdir = joinpath(wdir, "results")
    mextra.plotdir = joinpath(wdir, "plots")

    checkdir(mextra.ocdir, mextra.plotdir)

    debug_red(mextra.plotdir)

    cd(mgene["wdir"])

    progressfile = "_done.csv"            #progress file 

    # if !haskey(mextra, :rootdir) mextra.rootdir= "./" end
    # if !haskey(mextra, :wdir) mextra.wdir= "./" end
    #if !haskey(mextra, :plotdir) mextra.plotdir= "./plotSelect" end
    # if !haskey(mextra, :ocdir) mextra.ocdir= "./oc" end

    println(mrepro)
    println(mgene)

    if isfile(progressfile)
        dfp = CSV.File(progressfile, delim=",") |> DataFrame
    else
        dfp = DataFrame(votname=String[])
    end

    dfblck = get_blacklist(mgene)

    Kdeg = 57.69

    if mextra.optim == "no"
        optim = false
        println("## Using preprocessed file : $(mrepro["optsol"])")
        ## delim Guess
        fline = readline(mrepro["optsol"])
        if occursin(";", fline)
            delimiter = ";"
        else
            delimiter = ","
        end
        dfoptsol = CSV.File(mrepro["optsol"], delim=delimiter) |> DataFrame

        gaia = pyimport("astroquery.gaia")

        radone = []
        decdone = []
        for row in eachrow(dfoptsol)
            if mgene["getvot"] == "yes"
                rect = false
                ra = row.ra
                dec = row.dec
                push!(radone, ra)
                push!(decdone, dec)
                tol = mrepro["tol"]
                radius = mrepro["radius"]

                estangle = Kdeg * 25 / row.distance             # estimated field size in degree (if angle small for tangent)
                debug_red("est. angle $estangle deg")

                name = @sprintf("RA%.3fDec%.3f", ra, dec)
                votname = @sprintf("%s-%2.1fdeg.vot", name, radius)

                if votname in dfp.votname || votname in dfblck.votname
                    println("## $votname skipped..")
                else

                    mextra.votname = get_gaia_data_many(gaia, radius, tol, ra, dec, name, rect)

                    if haskey(row, :w3dm)
                        mextra.w3d = row.w3dm
                    else
                        mextra.w3d = row.w3d
                    end
                    if haskey(row, :wvelm)
                        mextra.wvel = row.wvelm
                    else
                        mextra.wvel = row.wvel
                    end
                    if haskey(row, :whrdm)
                        mextra.whrd = row.whrdm
                    else
                        mextra.whrd = row.whrd
                    end
                    if haskey(row, :epsm)
                        mextra.eps = row.epsm
                    else
                        mextra.eps = row.eps
                    end
                    if haskey(row, :mclm)
                        mextra.mcl = Int(floor(row.mclm))
                    else
                        mextra.mcl = Int(floor(row.mcl))
                    end
                    if haskey(row, :mneim)
                        mextra.mnei = Int(floor(row.mneim))
                    else
                        mextra.mnei = Int(floor(row.mnei))
                    end

                    #smextra.wvel= 8.0
                    debug_red("Weights : $(mextra.w3d) $(mextra.wvel) $(mextra.whrd)  ")


                    extra(mextra, optim)
                    push!(dfp, [mextra.votname])
                    CSV.write(progressfile, dfp, delim=";")

                    if mgene["rmvot"] == "yes"
                        rm(mextra.votname)
                        println("## votable $(mextra.votname) removed")
                    end
                end
                plot_sky(radone, decdone, radius=50, figname="reprocess-allsky.png")
            end
        end
    end
    println("## Reprocessing completed...")
    rm(progressfile)
end

#########################  randomfields function
function randomfields(meta)
    tstart = now()
    println(blue("## Processing random fields ..."))
    println(blue("## Starting at $tstart"))

    mrandom = meta["random"]
    mgene = meta["general"]
    mextra = read_params(mrandom["extrafile"], false)

    wdir = abspath(mgene["wdir"])
    cd(wdir)
    progressfile = "_done.csv"            #progress file 
    mextra.rootdir = wdir
    mextra.wdir = wdir
    mextra.plotdir = joinpath(wdir, "plots")
    mextra.ocdir = joinpath(wdir, "results")
    
    mkpath(mextra.plotdir)
    mkpath(mextra.ocdir)

    println(mrandom)
    println(mgene)

    if isfile(progressfile)
        dfp = CSV.File(progressfile, delim=",") |> DataFrame
    else
        dfp = DataFrame(votname=String[])
    end

    ndone = size(dfp)[1]
    nfields = mrandom["fields"]
    radius = mrandom["radius"]
    tol = mrandom["tol"]
    mode = mrandom["mode"]

    if mode == "galactic"
        bscale = mrandom["bscale"]
    else
        bscale = -1
    end
    rect = false   # conesearch

    if ndone < nfields
        notfinished = true
    else
        notfinished = false
    end

    gaia = pyimport("astroquery.gaia")

    if mextra.optim == "yes"
        optim = true
    elseif mextra.optim == "no"
        optim = false
    end

    while notfinished
        ra, dec = get_random_field(mode, bscale)
        name = @sprintf("RA%.3fDec%.3f", ra, dec)
        debug_red(name)

        mextra.votname = get_gaia_data_many(gaia, radius, tol, ra, dec, name, rect)
        extra(mextra, optim)

        push!(dfp, [mextra.votname])
        CSV.write(progressfile, dfp, delim=";")
        if mgene["rmvot"] == "yes"
            rm(mextra.votname)
            println("## votable $(mextra.votname) removed")
        end

        ndone += 1
        println("## $ndone random fields processed...")
        if ndone > nfields
            notfinished = false
        end
    end
    println("## Random fields completed...")
    rm(progressfile)
end

#########################  gridding function
function gridding(meta)

    tstart = now()
    println(blue("## Processing gridding fields ..."))
    println(blue("## Starting at $tstart"))

    mgrid = meta["gridding"]
    mgene = meta["general"]
    mextra = read_params(mgrid["extrafile"], false)

    wdir = abspath(mgene["wdir"])
    cd(wdir)
    progressfile = "_done.csv"            #progress file 
    mextra.rootdir = wdir
    mextra.wdir = wdir
    mextra.plotdir = joinpath(wdir, "plots")
    mextra.ocdir = joinpath(wdir, "results")
    
    # ensure directories exist so we don't get PyPlot write errors!
    mkpath(mextra.plotdir)
    mkpath(mextra.ocdir)

    println(mgrid)
    println(mgene)

    if isfile(progressfile)
        dfp = CSV.File(progressfile, delim=",") |> DataFrame
    else
        dfp = DataFrame(votname=String[])
    end

    gaia = pyimport("astroquery.gaia")

    radius = mgrid["radius"]
    tol = mgrid["tol"]
    ref = mgrid["ref"]
    rect = false   # conesearch

    if ref == "galactic"
        xiter = mgrid["lgal_grid"]
        yiter = mgrid["bgal_grid"]
    elseif ref == "equatorial"
        xiter = mgrid["ra_grid"]
        yiter = mgrid["dec_grid"]
    end

    if mextra.optim == "yes"
        optim = true
    elseif mextra.optim == "no"
        optim = false
    end

    debug_red(xiter)
    radone = []
    decdone = []

    for xx in xiter[1]:xiter[3]:xiter[2]
        for yy in yiter[1]:yiter[3]:yiter[2]
            if ref == "galactic"
                println("## Starting (l,b): $xx $yy")
                ra, dec = galactic2equatorial(xx, yy)
            elseif ref == "equatorial"
                println("## Starting (RA,Dec): $xx $yy")
                ra = xx
                dec = yy
            end
            push!(radone, ra)
            push!(decdone, dec)

            name = @sprintf("RA%.3fDec%.3f", ra, dec)
            votname = @sprintf("%s-%2.1fdeg.vot", name, radius)

            if votname in dfp.votname
                println("## $name done...")
            else
                mextra.votname = get_gaia_data_many(gaia, radius, tol, ra, dec, name, rect)
                extra(mextra, optim)

                push!(dfp, [mextra.votname])
                CSV.write(progressfile, dfp, delim=";")
                if mgene["rmvot"] == "yes"
                    rm(mextra.votname)
                    println("## votable $(mextra.votname) removed")
                end
            end
            plot_sky(radone, decdone, radius=20, figname="gridding-allsky.png")
        end
    end
end
#########################  merge function
function merge(meta)

    tstart = now()
    println(blue("## Merging catalog ..."))
    println(blue("## Starting at $tstart"))

    mmerge = meta["merge"]
    mgene = meta["general"]

    cd(mgene["wdir"])

    catalog = mmerge["catalog"]
    mergefile = name = @sprintf("%s.merge", catalog)
    mode = mmerge["mode"]

    debug_red(mergefile)

    if mode == "duplicate"
        println("### Merge, removing duplicated clusters...")
        debug_red(mmerge)

        toldeg = mmerge["toldeg"]
        toldist = mmerge["toldist"]
        tolndiff = mmerge["tolndiff"]
        metric = mmerge["metric"]
        println("### Merge, metric $metric")

        dfcat = CSV.File(catalog, delim=";") |> DataFrame

        dfmerge = rm_duplicated(dfcat, toldeg, toldist, tolndiff, metric)
        CSV.write(mergefile, dfmerge, delim=";")
        println("## Catalog $mergefile created.")
    end

    if mode == "simbad"
        println("### Merge, searching for objects with Simbad...")
        coord = pyimport("astropy.coordinates")
        Simbad = pyimport("astroquery.simbad")

        tolquery = mmerge["tolquery"]       ## arcmin radius for the object search

        dfcat = CSV.File(catalog, delim=";") |> DataFrame

        namecla = []
        for row in eachrow(dfcat)
            ra = row.ra
            dec = row.dec
            c = coord.SkyCoord(ra, dec, unit="deg")
            res = Simbad.Simbad.query_region(c, radius="$tolquery arcmin")
            namecl = "-"

            if res !== nothing
                for lobj in res
                    t = split(lobj[1])
                    if t[1] == "Cl*"
                        namecl = t[2] * "_" * t[3]
                        debug_red(namecl)
                    end
                end
            end
            push!(namecla, namecl)
        end
        dfcat[!, :name] = namecla

        CSV.write(mergefile, dfcat, delim=";")
        println("## Catalog $mergefile created.")
    end

end
#########################
function get_blacklist(m)
    if haskey(m, "blacklist")
        blackfile = m["blacklist"]
        if isfile(blackfile)
            dfblck = CSV.File(blackfile, delim=",") |> DataFrame
            println("### Blacklist $blackfile read")
        else
            dfblck = DataFrame(votname=[""])
        end
    else
        dfblck = DataFrame(votname=[""])
    end
    return (dfblck)
end
#################################### MAIN ########################
let
    # println(ARGS)
    println("############################")
    println("### Building Gaia results...")

    metabuild = TOML.parsefile(ARGS[1])

    key = collect(keys(metabuild))
    # println(key)
    for k in key
        # key1 = collect(keys(metabuild[k]))

        if k == "reprocess"
            reprocess(metabuild)
        elseif k == "fullsky"
            fullSky(metabuild)
        elseif k == "random"
            randomfields(metabuild)
        elseif k == "gridding"
            gridding(metabuild)
        elseif k == "merge"
            merge(metabuild)
        end
    end
end

