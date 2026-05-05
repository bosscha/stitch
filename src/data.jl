## from module gaiaClustering
## Functions to deal with GAIA data and to normalize
##

############
## Query the GAIA data towards the coord with a radius conesearch
#####
function query_gaia(coord, radius, dump=false)
    #####
    return (0)

end

function copy(s::Df)::Df
    c = Df(s.ndata, zeros(length(s.data[:, 1]), s.ndata), zeros(length(s.raw[:, 1]), s.ndata), zeros(length(s.err[:, 1]), s.ndata),
        zeros(length(s.sourceid[:, 1]), s.ndata))
    c.data[:, :] = s.data[:, :]
    c.raw[:, :] = s.raw[:, :]
    c.err[:, :] = s.err[:, :]
    c.sourceid[:, :] = s.sourceid[:, :]

    return (c)
end

## dummy ...
function copy1(s::Df)::Df
    c = Df(s.ndata, zeros(length(s.data[:, 1]), s.ndata), zeros(length(s.raw[:, 1]), s.ndata), zeros(length(s.err[:, 1]), s.ndata),
        zeros(length(s.sourceid[:, 1]), s.ndata))
    c.data[:, :] = s.data[:, :]
    c.raw[:, :] = s.raw[:, :]
    c.err[:, :] = s.err[:, :]
    c.sourceid[:, :] = s.sourceid[:, :]

    return (c)
end

#########
## function to create the df
######
function read_votable(voname::String)
    ######
    warnings = pyimport("warnings")
    warnings.filterwarnings("ignore")
    votable = pyimport("astropy.io.votable")
    vot = votable.parse(voname)
    data = vot.get_first_table()

    println("## Votable $voname read")

    ### return(data["array"]["data"])
    return (data.array.data)
end


#########
function filter_data(gaia, dist_range=[0., 2000], vra_range=[-250, 250],
    vdec_range=[-250., 250], mag_range=[-1e9, 1e9]; zpt=false)::Df
    ########

    ngaia = length(gaia)

    source_id = zeros(Int64, ngaia)
    lgal = zeros(ngaia)
    bgal = zeros(ngaia)
    ra = zeros(ngaia)
    dec = zeros(ngaia)
    distance = zeros(ngaia)
    pmra = zeros(ngaia)
    pmdec = zeros(ngaia)
    parallax = zeros(ngaia)
    vra = zeros(ngaia)
    vdec = zeros(ngaia)
    g = zeros(ngaia)
    rp = zeros(ngaia)
    bp = zeros(ngaia)
    parallax_error = zeros(ngaia)
    pmra_error = zeros(ngaia)
    pmdec_error = zeros(ngaia)
    radialvel = zeros(ngaia)
    ## Galactic proper motion and velocities
    pml = zeros(ngaia)
    pmb = zeros(ngaia)
    vl = zeros(ngaia)
    vb = zeros(ngaia)
    ## ZPT parameters
    nu_eff_used_in_astrometry = zeros(ngaia)
    pseudocolour = zeros(ngaia)
    ecl_lat = zeros(ngaia)
    astrometric_params_solved = zeros(ngaia)
    zcorr = zeros(ngaia)

    ## Extinction A_0, E(B-R),  iron abundance(dex)
    ag = zeros(ngaia)
    a0 = zeros(ngaia)
    ebmr = zeros(ngaia)
    mh = zeros(ngaia)

    for i in 1:ngaia

        source_id[i] = safe_pyconvert(Int64, get(gaia, i - 1).source_id, 0)
        lgal[i] = safe_pyconvert(Float64, get(gaia, i - 1).l, 0.0)
        bgal[i] = safe_pyconvert(Float64, get(gaia, i - 1).b, 0.0)
        ra[i] = safe_pyconvert(Float64, get(gaia, i - 1).ra, 0.0)
        dec[i] = safe_pyconvert(Float64, get(gaia, i - 1).dec, 0.0)
        
        parallax_val = safe_pyconvert(Float64, get(gaia, i - 1).parallax, 0.0)
        parallax[i] = parallax_val
        distance[i] = parallax_val == 0.0 ? Inf : 1000. / parallax_val
        
        pmra[i] = safe_pyconvert(Float64, get(gaia, i - 1).pmra, 0.0)
        pmdec[i] = safe_pyconvert(Float64, get(gaia, i - 1).pmdec, 0.0)
        vra[i] = 4.74e-3 * pmra[i] * distance[i]
        vdec[i] = 4.74e-3 * pmdec[i] * distance[i]

        ## Galactic Proper motions
        muG = PM_equatorial2galactic(pmra[i], pmdec[i], ra[i], dec[i], lgal[i])
        pml[i] = muG[1]
        pmb[i] = muG[2]
        vl[i] = 4.74e-3 * pml[i] * distance[i]
        vb[i] = 4.74e-3 * pmb[i] * distance[i]

        #fix for EDR3
        radialvel[i] = safe_pyconvert(Float64, get(gaia, i - 1).radial_velocity, 0.0)

        ### errors.
        parallax_error[i] = safe_pyconvert(Float64, get(gaia, i - 1).parallax_error, 0.0)
        pmra_error[i] = safe_pyconvert(Float64, get(gaia, i - 1).pmra_error, 0.0)
        pmdec_error[i] = safe_pyconvert(Float64, get(gaia, i - 1).pmdec_error, 0.0)

        g[i] = safe_pyconvert(Float64, get(gaia, i - 1).phot_g_mean_mag, 0.0)
        rp[i] = safe_pyconvert(Float64, get(gaia, i - 1).phot_rp_mean_mag, 0.0)
        bp[i] = safe_pyconvert(Float64, get(gaia, i - 1).phot_bp_mean_mag, 0.0)

        # extinction,  reddening, iron abundance
        ag[i] = safe_pyconvert(Float64, get(gaia, i - 1).ag_gspphot, NaN)
        a0[i] = safe_pyconvert(Float64, get(gaia, i - 1).azero_gspphot, NaN)
        ebmr[i] = safe_pyconvert(Float64, get(gaia, i - 1).ebpminrp_gspphot, NaN)
        mh[i] = safe_pyconvert(Float64, get(gaia, i - 1).mh_gspphot, NaN)

        ## for ZPT correction
        nu_eff_used_in_astrometry[i] = safe_pyconvert(Float64, get(gaia, i - 1).nu_eff_used_in_astrometry, 0.0)
        pseudocolour[i] = safe_pyconvert(Float64, get(gaia, i - 1).pseudocolour, 0.0)
        ecl_lat[i] = safe_pyconvert(Float64, get(gaia, i - 1).ecl_lat, 0.0)
        astrometric_params_solved[i] = safe_pyconvert(Float64, get(gaia, i - 1).astrometric_params_solved, 0.0)
    end

    if zpt
        try
            zpt_module = pyimport("zero_point.zpt")
            pyimport("warnings").filterwarnings("ignore", category=pyimport("builtins").UserWarning)
            zpt_module.load_tables()
            
            # Masking for valid astrometric solutions (31 or 95)
            # as the library raises ValueError otherwise.
            mask_zpt = [ (p == 31 || p == 95) for p in astrometric_params_solved ]
            params_zpt = Base.copy(astrometric_params_solved)
            for k in 1:length(params_zpt)
                if !mask_zpt[k] params_zpt[k] = 31 end
            end

            zcorr_raw = zpt_module.get_zpt(g, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, params_zpt)
            zcorr = zcorr_raw .* mask_zpt

            parallax = parallax .- zcorr
            distance = 1000. ./ parallax
            vra = 4.74e-3 .* pmra .* distance
            vdec = 4.74e-3 .* pmdec .* distance
            vl = 4.74e-3 .* pml .* distance
            vb = 4.74e-3 .* pmb .* distance

            debug_red("Cleaning...")
            GC.gc()
        catch e
            println("## Issues with the ZPT correction: ", e)
        end
    end

    ## Filtering ...
    i1 = distance .> dist_range[1]
    i2 = distance .< dist_range[2]
    i3 = vra .> vra_range[1]
    i4 = vra .< vra_range[2]
    i5 = vdec .> vdec_range[1]
    i6 = vdec .< vdec_range[2]
    i7 = g .> mag_range[1]
    i8 = g .< mag_range[2]
    i9 = rp .> mag_range[1]
    i10 = rp .< mag_range[2]
    i11 = bp .> mag_range[1]
    i12 = bp .< mag_range[2]

    ifinal = i1 .& i2 .& i3 .& i4 .& i5 .& i6 .& i7 .& i8 .& i9 .& i10 .& i11 .& i12

    ## G magnitude
    ##  gbar =  g[ifinal] - 5 .* log10.(distance[ifinal]) .+ 17.
    gbar = g[ifinal]

    ## Df of the filtered dat
    ndata = length(distance[ifinal])
    s = Df(ndata, zeros(8, ndata), zeros(17, ndata), zeros(8, ndata), zeros(1, ndata))

    s.data[1, :] = lgal[ifinal]
    s.data[2, :] = bgal[ifinal]
    s.data[3, :] = distance[ifinal]
    s.data[4, :] = vl[ifinal]
    s.data[5, :] = vb[ifinal]
    s.data[6, :] = gbar
    s.data[7, :] = g[ifinal] .- rp[ifinal]
    s.data[8, :] = bp[ifinal] .- g[ifinal]

    s.raw[1, :] = ra[ifinal]
    s.raw[2, :] = dec[ifinal]
    s.raw[3, :] = lgal[ifinal]
    s.raw[4, :] = bgal[ifinal]
    s.raw[5, :] = parallax[ifinal]
    s.raw[6, :] = pmra[ifinal]
    s.raw[7, :] = pmdec[ifinal]
    s.raw[8, :] = pml[ifinal]
    s.raw[9, :] = pmb[ifinal]
    s.raw[10, :] = g[ifinal]
    s.raw[11, :] = rp[ifinal]
    s.raw[12, :] = bp[ifinal]
    s.raw[13, :] = radialvel[ifinal]
    s.raw[14, :] = ag[ifinal]
    s.raw[15, :] = a0[ifinal]
    s.raw[16, :] = ebmr[ifinal]
    s.raw[17, :] = mh[ifinal]

    ## Errors ..
    s.err[1, :] = parallax_error[ifinal]
    s.err[2, :] = pmra_error[ifinal]
    s.err[3, :] = pmdec_error[ifinal]
    s.err[4, :] = zcorr[ifinal]

    ## GAIA source id
    s.sourceid[1, :] = source_id[ifinal]

    println("## Filtering done ...")
    println(yellow("## Stars selected: $ndata"))

    debug_red("Cleaning ZPT , data, etc...")
    zcorr = 0
    nu_eff_used_in_astrometry = 0
    pseudocolour = 0
    ecl_lat = 0
    astrometric_params_solved = 0
    GC.gc()

    return (s)
end

function equatorial2galactic(α, δ)
    ## NGP coordinates
    αG = 192.85948
    δG = 27.12825
    lascend = 32.93

    b = asind(cosd(δ) * cosd(δG) * cosd(α - αG) + sind(δ) * sind(δG))

    x = cosd(δ) * cosd(δG) * sind(α - αG)
    y = sind(δ) - sind(b) * sind(δG)
    l = atand(y, x) + lascend

    if x >= 0 && y <= 0
        l += 360.0
    end
    if x <= 0 && y < 0
        l += 360.0
    end
    #if x<0  && y>=0  l -= 180.0 end

    return (l, b)
end

## angle between two points on a sphere
function angle4sphere(long1, lat1, long2, lat2)
    dLon = long2 - long1
    cosang = cosd(lat1) * cosd(lat2) * cosd(dLon) + sind(lat1) * sind(lat2)
    ang = acosd(cosang)
end

## Transform PM from equatorial to galactic system.
## See Poleski 1997 / arXiv
## PM,,corr is from Conrad (2015)
function PM_equatorial2galactic(μα, μδ, α, δ, l)
    ## NGP coordinates
    αG = 192.85948
    δG = 27.12825

    C1 = sind(δG) * cosd(δ) - cosd(δG) * sind(δ) * cosd(α - αG)
    C2 = cosd(δG) * sind(α - αG)
    k = 1 / sqrt(C1^2 + C2^2)
    A = k * [C1 C2; -C2 C1]
    PMG = A * [μα; μδ]

    ## PM along gal. lat. corrected for differential velocity
    ## Oort constants. Not applied.
    #A , B = (14.5 , -13.) ./ 4.74
    # PMG[1] =  PMG[1] - (A*cosd(2l) + B)

    return (PMG)
end

## Compute the X,Y,Z galactic coordinates (centered on the Galactic Center)
## See Ellsworth-Bowers et al. (2013)
## Rgal was updated from Anderson et al. (2018)
## xg,yg,zg in pc
## Modified to be left-handed (yg is modified to -yg from original)
## see http://astro.utoronto.ca/~bovy/AST1420/notes/notebooks/A.-Coordinate-systems.html
##

function galXYZ(ra, dec, distance)
    # debug_red("CHANGED l b to RA Dec!!!!!!!!!!!!!!")
    # debug_red("$ra $dec $distance")

    Rgal = 8.23e3   # (Leung et al. 2022)
    zsun = 20.8     #  (Bennett & Bovy 2019)
    θ = asin(zsun / Rgal)
    η = 58.5986320306   # degrees

    raGC = 266.40506655
    decGC = -28.93616241    ## Galactic center sky coordinates
    ricrs = [distance * cosd(ra) * cosd(dec) distance * sind(ra) * cosd(dec) distance * sind(dec)]'

    xGC = [Rgal 0 0]'    # sun coordinates in the Galaxy

    # formulas from https://docs.astropy.org/en/stable/coordinates/galactocentric.html

    H = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]      # θ in radians
    R1 = [cosd(decGC) 0 sind(decGC); 0 1 0; -sind(decGC) 0 cosd(decGC)]
    R2 = [cosd(raGC) sind(raGC) 0; -sind(raGC) cosd(raGC) 0; 0 0 1]
    R3 = [1 0 0; 0 cosd(η) sind(η); 0 -sind(η) cosd(η)]

    R = R3 * R1 * R2

    rfull = R * ricrs - xGC
    rfull = H * rfull

    # debug_red(rfull)
    xg = rfull[1]
    yg = rfull[2]
    zg = rfull[3]

    # debug_red("$xg $yg $zg")

    return (xg, yg, zg)
end

### Correction of the radial velocity
### Conrad (3025)
function RVEL_corr(rvel, distance, l)
    ## Oort's constant
    A = 14.5
    rv = rvel - A * 1e-3 * distance * sind(2l)
    return (rv)
end

### compute galactic U V W from idlastro gal_uvw.pro
##
## ra, dec: degrees
## distance; pc
## pmra,pmdec: milliarcsec/yr
## vrad: km/s

## UVW: km/s
#      U - Velocity (km/s) positive toward the Galactic *anti*center
#      V - Velocity (km/s) positive in the direction of Galactic rotation
#      W - Velocity (km/s) positive toward the North Galactic Pole
#
function galUVW(ra, dec, distance, pmra, pmdec, vrad; LSR_vel=[-8.5; 13.38; 6.49])
    k = 4.74047     #Equivalent of 1 A.U/yr in km/s

    T = [0.0548756 0.873437 0.483835;
        0.494109 -0.44483 0.746982;
        -0.867666 -0.198076 0.455984]

    A1 = [cosd(ra) sind(ra) 0; sind(ra) -cosd(ra) 0; 0 0 -1]
    A2 = [cosd(dec) 0 -sind(dec); 0 -1 0; -sind(dec) 0 -cosd(dec)]
    vec1 = vrad
    vec2 = k * pmra * 1e-3 * distance
    vec3 = k * pmdec * 1e-3 * distance
    v = [vec1; vec2; vec3]

    uvw = T * A1 * A2 * v + LSR_vel

    return (uvw)
end


######
function add_cartesian(s::Df, centering=true)::Df
    ######
    dfresult = copy(s)
    off = zeros(2)

    if centering
        off[1] = mean(s.data[1, :])
        off[2] = mean(s.data[2, :])
    end

    lgal = DEG2RAD .* (s.data[1, :] .- off[1])       ## to be checked
    bgal = DEG2RAD .* (s.data[2, :] .- off[2])

    dfresult.data[1, :] = s.data[3, :] .* cos.(bgal) .* cos.(lgal)
    dfresult.data[2, :] = s.data[3, :] .* cos.(bgal) .* sin.(lgal)
    dfresult.data[3, :] = s.data[3, :] .* sin.(bgal)

    println("## Cartesian transformation done ...")

    return (dfresult)
end


######
function normalization_PerBlock(s::Df, block, weightblock, norm, density=false, verbose=true)
    ######
    dfresult = copy(s)
    ndf = size(s.data)
    scale8d = zeros(ndf[1])
    vector8d = 0.

    ind = 1
    for aw in zip(block, weightblock)
        weight = aw[2]
        for ak in aw[1]
            normK = normalizationVector(norm, density, dfresult.data[ak, :])
            # normK[2] = normK[2] * totalWeight
            dfresult.data[ak, :] = weight .* (s.data[ak, :] .- normK[1]) ./ normK[2]
            scale8d[ind] = weight / normK[2]
            vector8d += scale8d[ind]^2
            ind += 1
        end
    end

    vector8d = sqrt(vector8d)
    scale8d[:] = scale8d[:] ./ vector8d
    dfresult.data[:, :] = dfresult.data[:, :] ./ vector8d

    if verbose
        println("## Normalization $norm done...")
        println("### [1pc,1pc,1pc,1km/s,1km/s,1mag,1mag,1mag] equivalent to $scale8d")
        println("##")
    end

    return (dfresult, scale8d)
end

######
function normalizationVector(norm, density, arr)
    ######
    vecNorm = [0.0, 1.0]

    if norm == "identity"
        vecNorm = [0.0, 1.0]

    elseif norm == "normal"
        stdArr = std(arr)
        meanArr = mean(arr)
        vecNorm = [meanArr, stdArr]

    elseif norm == "minmax"
        minarr = minimum(vcat(arr...))
        maxarr = maximum(vcat(arr...))
        vecNorm = [minarr, maxarr - minarr]

    end

    if density
        vecNorm[2] = vecNorm[2] * length(arr)
    end

    return (vecNorm)
end

######
function subsetDf(df::Df, indx)::Df
    ######

    ndat = length(indx)
    subset = Df(ndat, df.data[:, indx], df.raw[:, indx], df.err[:, indx], df.sourceid[:, indx])

    return (subset)
end


####
# Create the DataFrame to save the cluster...
##
function export_df(votname, ocdir, df, dfcart, labels, labelmax, pc, m::GaiaClustering.meta; save=true, cluster_id="")
    ra = df.raw[1, labels[labelmax]]
    dec = df.raw[2, labels[labelmax]]
    l = df.data[1, labels[labelmax]]
    b = df.data[2, labels[labelmax]]
    parallax = df.raw[5, labels[labelmax]]
    d = df.data[3, labels[labelmax]]
    pmra = df.raw[6, labels[labelmax]]
    pmdec = df.raw[7, labels[labelmax]]
    X = dfcart.data[1, labels[labelmax]]
    Y = dfcart.data[2, labels[labelmax]]
    Z = dfcart.data[3, labels[labelmax]]
    vl = df.data[4, labels[labelmax]]
    vb = df.data[5, labels[labelmax]]
    vrad = df.raw[13, labels[labelmax]]
    gbar = df.raw[10, labels[labelmax]]
    rp = df.raw[11, labels[labelmax]]
    bp = df.raw[12, labels[labelmax]]
    ag = df.raw[14, labels[labelmax]]
    a0 = df.raw[15, labels[labelmax]]
    ebmr = df.raw[16, labels[labelmax]]
    mh = df.raw[17, labels[labelmax]]
    parallax_err = df.err[1, labels[labelmax]]

    source_id = df.sourceid[1, labels[labelmax]]


    #####
    maxy = maximum(Y)
    maxz = maximum(Z)
    #####

    # galactocentric coordinates...
    s = size(ra)
    xg = zeros(s[1])
    yg = zeros(s[1])
    zg = zeros(s[1])

    ind = 1
    for i in 1:s[1]
        xx, yy, zz = galXYZ(ra[i], dec[i], d[i])
        xg[i] = xx
        yg[i] = yy
        zg[i] = zz
    end

    oc = DataFrame(sourceid=source_id, ra=ra, dec=dec, l=l, b=b, parallax=parallax, parallax_err=parallax_err, distance=d,
        pmra=pmra, pmdec=pmdec, X=X, Y=Y, Z=Z, vl=vl,
        vb=vb, vrad=vrad, Xg=xg, Yg=yg, Zg=zg, gbar=gbar, rp=rp, bp=bp, ag=ag, a0=a0, ebmr=ebmr, mh=mh)

    # spc = size(pc)
    # if m.pca == "yes" && length(spc) >= 2 && s[1] == spc[2]
    #     for i in 1:spc[1]
    #         colname = "PC$i"
    #         oc[!, colname] = pc[i, :]
    #     end
    # else
    #     println("### Warning : PCA size and star number are not equal....")
    # end

    ## add type core(1)/tail(2)
    st = size(oc[!, :ra])
    type = ones(st[1])  ## default to 1
    if m.tail == "yes"
        oc[!, :type] = type
        for i in 1:st[1]
            if oc[i, :sourceid] in df.sourceid[1, labels[3]]
                oc[i, :type] = 2    ### core==step2
            end
        end
    else
        ### core==step1
        oc[!, :type] = type
    end

    oc[!, :type] = convert.(Int8, oc[!, :type])
    if cluster_id != ""
        oc[!, :cluster_id] .= cluster_id
    else
        oc[!, :cluster_id] .= votname
    end

    name = split(votname, ".")
    infix = ""
    for iname in name
        if iname != "vot"
            infix *= iname * "."
        end
    end
    infix *= "oc.csv"
    filename = @sprintf("%s/%s", ocdir, infix)
    if save
        if m.savedb == "yes"
            println("### Saving cluster members into PostgreSQL table $(m.dbtable)...")

            # Construct connection string
            # Handle empty password string properly for Postgres
            pwd_arg = (m.dbpass != "") ? "password=$(m.dbpass)" : ""
            conn_str = "host=$(m.dbhost) user=$(m.dbuser) dbname=$(m.dbname) $pwd_arg"

            conn = LibPQ.Connection(conn_str)
            try
                execute(
                    conn,
                    """
      CREATE TABLE IF NOT EXISTS $(m.dbtable) (
          sourceid BIGINT, ra FLOAT, dec FLOAT, l FLOAT, b FLOAT, parallax FLOAT, parallax_err FLOAT, distance FLOAT,
          pmra FLOAT, pmdec FLOAT, X FLOAT, Y FLOAT, Z FLOAT, vl FLOAT, vb FLOAT, vrad FLOAT, Xg FLOAT, Yg FLOAT, Zg FLOAT,
          gbar FLOAT, rp FLOAT, bp FLOAT, ag FLOAT, a0 FLOAT, ebmr FLOAT, mh FLOAT, type SMALLINT, cluster_id VARCHAR(255)
      )
      """
                )

                # Make sure the column exists incase the table was already created
                execute(conn, "ALTER TABLE $(m.dbtable) ADD COLUMN IF NOT EXISTS cluster_id VARCHAR(255);")

                # if m.pca == "yes" && s[1] == spc[2]
                #     for i in 1:spc[1]
                #         execute(conn, "ALTER TABLE $(m.dbtable) ADD COLUMN IF NOT EXISTS PC$i FLOAT;")
                #     end
                # end

                ncol = size(oc)[2]
                val_placeholders = join(["\$$i" for i in 1:ncol], ", ")
                col_names = join(names(oc), ", ")
                LibPQ.load!(oc, conn, "INSERT INTO $(m.dbtable) ($col_names) VALUES ($val_placeholders);")
                @printf("### Database insertion completed for %d rows in table %s \n", size(oc)[1], m.dbtable)
            catch e
                println("### DB ERROR: ", e)
            finally
                close(conn)
            end
        end

        if isdir(ocdir)
            CSV.write(filename, oc, delim=';')
            @printf("### %s created  in %s \n", filename, ocdir)
        else
            print("\n### Error, result (oc) directory $ocdir not found... \n")
            exit()
        end
    end

    return (oc)
end

#######################################
## export cluster metadata to PostgreSQL
#######################################
function export_sc_db(sc::DataFrame, m::GaiaClustering.meta)
    table_name = "$(m.dbtable)_metadata"
    println("### Saving cluster metadata into PostgreSQL table $table_name...")

    pwd_arg = (m.dbpass != "") ? "password=$(m.dbpass)" : ""
    conn_str = "host=$(m.dbhost) user=$(m.dbuser) dbname=$(m.dbname) $pwd_arg"

    conn = LibPQ.Connection(conn_str)
    try
        col_defs = String[]
        for col in propertynames(sc)
            T = eltype(sc[!, col])
            type_str = ""
            if T <: Integer
                type_str = "BIGINT"
            elseif T <: AbstractFloat
                type_str = "FLOAT"
            else
                type_str = "VARCHAR(255)"
            end
            push!(col_defs, "$col $type_str")
        end

        create_sql = "CREATE TABLE IF NOT EXISTS $table_name (" * join(col_defs, ", ") * ");"
        execute(conn, create_sql)

        for coldef in col_defs
            execute(conn, "ALTER TABLE $table_name ADD COLUMN IF NOT EXISTS $coldef;")
        end

        ncol = size(sc)[2]
        val_placeholders = join(["\$$i" for i in 1:ncol], ", ")
        col_names = join(names(sc), ", ")

        LibPQ.load!(sc, conn, "INSERT INTO $table_name ($col_names) VALUES ($val_placeholders);")

        println("### Cluster metadata insertion completed for table $table_name")
    catch e
        println("### DB ERROR in export_sc_db: ", e)
    finally
        close(conn)
    end
end
#######################################
## a built-in version of getdata
function get_data(m::GaiaClustering.meta)
    println("## Distance cut : $(m.mindist) $(m.maxdist) pc")

    if m.zpt == "yes"
        zoff = true
        println("## Applying Zero Point offset correction on parallax...")
    else
        zoff = false
    end

    data = read_votable(m.votdir * "/" * m.votname)
    df = filter_data(data, [m.mindist, m.maxdist], zpt=zoff)
    dfcart = add_cartesian(df)
    blck = [[1, 2, 3], [4, 5], [6, 7, 8]]
    wghtblck = [4.0, 5.0, 1.0]
    norm = "identity"

    dfcartnorm, scale8 = normalization_PerBlock(dfcart, blck, wghtblck, norm, false)

    return (df, dfcart, dfcartnorm)
end

#######################################
## a postgresql version of getdata
function filter_pg_data(df_pg::DataFrame, dist_range, vra_range, vdec_range, mag_range; zpt=false)
    ngaia = size(df_pg)[1]

    source_id = zeros(Int64, ngaia)
    lgal = zeros(ngaia)
    bgal = zeros(ngaia)
    ra = zeros(ngaia)
    dec = zeros(ngaia)
    distance = zeros(ngaia)
    pmra = zeros(ngaia)
    pmdec = zeros(ngaia)
    parallax = zeros(ngaia)
    vra = zeros(ngaia)
    vdec = zeros(ngaia)
    g = zeros(ngaia)
    rp = zeros(ngaia)
    bp = zeros(ngaia)
    parallax_error = zeros(ngaia)
    pmra_error = zeros(ngaia)
    pmdec_error = zeros(ngaia)
    radialvel = zeros(ngaia)
    pml = zeros(ngaia)
    pmb = zeros(ngaia)
    vl = zeros(ngaia)
    vb = zeros(ngaia)
    nu_eff_used_in_astrometry = zeros(ngaia)
    pseudocolour = zeros(ngaia)
    ecl_lat = zeros(ngaia)
    astrometric_params_solved = zeros(ngaia)
    zcorr = zeros(ngaia)

    ag = zeros(ngaia)
    a0 = zeros(ngaia)
    ebmr = zeros(ngaia)
    mh = zeros(ngaia)

    has_pmra_error = hasproperty(df_pg, :pmra_error)
    has_pmdec_error = hasproperty(df_pg, :pmdec_error)

    for i in 1:ngaia
        source_id[i] = ismissing(df_pg.source_id[i]) ? 0 : df_pg.source_id[i]
        lgal[i] = ismissing(df_pg.l[i]) ? 0.0 : df_pg.l[i]
        bgal[i] = ismissing(df_pg.b[i]) ? 0.0 : df_pg.b[i]
        ra[i] = ismissing(df_pg.ra[i]) ? 0.0 : df_pg.ra[i]
        dec[i] = ismissing(df_pg.dec[i]) ? 0.0 : df_pg.dec[i]
        parallax[i] = ismissing(df_pg.parallax[i]) ? 0.0 : df_pg.parallax[i]
        distance[i] = parallax[i] == 0.0 ? Inf : 1000. / parallax[i]
        pmra[i] = ismissing(df_pg.pmra[i]) ? 0.0 : df_pg.pmra[i]
        pmdec[i] = ismissing(df_pg.pmdec[i]) ? 0.0 : df_pg.pmdec[i]
        vra[i] = 4.74e-3 * pmra[i] * distance[i]
        vdec[i] = 4.74e-3 * pmdec[i] * distance[i]

        muG = GaiaClustering.PM_equatorial2galactic(pmra[i], pmdec[i], ra[i], dec[i], lgal[i])
        pml[i] = muG[1]
        pmb[i] = muG[2]
        vl[i] = 4.74e-3 * pml[i] * distance[i]
        vb[i] = 4.74e-3 * pmb[i] * distance[i]

        radialvel[i] = ismissing(df_pg.radial_velocity[i]) ? 0.0 : df_pg.radial_velocity[i]

        parallax_error[i] = ismissing(df_pg.parallax_error[i]) ? 0.0 : df_pg.parallax_error[i]
        pmra_error[i] = has_pmra_error && !ismissing(df_pg.pmra_error[i]) ? df_pg.pmra_error[i] : 0.0
        pmdec_error[i] = has_pmdec_error && !ismissing(df_pg.pmdec_error[i]) ? df_pg.pmdec_error[i] : 0.0

        g[i] = ismissing(df_pg.phot_g_mean_mag[i]) ? 0.0 : df_pg.phot_g_mean_mag[i]
        rp[i] = ismissing(df_pg.phot_rp_mean_mag[i]) ? 0.0 : df_pg.phot_rp_mean_mag[i]
        bp[i] = ismissing(df_pg.phot_bp_mean_mag[i]) ? 0.0 : df_pg.phot_bp_mean_mag[i]

        ag[i] = ismissing(df_pg.ag_gspphot[i]) ? NaN : df_pg.ag_gspphot[i]
        a0[i] = ismissing(df_pg.azero_gspphot[i]) ? NaN : df_pg.azero_gspphot[i]
        ebmr[i] = ismissing(df_pg.ebpminrp_gspphot[i]) ? NaN : df_pg.ebpminrp_gspphot[i]
        mh[i] = ismissing(df_pg.mh_gspphot[i]) ? NaN : df_pg.mh_gspphot[i]

        nu_eff_used_in_astrometry[i] = ismissing(df_pg.nu_eff_used_in_astrometry[i]) ? 0.0 : df_pg.nu_eff_used_in_astrometry[i]
        pseudocolour[i] = ismissing(df_pg.pseudocolour[i]) ? 0.0 : df_pg.pseudocolour[i]
        ecl_lat[i] = ismissing(df_pg.ecl_lat[i]) ? 0.0 : df_pg.ecl_lat[i]
        astrometric_params_solved[i] = ismissing(df_pg.astrometric_params_solved[i]) ? 0 : df_pg.astrometric_params_solved[i]
    end

    if zpt
        try
            zpt_module = pyimport("zero_point.zpt")
            pyimport("warnings").filterwarnings("ignore", category=pyimport("builtins").UserWarning)
            zpt_module.load_tables()
            
            # Masking for valid astrometric solutions (31 or 95)
            mask_zpt = [ (p == 31 || p == 95) for p in astrometric_params_solved ]
            params_zpt = Base.copy(astrometric_params_solved)
            for k in 1:length(params_zpt)
                if !mask_zpt[k] params_zpt[k] = 31 end
            end

            zcorr_raw = zpt_module.get_zpt(g, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, params_zpt)
            zcorr = zcorr_raw .* mask_zpt


            parallax = parallax .- zcorr
            distance = 1000. ./ parallax
            vra = 4.74e-3 .* pmra .* distance
            vdec = 4.74e-3 .* pmdec .* distance
            vl = 4.74e-3 .* pml .* distance
            vb = 4.74e-3 .* pmb .* distance

            GC.gc()
        catch e
            println("## Issues with the ZPT correction: ", e)
        end
    end

    i1 = distance .> dist_range[1]
    i2 = distance .< dist_range[2]
    i3 = vra .> vra_range[1]
    i4 = vra .< vra_range[2]
    i5 = vdec .> vdec_range[1]
    i6 = vdec .< vdec_range[2]
    i7 = g .> mag_range[1]
    i8 = g .< mag_range[2]
    i9 = rp .> mag_range[1]
    i10 = rp .< mag_range[2]
    i11 = bp .> mag_range[1]
    i12 = bp .< mag_range[2]

    ifinal = i1 .& i2 .& i3 .& i4 .& i5 .& i6 .& i7 .& i8 .& i9 .& i10 .& i11 .& i12

    gbar = g[ifinal]

    ndata = length(distance[ifinal])
    s = GaiaClustering.Df(ndata, zeros(8, ndata), zeros(17, ndata), zeros(8, ndata), zeros(1, ndata))

    s.data[1, :] = lgal[ifinal]
    s.data[2, :] = bgal[ifinal]
    s.data[3, :] = distance[ifinal]
    s.data[4, :] = vl[ifinal]
    s.data[5, :] = vb[ifinal]
    s.data[6, :] = gbar
    s.data[7, :] = g[ifinal] .- rp[ifinal]
    s.data[8, :] = bp[ifinal] .- g[ifinal]

    s.raw[1, :] = ra[ifinal]
    s.raw[2, :] = dec[ifinal]
    s.raw[3, :] = lgal[ifinal]
    s.raw[4, :] = bgal[ifinal]
    s.raw[5, :] = parallax[ifinal]
    s.raw[6, :] = pmra[ifinal]
    s.raw[7, :] = pmdec[ifinal]
    s.raw[8, :] = pml[ifinal]
    s.raw[9, :] = pmb[ifinal]
    s.raw[10, :] = g[ifinal]
    s.raw[11, :] = rp[ifinal]
    s.raw[12, :] = bp[ifinal]
    s.raw[13, :] = radialvel[ifinal]
    s.raw[14, :] = ag[ifinal]
    s.raw[15, :] = a0[ifinal]
    s.raw[16, :] = ebmr[ifinal]
    s.raw[17, :] = mh[ifinal]

    s.err[1, :] = parallax_error[ifinal]
    s.err[2, :] = pmra_error[ifinal]
    s.err[3, :] = pmdec_error[ifinal]
    s.err[4, :] = zcorr[ifinal]

    s.sourceid[1, :] = source_id[ifinal]

    return (s)
end

function get_data_pg(m::GaiaClustering.meta, pixels::Vector{Int})
    println("## Fetching HEALPix pixels $pixels from postgres database...")

    pwd_arg = (m.dbpass != "") ? "password=$(m.dbpass)" : ""
    conn_str = "host=$(m.dbhost) user=$(m.dbuser) dbname=$(m.dbname) $pwd_arg"

    conn = LibPQ.Connection(conn_str)

    pixels_str = join(pixels, ", ")
    query = """
    SELECT * FROM gaia_source 
    WHERE (source_id >> 49) IN ($pixels_str)
    """

    result = execute(conn, query)
    df_pg = DataFrame(result)
    close(conn)

    println("## Retrieved $(size(df_pg)[1]) stars.")

    if size(df_pg)[1] == 0
        error("No data found for pixels $pixels")
    end

    println("## Distance cut : $(m.mindist) $(m.maxdist) pc")

    if m.zpt == "yes"
        zoff = true
        println("## Applying Zero Point offset correction on parallax...")
    else
        zoff = false
    end

    df = filter_pg_data(df_pg, [m.mindist, m.maxdist], [m.minvra, m.maxvra], [m.minvdec, m.maxvdec], [m.ming, m.maxg], zpt=zoff)
    dfcart = add_cartesian(df)
    blck = [[1, 2, 3], [4, 5], [6, 7, 8]]
    wghtblck = [4.0, 5.0, 1.0]
    norm = "identity"

    dfcartnorm, scale8 = normalization_PerBlock(dfcart, blck, wghtblck, norm, false)

    return (df, dfcart, dfcartnorm)
end

