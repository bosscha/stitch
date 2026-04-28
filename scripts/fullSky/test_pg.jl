using DataFrames, LibPQ
using PyCall

rootdir = ENV["GAIA_ROOT"]
push!(LOAD_PATH, "$rootdir/run/src")
using GaiaClustering

# we need the same filtering logic but from a DataFrame
function filter_pg_data(df_pg::DataFrame, dist_range, vra_range, vdec_range, mag_range ; zpt=false)
    ngaia = size(df_pg)[1]
    
    source_id = zeros(Int64,ngaia)
    lgal = zeros(ngaia)
    bgal = zeros(ngaia)
    ra = zeros(ngaia)
    dec= zeros(ngaia)
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
    pmra_error     = zeros(ngaia)
    pmdec_error    = zeros(ngaia)
    radialvel      = zeros(ngaia)
    pml = zeros(ngaia)
    pmb = zeros(ngaia)
    vl = zeros(ngaia)
    vb = zeros(ngaia)
    nu_eff_used_in_astrometry= zeros(ngaia)
    pseudocolour= zeros(ngaia)
    ecl_lat= zeros(ngaia)
    astrometric_params_solved= zeros(ngaia)
    zcorr= zeros(ngaia)

    ag= zeros(ngaia)
    a0= zeros(ngaia)
    ebmr= zeros(ngaia)
    mh= zeros(ngaia)

    for i in 1:ngaia
        source_id[i]= df_pg.source_id[i]
        lgal[i]     = ismissing(df_pg.l[i]) ? 0.0 : df_pg.l[i]
        bgal[i]     = ismissing(df_pg.b[i]) ? 0.0 : df_pg.b[i]
        ra[i]       = df_pg.ra[i]
        dec[i]      = df_pg.dec[i]
        parallax[i] = df_pg.parallax[i]
        distance[i] = 1000. / parallax[i]
        pmra[i]     = df_pg.pmra[i]
        pmdec[i]    = df_pg.pmdec[i]
        vra[i]      = 4.74e-3 * pmra[i]  * distance[i]
        vdec[i]     = 4.74e-3 * pmdec[i] * distance[i]

        muG = GaiaClustering.PM_equatorial2galactic(pmra[i],pmdec[i]  , ra[i] , dec[i] , lgal[i])
        pml[i] = muG[1]
        pmb[i] = muG[2]
        vl[i]  = 4.74e-3 * pml[i]  * distance[i]
        vb[i]  = 4.74e-3 * pmb[i]  * distance[i]

        radialvel[i]    = ismissing(df_pg.radial_velocity[i]) ? 0.0 : df_pg.radial_velocity[i]

        parallax_error[i]  = df_pg.parallax_error[i]
        pmra_error[i]  = ismissing(df_pg.pmra_error[i]) ? 0.0 : df_pg.pmra_error[i]
        pmdec_error[i] = ismissing(df_pg.pmdec_error[i]) ? 0.0 : df_pg.pmdec_error[i]

        g[i]        = ismissing(df_pg.phot_g_mean_mag[i]) ? 0.0 : df_pg.phot_g_mean_mag[i]
        rp[i]       = ismissing(df_pg.phot_rp_mean_mag[i]) ? 0.0 : df_pg.phot_rp_mean_mag[i]
        bp[i]       = ismissing(df_pg.phot_bp_mean_mag[i]) ? 0.0 : df_pg.phot_bp_mean_mag[i]

        ag[i]       = ismissing(df_pg.ag_gspphot[i]) ? NaN : df_pg.ag_gspphot[i]
        a0[i]       = ismissing(df_pg.azero_gspphot[i]) ? NaN : df_pg.azero_gspphot[i]
        ebmr[i]     = ismissing(df_pg.ebpminrp_gspphot[i]) ? NaN : df_pg.ebpminrp_gspphot[i]
        mh[i]       = ismissing(df_pg.mh_gspphot[i]) ? NaN : df_pg.mh_gspphot[i]

        nu_eff_used_in_astrometry[i] = ismissing(df_pg.nu_eff_used_in_astrometry[i]) ? 0.0 : df_pg.nu_eff_used_in_astrometry[i]
        pseudocolour[i]= ismissing(df_pg.pseudocolour[i]) ? 0.0 : df_pg.pseudocolour[i]
        ecl_lat[i]= ismissing(df_pg.ecl_lat[i]) ? 0.0 : df_pg.ecl_lat[i]
        astrometric_params_solved[i]= ismissing(df_pg.astrometric_params_solved[i]) ? 0 : df_pg.astrometric_params_solved[i]
    end

    if zpt
        try
            zpt_module = pyimport("zero_point.zpt")
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
            vra      = 4.74e-3 .* pmra  .* distance
            vdec     = 4.74e-3 .* pmdec .* distance
            vl       = 4.74e-3 .* pml  .* distance
            vb       = 4.74e-3 .* pmb  .* distance
            
            GC.gc()
        catch e
            println("## Issues with the ZPT correction: ", e)
        end
    end

    i1 =  distance .> dist_range[1]
    i2 =  distance .< dist_range[2]
    i3 =  vra .> vra_range[1]
    i4 =  vra .< vra_range[2]
    i5 =  vdec .> vdec_range[1]
    i6 =  vdec .< vdec_range[2]
    i7 = g  .> mag_range[1]
    i8 = g  .< mag_range[2]
    i9  = rp  .> mag_range[1]
    i10 = rp  .< mag_range[2]
    i11  = bp  .> mag_range[1]
    i12  = bp  .< mag_range[2]

    ifinal = i1 .& i2 .& i3 .& i4 .& i5 .& i6 .& i7 .& i8 .& i9 .& i10 .& i11 .& i12

    gbar = g[ifinal]

    ndata = length(distance[ifinal])
    s = GaiaClustering.Df(ndata, zeros(8,ndata), zeros(17,ndata) , zeros(8,ndata) , zeros(1,ndata))

    s.data[1,:] = lgal[ifinal]
    s.data[2,:] = bgal[ifinal]
    s.data[3,:] = distance[ifinal]
    s.data[4,:] = vl[ifinal]
    s.data[5,:] = vb[ifinal]
    s.data[6,:] = gbar
    s.data[7,:] = g[ifinal] .- rp[ifinal]
    s.data[8,:] = bp[ifinal] .- g[ifinal]

    s.raw[1,:] = ra[ifinal]
    s.raw[2,:] = dec[ifinal]
    s.raw[3,:] = lgal[ifinal]
    s.raw[4,:] = bgal[ifinal]
    s.raw[5,:] = parallax[ifinal]
    s.raw[6,:] = pmra[ifinal]
    s.raw[7,:] = pmdec[ifinal]
    s.raw[8,:] = pml[ifinal]
    s.raw[9,:] = pmb[ifinal]
    s.raw[10,:] = g[ifinal]
    s.raw[11,:] = rp[ifinal]
    s.raw[12,:] = bp[ifinal]
    s.raw[13,:] = radialvel[ifinal]
    s.raw[14,:] = ag[ifinal]
    s.raw[15,:] = a0[ifinal]
    s.raw[16,:] = ebmr[ifinal]
    s.raw[17,:] = mh[ifinal]

    s.err[1,:] = parallax_error[ifinal]
    s.err[2,:] = pmra_error[ifinal]
    s.err[3,:] = pmdec_error[ifinal]
    s.err[4,:] = zcorr[ifinal]

    s.sourceid[1,:] = source_id[ifinal]

    return(s)
end

println("Successfully built function filter_pg_data")
