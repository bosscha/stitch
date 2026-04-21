-- Create necessary tables for Gaia clustering

CREATE TABLE IF NOT EXISTS gaia_source (
    source_id BIGINT,
    ra DOUBLE PRECISION,
    dec DOUBLE PRECISION,
    parallax DOUBLE PRECISION,
    parallax_error DOUBLE PRECISION,
    pmra DOUBLE PRECISION,
    pmdec DOUBLE PRECISION,
    l DOUBLE PRECISION,
    b DOUBLE PRECISION,
    phot_g_mean_mag DOUBLE PRECISION,
    phot_bp_mean_mag DOUBLE PRECISION,
    phot_rp_mean_mag DOUBLE PRECISION,
    phot_g_mean_flux DOUBLE PRECISION,
    phot_g_mean_flux_error DOUBLE PRECISION,
    phot_bp_mean_flux DOUBLE PRECISION,
    phot_bp_mean_flux_error DOUBLE PRECISION,
    phot_rp_mean_flux DOUBLE PRECISION,
    phot_rp_mean_flux_error DOUBLE PRECISION,
    ruwe DOUBLE PRECISION,
    astrometric_excess_noise DOUBLE PRECISION,
    astrometric_params_solved BIGINT,
    radial_velocity DOUBLE PRECISION,
    radial_velocity_error DOUBLE PRECISION,
    ag_gspphot DOUBLE PRECISION,
    azero_gspphot DOUBLE PRECISION,
    ebpminrp_gspphot DOUBLE PRECISION,
    mh_gspphot DOUBLE PRECISION,
    nu_eff_used_in_astrometry DOUBLE PRECISION,
    pseudocolour DOUBLE PRECISION,
    ecl_lat DOUBLE PRECISION
);

CREATE TABLE IF NOT EXISTS import_log (
    file_name VARCHAR(255)
);

CREATE TABLE IF NOT EXISTS clusters (
    sourceid BIGINT, 
    ra FLOAT, 
    dec FLOAT, 
    l FLOAT, 
    b FLOAT, 
    parallax FLOAT, 
    parallax_err FLOAT, 
    distance FLOAT,
    pmra FLOAT, 
    pmdec FLOAT, 
    X FLOAT, 
    Y FLOAT, 
    Z FLOAT, 
    vl FLOAT, 
    vb FLOAT, 
    vrad FLOAT, 
    Xg FLOAT, 
    Yg FLOAT, 
    Zg FLOAT,
    gbar FLOAT, 
    rp FLOAT, 
    bp FLOAT, 
    ag FLOAT, 
    a0 FLOAT, 
    ebmr FLOAT, 
    mh FLOAT, 
    type SMALLINT, 
    cluster_id VARCHAR(255)
);
