using LinearAlgebra

function galXYZ(ra, dec, distance)
    Rgal = 8.23e3   # (Leung et al. 2022)
    zsun = 20.8     #  (Bennett & Bovy 2019)
    θ = asin(zsun / Rgal)
    η = 58.5986320306   # degrees

    raGC = 266.40506655
    decGC = -28.93616241    ## Galactic center sky coordinates
    
    # cosd/sind replacement
    cosd(x) = cos(deg2rad(x))
    sind(x) = sin(deg2rad(x))

    ricrs = [distance * cosd(ra) * cosd(dec), distance * sind(ra) * cosd(dec), distance * sind(dec)]

    xGC = [Rgal, 0, 0]    # sun coordinates in the Galaxy

    H = [cos(θ) 0 sin(θ); 0 1 0; -sin(θ) 0 cos(θ)]      # θ in radians
    R1 = [cosd(decGC) 0 sind(decGC); 0 1 0; -sind(decGC) 0 cosd(decGC)]
    R2 = [cosd(raGC) sind(raGC) 0; -sind(raGC) cosd(raGC) 0; 0 0 1]
    R3 = [1 0 0; 0 cosd(η) sind(η); 0 -sind(η) cosd(η)]

    R = R3 * R1 * R2

    rfull = R * ricrs - xGC
    rfull = H * rfull

    return rfull
end

# Test cases
println("Sun (dist=0): ", galXYZ(0, 0, 0))
println("GC (ra=266.4, dec=-28.9, dist=8230): ", galXYZ(266.40506655, -28.93616241, 8230))
println("NGP (ra=192.859, dec=27.128, dist=1000): ", galXYZ(192.85948, 27.12825, 1000))
