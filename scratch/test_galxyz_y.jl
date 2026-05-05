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

# Test case for l=90, b=0 (approximate RA/Dec)
# l=90, b=0 is roughly RA=317.7, Dec=48.3 (Cygnus area)
println("l=90, b=0 (dist=1000): ", galXYZ(317.7, 48.3, 1000))
