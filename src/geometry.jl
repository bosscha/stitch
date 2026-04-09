### Function to analyze the geometry
###
function voronoi_perimeter(ver, region)
    # compute the perimetet  in the voronoi region and vertices

    perimeter = 0.
    sizeRegion = length(region)

    if -1 in region
        perimeter = 0.
    else
        for index in 1:sizeRegion
            index2 = (index%sizeRegion) + 1

            if region[index]!= -1 || region[index2] != -1

                x0 = ver[region[index]+1,1]
                y0 = ver[region[index]+1,2]

                x1 = ver[region[index2]+1,1]
                y1 = ver[region[index2]+1,2]
                perimeter += sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0))
            end
        end
    end
    return(perimeter)
end
##
function voronoi_area(ver, region)
    # compute the area   in the voronoi region and vertices

    area = 0.
    sizeRegion = length(region)

    if -1 in region
        area = 0.
    else
        for index in 1:sizeRegion
            x0 = ver[region[index]+1,1]
            y0 = ver[region[index]+1,2]

            index2 = (index%sizeRegion) + 1
            x1 = ver[region[index2]+1,1]
            y1 = ver[region[index2]+1,2]

            area += 0.5 * (x0*y1 - x1*y0)
        end
    end
    return(abs(area))
end
## Voronoi tesselation using the scipy function
##
function voronoi_python(pts, verbose = true)
    let
        ndat = length(pts)
        peri = zeros(ndat)
        area = zeros(ndat)

        vor = 0
        spatial = pyimport("scipy.spatial")

        try
            # vor = spatial[:Voronoi](pts)
            vor = spatial.Voronoi(pts)
        catch
            println("## Voronoi error...")
            return([1000,1000],[1000,1000])   ## arbitrary values for peri and area
        end

        if verbose println("## Voronoi tesselation done.") end

        ver = vor.vertices
        reg = vor.regions
        pt  = vor.point_region

        for i in 1:ndat
            region  =  reg[pt[i]+1]
            peri[i] =  voronoi_perimeter(ver,region)
            area[i] =  voronoi_area(ver,region)
        end
        println("ending voronoi ...")
        return(peri , area)
    end
end
## Gemini fix for voronoi....

function voronoi(pts, verbose = false)
    ndat = length(pts)
    
    if verbose
        println("Processing $ndat points...")
    end

    xx = Float64[p[1] for p in pts]
    yy = Float64[p[2] for p in pts]

    min_x, max_x = minimum(xx), maximum(xx)
    min_y, max_y = minimum(yy), maximum(yy)
    
    range_x = max_x - min_x
    range_y = max_y - min_y
    
    range_x = range_x == 0.0 ? 1.0 : range_x
    range_y = range_y == 0.0 ? 1.0 : range_y

    scale_x = 0.8 / range_x
    scale_y = 0.8 / range_y
    
    v_pts = VoronoiCells.IndexablePoint2D[]
    for i in 1:ndat
        # Keep the microscopic jitter just to be safe against exact duplicates
        jitter_x = 1e-6 * cos(Float64(i))
        jitter_y = 1e-6 * sin(Float64(i))
        
        # Data maps safely between 1.1 and 1.9
        nx = 1.1 + (xx[i] - min_x) * scale_x + jitter_x
        ny = 1.1 + (yy[i] - min_y) * scale_y + jitter_y
        push!(v_pts, VoronoiCells.IndexablePoint2D(nx, ny, i))
    end

    # --- THE FIX: ADD 4 GHOST CORNER POINTS ---
    # These sit outside the data (at 1.01 and 1.99) and act as a forced bounding box.
    # This guarantees NO stars (like point #21) have infinite boundaries.
    push!(v_pts, VoronoiCells.IndexablePoint2D(1.01, 1.01, ndat + 1))
    push!(v_pts, VoronoiCells.IndexablePoint2D(1.01, 1.99, ndat + 2))
    push!(v_pts, VoronoiCells.IndexablePoint2D(1.99, 1.01, ndat + 3))
    push!(v_pts, VoronoiCells.IndexablePoint2D(1.99, 1.99, ndat + 4))

    # Run the tessellation 
    tess = voronoicells(v_pts)
    
    # Compute scaled areas (This will no longer crash!)
    all_scaled_areas = voronoiarea(tess)

    # Throw away the 4 ghost points, keep only your real data (1 to ndat)
    real_scaled_areas = all_scaled_areas[1:ndat]

    # Reverse the mathematical scaling so DBSCAN gets accurate physical weights
    actual_areas = real_scaled_areas ./ (scale_x * scale_y)

    if verbose
        println("End voronoi...")
    end
    
    return (actual_areas, actual_areas)
end