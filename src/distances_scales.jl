# function gridHorzDistance(GLODAP_latitudes::Vector{Float64} ,GLODAP_longitudes::Vector{Float64} ,latlonGrid)
function grid_horz_dist(GLODAP_latitudes::AbstractVector{<:Real}, GLODAP_longitudes::AbstractVector{<:Real}, latlonGrid)
    # Compute the mean distance between each station in a cruise

    # This needs to be completely redone. What we want is to somehow identify 
    # each unique station 

    uniqueLocations = unique(zip(GLODAP_latitudes, GLODAP_longitudes))
    lonRange = maximum(GLODAP_longitudes) - minimum(GLODAP_longitudes)
    latRange = maximum(GLODAP_latitudes) - minimum(GLODAP_latitudes)
    lonRange < latRange ? sort!(uniqueLocations, by=x -> x[1]) :
    sort!(uniqueLocations, by=x -> x[2])
    GLODAP_longitudes = [uniqueLocations[i][1] for i = 1:length(uniqueLocations)]
    GLODAP_latitudes = [uniqueLocations[i][2] for i = 1:length(uniqueLocations)]

    dLon = central_diff(GLODAP_longitudes)
    dLat = central_diff(GLODAP_latitudes)

    dLat_m = dLat * 111.2
    dLon_m = dLon * 111.2 .* cos.(GLODAP_latitudes * pi / 180)

    horzDist_kilometres = sum(sqrt.(dLat_m .^ 2 + dLon_m .^ 2)) / length(latlonGrid)
    horzDist_kilometres = fill(horzDist_kilometres, size(latlonGrid))

    return horzDist_kilometres
end

# function gridVertDistance(pressureGrid ,GLODAP_pressures=nothing)
function grid_vert_dist(pgrid)::Vector{Int32}
    # Because GO-SHIP Easy Ocean grids everything onto a 10m vertical grid, all
    # we need to do is return a distance of 10m
    vertDist_metres = fill(10, size(pgrid))
    return vertDist_metres
end

# For some reason which escapes me, I can't type either of the grid___Distance
# functions without breaking them. Something to understand

# function createSigmaGrid(sigmaVals::Vector{Float64},numLevels::Int64=600)
function create_sigma_grid(sigma_vals::AbstractVector{<:AbstractFloat}, numLevels::Int32=600)
    # Creates a 600 level (default) evenly spaced grid in density space
    sigma_vals = unique(sort(filter(!isnan, sigma_vals)))
    sigmaStep = convert(Int64, ceil(length(sigma_vals) / numLevels))
    sigmaGrid = sigma_vals[1:sigmaStep:end]
    return sigmaGrid
end

# function gridSigDistance(sigmaGrid::Vector{Float64})
function grid_sigma_distance(sigma_grid::AbstractVector{Float64})::Vector{Float64}
    # Work out the characteristic mean sigma grid distance. Could probably merge
    # this with the previous function
    sigmaMeanDist = mean(central_diff(sigma_grid))
    sigmaMeanDist = fill(sigmaMeanDist, size(sigma_grid))

    return sigmaMeanDist
end

# function calcScaleFactors(verticalDistance::Vector,horizontalDistance::Vector ;printScales=true)
function calc_scale_factors(vert_dist::Vector{<:Real}, horz_dist::Vector{<:Real} ; print_scales=true)::Tuple{Vector,Vector}
    # Calculate vertical and horizontal scale factors, print them out (if desired)
    # and return them
    dP_grid, dL_grid = ndgrid(vert_dist, horz_dist)
    scaleVert = ones(size(dP_grid)) ./ dP_grid
    scaleHorz = ones(size(dL_grid)) ./ dL_grid
    if print_scales
        println("Horizontal scale factor: ", mean(scaleHorz))
        println("Vertical scale factor: ", mean(scaleVert))
    end
    return scaleVert, scaleHorz
end