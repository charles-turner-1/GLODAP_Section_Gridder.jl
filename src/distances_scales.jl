# function gridHorzDistance(GLODAP_latitudes::Vector{Float64} ,GLODAP_longitudes::Vector{Float64} ,latlonGrid)
function grid_horz_dist(
    GLODAP_latitudes::AbstractVector{<:Real},
    GLODAP_longitudes::AbstractVector{<:Real}, 
    ll_grid::AbstractVector{<:Real},
)::AbstractVector{<:AbstractFloat}
    # Compute the mean distance between each station in a cruise

    # This needs to be completely redone. What we want is to somehow identify 
    # each unique station 

    uniq_locs = unique(zip(GLODAP_latitudes, GLODAP_longitudes))


    lon_range = maximum(GLODAP_longitudes) - minimum(GLODAP_longitudes)
    lat_range = maximum(GLODAP_latitudes) - minimum(GLODAP_latitudes)
    lon_range < lat_range ?   # Probably redundant, look into removing
        sort!(uniq_locs, by=x -> x[1]) : 
        sort!(uniq_locs, by=x -> x[2])

    GLODAP_longitudes = [uniq_locs[i][1] for i = 1:length(uniq_locs)]
    GLODAP_latitudes = [uniq_locs[i][2] for i = 1:length(uniq_locs)]

    d_lon = central_diff(GLODAP_longitudes)
    d_lat = central_diff(GLODAP_latitudes)

    d_lat_m = d_lat * 111.2
    d_lon_m = d_lon * 111.2 .* cos.(GLODAP_latitudes * pi / 180)

    horzdist_km = sum(sqrt.(d_lat_m .^ 2 + d_lon_m .^ 2)) / length(ll_grid)
    horzdist_km = fill(horzdist_km, size(ll_grid))

    return horzdist_km
end

# function gridVertDistance(pressureGrid ,GLODAP_pressures=nothing)
"""
    grid_vert_dist(pgrid::Vector{<:Real}, GLODAP_pressures=nothing) -> Vector{<:Integer}
Calculate the vertical distance for a given pressure grid. GO-SHIP Easy Ocean grids 
everything onto a 10m vertical grid, so this function just returns a vector of 10s.

Provided for consistency & convenience
"""
function grid_vert_dist(pgrid)::Vector{<:Integer}
    return fill(10, size(pgrid))
end


# function createSigmaGrid(sigmaVals::Vector{Float64},numLevels::Int64=600)
function create_sigma_grid(sigma_vals::AbstractVector{<:AbstractFloat}, numLevels::Integer=600)
    # Creates a 600 level (default) evenly spaced grid in density space
    sigma_vals = unique(sort(filter(!isnan, sigma_vals)))
    step = convert(Int64, ceil(length(sigma_vals) / numLevels))
    sigma_grid = sigma_vals[1:step end]
    return sigma_grid
end

# function gridSigDistance(sigmaGrid::Vector{Float64})
function grid_sigma_distance(
    sigma_grid::AbstractVector{<:AbstractFloat}
)::Vector{<:AbstractFloat}
    minval, maxval = extrema(sigma_grid)
    mean_dist = (maxval - minval) / length(sigma_grid)
    sigma_mean_dist = fill(mean_dist, size(sigma_grid))

    return sigma_mean_dist
end

# function calcScaleFactors(verticalDistance::Vector,horizontalDistance::Vector ;printScales=true)
"""
    calc_scale_factors(vert_dist::Vector{<:Real}, horz_dist::Vector{<:Real}; print_scales=true) -> ScaleFactors

Calculate vertical and horizontal scale factors based on the provided vertical and horizontal distances.

# Arguments
- `vert_dist::Vector{<:Real}`: A vector of vertical distances.
- `horz_dist::Vector{<:Real}`: A vector of horizontal distances.
- `print_scales::Bool` (optional): If `true`, prints the mean horizontal and vertical scale factors. Defaults to `true`.

# Returns
- `ScaleFactors`: A custom type containing the calculated vertical and horizontal scale factors as arrays.

# Notes
- The function computes the scale factors as the reciprocal of the distances.
- The `ndgrid` function is used to create grids of the input distances for computation.
"""
function calc_scale_factors(
    vert_dist::Vector{<:Real},
    horz_dist::Vector{<:Real};
    print_scales=true
)::ScaleFactors
    # Calculate vertical and horizontal scale factors, print them out (if desired)
    # and return them
    dP_grid, dL_grid = ndgrid(vert_dist, horz_dist)
    scaleVert = ones(size(dP_grid)) ./ dP_grid
    scaleHorz = ones(size(dL_grid)) ./ dL_grid
    if print_scales
        println("Horizontal scale factor: ", mean(scaleHorz))
        println("Vertical scale factor: ", mean(scaleVert))
    end
    return ScaleFactors(scaleVert, scaleHorz)
end

"""
Contains the scale factors for the vertical and horizontal distances
"""
struct ScaleFactors
    vert::Vector{<:Real}
    horz::Vector{<:Real}
end
