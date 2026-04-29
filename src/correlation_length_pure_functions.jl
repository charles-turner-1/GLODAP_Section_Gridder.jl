"""
    corrlen_scale_factor(; horzCoordinate, latitudes=nothing, meanLatitude=nothing)

Return the legacy horizontal correlation-length scale factor in km per degree.
If `latitudes` is provided it intentionally overwrites `meanLatitude`, matching the
existing behaviour in `horzCorrDistanceKilometres`.
"""
function corrlen_scale_factor(; 
    horzCoordinate::Union{String,Nothing}=nothing,
    latitudes::Union{AbstractVector{<:Real},Nothing}=nothing,
    meanLatitude::Union{Real,Nothing}=nothing,
)
    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    if latitudes !== nothing
        meanLatitude = mean(latitudes)
    end

    if horzCoordinate == "longitude"
        return 111.2 .* cos.(meanLatitude * pi / 180)
    elseif horzCoordinate == "latitude"
        return 111.2
    end
end

"""
    corrlen_mean_latitude(obsLat)

Return the mean latitude after applying the legacy non-NaN filter.
"""
function corrlen_mean_latitude(obsLat::AbstractVector{<:Real})
    return mean(filter(!isnan, obsLat))
end

"""
    corrlen_finalize_lenz(lenz, lenx)

Apply the legacy sentinel expansion and 1000-unit cap to vertical correlation lengths.
"""
function corrlen_finalize_lenz(
    lenz::Union{Real,AbstractVector{<:Real}},
    lenx::AbstractVector{<:Real},
)
    if lenz == 10_000
        lenz = fill(lenz, size(lenx))
    end

    return min.(lenz, 1000)
end

"""
    corrlen_interpolate_to_grid(sampleGrid, sampledValues, fullGrid; extrapolation_bc)

Interpolate sampled correlation lengths back onto the full grid.
"""
function corrlen_interpolate_to_grid(
    sampleGrid::AbstractVector{<:Real},
    sampledValues::AbstractVector{<:Real},
    fullGrid::AbstractVector{<:Real};
    extrapolation_bc,
)
    interpolant = LinearInterpolation(sampleGrid, sampledValues, extrapolation_bc=extrapolation_bc)
    return interpolant.(fullGrid)
end

"""
    corrlen_retry_search(searchz, multiplier)

Return the legacy retry-wrapped search function used when fitting horizontal
correlation lengths with progressively larger search windows.
"""
function corrlen_retry_search(searchz::Function, multiplier::Integer)
    function search_z_func(z)
        return multiplier * searchz(z)
    end

    return search_z_func
end
