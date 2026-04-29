"""
    legacy_obs_xvals(horzCoordinate, horzGrid, obsLat, obsLon)

Return the horizontal observation coordinate used by the legacy camelCase pipeline.
This intentionally preserves the old branch structure, including the unused
longitude normalisation assignment, so behaviour stays unchanged while the logic
becomes easier to test.
"""
function legacy_obs_xvals(
    horzCoordinate::String,
    horzGrid::AbstractVector{<:Real},
    obsLat::AbstractVector{<:Real},
    obsLon::AbstractVector{<:Real},
)
    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    if horzCoordinate == "longitude"
        obsXval = obsLon
        obsLon = modulo_lon(horzGrid, obsLon)
    else
        obsXval = obsLat
    end

    return obsXval
end

"""
    legacy_partial_cruise_deltas(horzGrid, obsXval)

Return the legacy min/max section-minus-observation deltas used to detect partial cruises.
"""
function legacy_partial_cruise_deltas(
    horzGrid::AbstractVector{<:Real},
    obsXval::AbstractVector{<:Real},
)
    maxGridVal = maximum(horzGrid)
    minGridVal = minimum(horzGrid)

    maxObsVal = maximum(obsXval)
    minObsVal = minimum(obsXval)

    ΔminVal = minGridVal - minObsVal
    ΔmaxVal = maxGridVal - maxObsVal

    return (ΔminVal=ΔminVal, ΔmaxVal=ΔmaxVal, minObsVal=minObsVal, maxObsVal=maxObsVal)
end

"""
    legacy_truncate_mask(mask, horzGrid, minObsVal, maxObsVal)

Return a copy of the legacy section mask truncated to the observed horizontal span.
"""
function legacy_truncate_mask(
    mask::AbstractMatrix{Bool},
    horzGrid::AbstractVector{<:Real},
    minObsVal::Real,
    maxObsVal::Real,
)
    llGrid = Matrix{Float64}(undef, size(mask'))
    [llGrid[:, i] = horzGrid for i in 1:size(mask, 1)]
    llGrid = llGrid'

    truncatedMask = copy(mask)
    truncatedMask[llGrid .< minObsVal] .= false
    truncatedMask[llGrid .> maxObsVal] .= false

    return truncatedMask
end

"""
    legacy_variable_exception_good_indices(variable_length, station, pressure, stationList, minPressureList, maxPressureList, colIdx)

Apply the legacy manual exclusion rules and return the surviving observation indices.
"""
function legacy_variable_exception_good_indices(
    variable_length::Integer,
    station::AbstractVector{<:Real},
    pressure::AbstractVector{<:Real},
    stationList::AbstractVector{<:Integer},
    minPressureList::AbstractVector{<:Integer},
    maxPressureList::AbstractVector{<:Integer},
    colIdx::AbstractVector{<:Integer},
)
    goodIdx = collect(1:variable_length)

    for idx in colIdx
        criteriaMatched = fill(0, size(goodIdx))
        stnVal = stationList[idx]
        maxPrsVal = maxPressureList[idx]
        minPrsVal = minPressureList[idx]
        criteriaMatched[station .== stnVal] .+= 1
        criteriaMatched[pressure .< maxPrsVal] .+= 1
        criteriaMatched[pressure .> minPrsVal] .+= 1

        goodIdx[criteriaMatched .== 3] .= -1
    end

    return goodIdx[goodIdx .> 0]
end

"""
    legacy_horzlen_factor(exceptionDataFrame, expocode, variableName, griddingType)

Resolve the legacy manual horizontal correlation length override from a preloaded table.
"""
function legacy_horzlen_factor(
    exceptionDataFrame::DataFrame,
    expocode::AbstractString,
    variableName::AbstractString,
    griddingType::AbstractString,
)
    EDFsubset = exceptionDataFrame[
        (exceptionDataFrame.Expocode .== expocode) .&
        (exceptionDataFrame.Variable .== variableName) .&
        (exceptionDataFrame.Gridding .== griddingType),
        :,
    ]

    if prod(size(EDFsubset)) > 0
        println("Expocode contains horizonal correlation length exception for selected variable and gridding")
        factor = EDFsubset[!, "Factor"]
        return factor[1]
    elseif griddingType == "isobaric"
        return 1
    elseif griddingType == "isopycnic"
        return 1.0
    end
end
