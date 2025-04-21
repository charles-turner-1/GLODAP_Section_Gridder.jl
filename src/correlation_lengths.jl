function horzCorrDistanceKilometres(horzLengthDegrees::Vector{Float64}
                               ;latitudes::Union{Vector{Float64},Nothing}=nothing
                               ,meanLatitude::Union{Float64,Nothing}=nothing
                               ,horzCoordinate::Union{String,Nothing}=nothing)
    # Calculate horizontal correlation length in Kilometres
    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    if latitudes !== nothing
        meanLatitude = mean(latitudes)
        # This will overwrite meanLatitude if both are given which I think is probably
        # the optimal behaviour
    end

    if horzCoordinate == "longitude"
        scaleFactor = 111.2 .* cos.(meanLatitude*pi/180)
    elseif horzCoordinate == "latitude"
        scaleFactor = 111.2
    end

    return scaleFactor * horzLengthDegrees
end

function calcCorrLengths(;variable::Vector{Float64},
                        obsLat::Vector{Float64},
                        obsLon::Vector{Float64},
                        obsPres::Vector{Float64},
                        presGrid::Vector{Float64},
                        pressureStepNumber::Integer=10,
                        verticalSearchRange::Number=100,
                        lenxFactor::Number=1,
    )
    # Calculate vertical and horizontal correlation lengths
    goodIdx = non_nan_indices(variable)

    length(goodIdx) < 1 ? error("No observations found") : nothing
    if obsLon == variable || obsLat == variable
        lenz = 10_000
    else
        lenz, _ = fitvertlen((obsLon[goodIdx],obsLat[goodIdx],obsPres[goodIdx])
        ,variable[goodIdx],presGrid[1:pressureStepNumber:end],searchz=verticalSearchRange)
    end

    lenx = nothing
    for rangeFactor = 1:10
        printstyled("Trying to fit horizontal correlation length: attempt "* string(rangeFactor) * "\n",color=:yellow)
        try
            lenx, _ = fithorzlen((obsLon[goodIdx],obsLat[goodIdx],obsPres[goodIdx])
            ,variable[goodIdx],presGrid[1:pressureStepNumber:end],searchz=rangeFactor^2*verticalSearchRange)
        catch
        end
        lenx !== nothing ? break : nothing
    end

    lenz == 10_000 ? lenz = fill(lenz, size(lenx)) : nothing


    lenz = min.(lenz,1000)

    meanLat = mean(filter(!isnan,obsLat))

    lenxKM = horzCorrDistanceKilometres(lenx,meanLatitude=meanLat,horzCoordinate="longitude")

    if pressureStepNumber > 1
        x_interpolant = LinearInterpolation(presGrid[1:pressureStepNumber:end],lenxKM,extrapolation_bc=NaN)
        z_interpolant = LinearInterpolation(presGrid[1:pressureStepNumber:end],lenz,extrapolation_bc=NaN)
    end

    lenz= z_interpolant.(presGrid)
    lenxKM = x_interpolant.(presGrid)

    lenzSmthd = DIVAnd.smoothfilter(presGrid,lenz,400) # Should presGrid have indices ie. [1:stepno:end] ?
    lenxKM = lenxFactor * DIVAnd.smoothfilter(presGrid,lenxKM,400)

    return lenzSmthd, lenxKM
end

function corrLengthVectorToMatrix(lenX::Vector{Float64},lenZ::Vector{Float64}
                                 ,presGrid::Vector{Float64}, llGrid::Vector{Float64})
    # Take the correlation length vectors, and turn them into matrices the same 
    # size as our mask
    sz = (length(presGrid),length(llGrid))
    lenXgrid = fill(1.0,sz)
    lenZgrid = fill(1.0,sz)
    for i = 1:length(llGrid)
        lenXgrid[:,i] = lenX
        lenZgrid[:,i] = lenZ
    end
    lenXgrid = convert(Matrix{Float64},lenXgrid)
    lenZgrid = convert(Matrix{Float64},lenZgrid)

    lenZgrid = min.(lenZgrid,1000)

    return lenZgrid, lenXgrid
end

function calcDensityCorrLengths(variable::Vector{Float64}
                                ;obsLat::Vector{Float64}
                                ,obsLon::Vector{Float64}
                                ,obsSigma::Vector{Float64}
                                ,sigGrid::Vector{Float64}
                                ,lenxPrescribed::Float64
                                ,sigmaStepNumber::Integer=10
                                ,verticalSearchRange::Float64=0.0001)
    # Calculate correlation lengths in density space.
    goodIdx = non_nan_indices(variable,obsSigma)
    lenz, _ = fitvertlen((obsLon[goodIdx], obsLat[goodIdx], obsSigma[goodIdx])
    ,variable[goodIdx],sigGrid[1:sigmaStepNumber:end],searchz=verticalSearchRange
    ,smoothz=0.1)

    ####
    lenx = nothing
    for rangeFactor = 1:10
        printstyled("Trying to fit horizontal correlation length in density space: attempt "* string(rangeFactor) * "\n",color=:yellow)
        try
            lenx, _ = fithorzlen((obsLon[goodIdx], obsLat[goodIdx], obsSigma[goodIdx])
            ,variable[goodIdx],sigGrid[1:sigmaStepNumber:end],searchz=0.1*rangeFactor
            ,smoothz=0.1)
        catch
        end
        if lenx !== nothing
            break
        end
    end

    if lenx === nothing
        lenx = fill(lenxPrescribed,size(lenz))
        # Give it a prescribed horizontal correlation length
    end

    meanLat = mean(filter(!isnan,obsLat))

    lenxKM = horzCorrDistanceKilometres(lenx,meanLatitude=meanLat,horzCoordinate="longitude")

    if sigmaStepNumber > 1
        x_interpolant = LinearInterpolation(sigGrid[1:sigmaStepNumber:end],lenxKM,extrapolation_bc=Flat())
        z_interpolant = LinearInterpolation(sigGrid[1:sigmaStepNumber:end],lenz,extrapolation_bc=Flat())
    end

    lenz= z_interpolant.(sigGrid)
    lenxKM = x_interpolant.(sigGrid)
    return lenz, lenxKM
end


"""
    fit_lengths(
        vars::DataFrame,
        data_residual::Vector{Float64},
        pr_grid::Vector{Real}, 
        search_z_func::Function
    )::Tuple{Vector{Float64}, Vector{Float64}}

Fits horizontal and vertical correlation lengths for a given dataset using residuals and pressure grid.

# Arguments
- `vars::DataFrame`: A DataFrame containing longitude, latitude, and pressure.
- `data_residual::Vector{Float64}`: The residuals of the data to be used for fitting.
- `pr_grid::Vector{Real}`: The pressure grid over which the correlation lengths are computed.
- `search_z_func::Function`: A function used to search for the vertical correlation length.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: A tuple containing:
  - `lenx::Vector{Float64}`: The horizontal correlation lengths.
  - `lenz::Vector{Float64}`: The vertical correlation lengths.

# Workflow
1. Computes horizontal correlation lengths using the `fithorzlen` function from DIVAnd.
2. Computes vertical correlation lengths using the `fitvertlen` function from DIVAnd.

# Notes
- The horizontal correlation lengths are computed based on longitude, latitude, and pressure.
- The vertical correlation lengths are constrained to a specific range using a limiting function.
- The `search_z_func` parameter customizes the search behavior for vertical correlation lengths.
- It is assumed that only residual data, not the original data, is passed to this function. 
    If original data is used, correlation length fitting may fail.
- Applies a limiting function to constrain vertical correlation lengths between 10 and 1000.
"""
function fit_lengths(
    vars::DataFrame,
    data_residual::Vector{Float64},
    pr_grid::Vector{Real}, 
    search_z_func::Function
)::Tuple{Vector{Float64}, Vector{Float64}}
    @info "Computing correlation lengths: horizontal"
    _x = (vars[!, "G2longitude"], vars[!, "G2latitude"], vars[!, "G2pressure"])
    lenx, _dbinfo = fithorzlen(_x, data_residual, pr_grid, searchz=search_z_func)

    @info "Computing correlation lengths: vertical"
    lenz, _dbinfo = fitvertlen(
        _x, 
        data_residual, 
        pr_grid,
        searchz=search_z_func, 
        limitfun= (z, len) -> max(min(len, 1000), 10)
)

    return lenx, lenz
end
