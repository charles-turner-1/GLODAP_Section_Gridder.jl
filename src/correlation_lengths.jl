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

function calcCorrLengths(variable::Vector{Float64}
                                ;obsLat::Vector{Float64}
                                ,obsLon::Vector{Float64}
                                ,obsPres::Vector{Float64}
                                ,presGrid::Vector{Float64}
                                ,pressureStepNumber::Integer=10
                                ,verticalSearchRange::Number=100
                                ,lenxFactor::Number=1)
    # Calculate vertical and horizontal correlation lengths
    goodIdx = findNonNaNIndices(variable)

    length(goodIdx) < 1 ? error("No observations found") : nothing
    if obsLon == variable || obsLat == variable
        lenz = 10_000
    else
        lenz, infoz = fitvertlen((obsLon[goodIdx],obsLat[goodIdx],obsPres[goodIdx])
        ,variable[goodIdx],presGrid[1:pressureStepNumber:end],searchz=verticalSearchRange)
    end

    lenx = nothing
    for rangeFactor = 1:10
        printstyled("Trying to fit horizontal correlation length: attempt "* string(rangeFactor) * "\n",color=:yellow)
        try
            lenx, dbinfo = fithorzlen((obsLon[goodIdx],obsLat[goodIdx],obsPres[goodIdx])
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
    goodIdx = findNonNaNIndices(variable,obsSigma)
    lenz, infoz = fitvertlen((obsLon[goodIdx], obsLat[goodIdx], obsSigma[goodIdx])
    ,variable[goodIdx],sigGrid[1:sigmaStepNumber:end],searchz=verticalSearchRange
    ,smoothz=0.1)

    ####
    lenx = nothing
    for rangeFactor = 1:10
        printstyled("Trying to fit horizontal correlation length in density space: attempt "* string(rangeFactor) * "\n",color=:yellow)
        try
            lenx, dbinfo = fithorzlen((obsLon[goodIdx], obsLat[goodIdx], obsSigma[goodIdx])
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

