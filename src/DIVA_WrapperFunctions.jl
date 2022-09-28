#= The aim of this file is to create a bunch of functions which we are going to
use in order to generalise `Diva_RunScript_Simple.jl` to a more functional,
general purpose tool for interpolating any GLODAP data onto a regular grid, as
defined in the GO-SHIP Easy Ocean toolbox. As such, this file itself is going to
be just a bunch of functions which can be called when you know what variables
you want to interpolate from what sections.    =#


function easyDIVAGrid(;variable::Vector{Float64}
                      ,vertVar::Vector{Float64}
                      ,latLon::Vector{Float64}
                      ,horzCoordinate::String
                      ,meanValue::String="scalar"
                      ,bgField::Union{Matrix{Float64},Nothing}=nothing
                      ,vertGrid::Vector{Float64}
                      ,horzGrid::Vector{Float64}
                      ,horzScale::Matrix{Float64}
                      ,vertScale::Matrix{Float64}
                      ,mask::Matrix{Bool}
                      ,Epsilon::Float64
                      ,horzCorrLength::Vector{Float64}
                      ,vertCorrLength::Vector{Float64})
    # This is the wrapped that calls DIVAnd.

    if horzCoordinate == "longitude"
        latLon = matchLonConventions(horzGrid,latLon)
    elseif horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    if meanValue == "horzMean"
        varMean, varAnom = splitMeanAnom(variable,vertVar,vertGrid,horzGrid)
    elseif meanValue == "scalar"
        varMean, varAnom = splitMeanAnom(variable)
    elseif meanValue == "climatology"
        if bgField === nothing
            error("If meanValue is specified to be \"climatology\" a background
            field must be provided")
        else
            varAnom = splitMeanAnom(;obsVariable=variable,obsPres=vertVar
                                    ,obsLatLon=latLon,bgField=bgField
                                    ,vertGrid=vertGrid,horzGrid=horzGrid)
            varMean = bgField
        end
    end

    # I think it's necessary to filter out NaN indices or we wind up with shitty
    # fields coming out
    goodIdx = findNonNaNIndices(variable)

    P_grid,L_grid = ndgrid(vertGrid,horzGrid)
    corrLengths = corrLengthVectorToMatrix(horzCorrLength,vertCorrLength,vertGrid,horzGrid)

    fi,s = DIVAndrun(mask, (vertScale,horzScale), (P_grid, L_grid), (vertVar[goodIdx],latLon[goodIdx]), varAnom[goodIdx], corrLengths, Epsilon)

    meanValue == "scalar" ? griddedVariable = fi .+ varMean : griddedVariable = fi + varMean
    return griddedVariable
end

function createDensitySpaceMask(obsSigma::Vector{Float64}, obsLonLat::Vector{Float64}
                               ,lonLatGrid::Vector{Float64}, sigmaGrid::Vector{Float64}
                               ,pressureSpaceMask::Matrix{Bool})
    # Create a mask for gridding observations in density space. This will stop 
    # flat grids gridding through ridges and whatnot. Might not actually be necessary.
    uniqueStations = unique(obsLonLat)
    sort!(uniqueStations)

    stationMaxSigma = fill(NaN,size(uniqueStations))
    stationMinSigma = fill(NaN,size(uniqueStations))

    stationMinSigmaIdx = fill(NaN,size(uniqueStations))
    stationMaxSigmaIdx = fill(NaN,size(uniqueStations))

    for i in eachindex(uniqueStations)
        stationSigma = obsSigma[findall(obsLonLat .== uniqueStations[i])]
        stationSigma = filter(!isnan,stationSigma)
        try # Try/catch is to avoid max/min coming out as nan, causing issues
            stationMaxSigma[i] = stationSigma[findmax(stationSigma,climTemp_col)[2]]
            stationMinSigma[i] = stationSigma[findmin(stationSigma)[2]]
            stationMinSigmaIdx[i] = findmin(abs.(fill(stationMinSigma[i],size(sigmaGrid)) - sigmaGrid))[2]
            stationMaxSigmaIdx[i] = findmin(abs.(fill(stationMaxSigma[i],size(sigmaGrid)) - sigmaGrid))[2]
        catch
            stationMaxSigma[i] = NaN
            stationMinSigma[i] = NaN
            stationMinSigmaIdx[i] = NaN
            stationMaxSigmaIdx[i] = NaN
        end
    end

    sigmaMask = fill(1,length(sigmaGrid),length(lonLatGrid))
    badCols = sum(pressureSpaceMask,dims=1) ./ sum(pressureSpaceMask,dims=1)
    replace!(badCols,NaN => 0)
    badCols = vec(convert(Matrix{Bool},badCols))


    siMax = LinearInterpolation(uniqueStations,stationMaxSigmaIdx,extrapolation_bc=NaN)
    siMin = LinearInterpolation(uniqueStations,stationMinSigmaIdx,extrapolation_bc=NaN)

    sigMinGrid = round.(siMin.(lonLatGrid))
    sigMaxGrid = round.(siMax.(lonLatGrid))

    for i = 1:length(lonLatGrid)
        sigColumn = sigmaMask[:,i]
        try
            minVal = convert(Int64,sigMinGrid[i])
            maxVal = convert(Int64,sigMaxGrid[i])
            sigColumn[1:minVal-1] .= 0
            sigColumn[maxVal+1:end] .= 0
        catch
            sigColumn .=  0
        end
        sigmaMask[:,i] = sigColumn
    end

    rmIndices = findall(badCols .==0)
    sigmaMask[:,rmIndices] .=0
    sigmaMask = convert(Matrix{Bool},sigmaMask)
    return sigmaMask
end

function easyDIVAisopycnal(;obsVariable::Vector{Float64}
                      ,obsSigma::Vector{Float64}
                      ,obsPressure::Vector{Float64}
                      ,obsLat::Vector{Float64}
                      ,obsLon::Vector{Float64}
                      ,latLon::Vector{Float64}
                      ,vertGrid::Vector{Float64}
                      ,horzGrid::Vector{Float64}
                      ,horzScale::Matrix{Float64}
                      ,vertScale::Matrix{Float64}
                      ,pressureMask::Matrix{Bool}
                      ,Epsilon::Float64
                      ,horzCorrLength::Vector{Float64}
                      ,horzCorrLengthPrescribed::Float64
                      ,vertCorrLength::Vector{Float64})
    # Applies DIVA gridding on density surfaces.
    if latLon == obsLon
        obsLon = matchLonConventions(horzGrid,obsLon)
        horzCoordinate = "longitude"
    else
        horzCoordinate = "latitude"
    end
    horzDist = gridHorzDistance(obsLat,obsLon,horzGrid)

    if horzCorrLengthPrescribed != 0
        horzCorrLength = fill(horzCorrLengthPrescribed,size(horzCorrLength))
    end

    sigGrid = createSigmaGrid(obsSigma)
    sigDist = gridSigDistance(sigGrid)
    scaleDens, scaleHorz = calcScaleFactors(sigDist,horzDist)

    sigmaMask = createDensitySpaceMask(obsSigma,latLon,horzGrid,sigGrid,pressureMask)
    sigmaMask .= true
    # The above fix is a complete hack and needs sorting at some point


    lenSig = calcDensityCorrLengths(obsVariable,obsLat=obsLat,obsLon=obsLon
    ,obsSigma=obsSigma,sigGrid=sigGrid,sigmaStepNumber=10,verticalSearchRange=0.01
    ,lenxPrescribed=horzCorrLengthPrescribed)


    S_grid,L_grid = ndgrid(sigGrid,horzGrid)

    griddedVariableEasySigma = easyDIVAGrid(variable=obsVariable,vertVar=obsSigma
    ,horzCoordinate=horzCoordinate,latLon=latLon,vertGrid=sigGrid,horzGrid=horzGrid
    ,horzScale=scaleHorz,vertScale=scaleDens,mask=sigmaMask,Epsilon=0.01
    ,horzCorrLength=lenSig[2],vertCorrLength=lenSig[1])

    griddedPressureEasySigma = easyDIVAGrid(variable=obsPressure,vertVar=obsSigma
    ,horzCoordinate=horzCoordinate,latLon=latLon,vertGrid=sigGrid,horzGrid=horzGrid
    ,horzScale=scaleHorz,vertScale=scaleDens,mask=sigmaMask,Epsilon=0.01
    ,horzCorrLength=5*lenSig[2],vertCorrLength=lenSig[1])

    griddedVarSigmaVector = griddedVariableEasySigma[:]
    griddedPresSigmaVector = griddedPressureEasySigma[:]
    S_gridVector = S_grid[:]
    L_gridVector = L_grid[:]

    horzDist = gridHorzDistance(obsLat,obsLon,horzGrid)
    vertDist = gridVertDistance(vertGrid)
    scaleVert, scaleHorz = calcScaleFactors(vertDist,horzDist)
    # We are getting some weird issues when reinterpolating back to pressure space
    # from density space. Presumably related to gridding?
    #sigmaGriddedVariable = easyDIVAGrid(variable=griddedVarSigmaVector,vertVar=griddedPresSigmaVector,latLon=L_gridVector,vertGrid=prGrid,horzGrid=lonGrid,horzScale=scaleHorz,vertScale=scaleVert,mask=pressureMask,Epsilon=0.01,horzCorrLength=len[2],vertCorrLength=len[1])

    sigmaGriddedVariable = easyDIVAGrid(variable=griddedVarSigmaVector,horzCoordinate=horzCoordinate,
    vertVar=griddedPresSigmaVector,latLon=L_gridVector,vertGrid=vertGrid,horzGrid=horzGrid,
    horzScale=scaleHorz,vertScale=scaleVert,mask=pressureMask,Epsilon=0.01,
    horzCorrLength=horzCorrLength,vertCorrLength=vertCorrLength)

    return sigmaGriddedVariable
end



function easyDIVACrossValidate(;variable::Vector{Float64}
                               ,vertVar::Vector{Float64}
                               ,latLon::Vector{Float64}
                               ,horzCoordinate::String
                               ,vertGrid::Vector{Float64}
                               ,horzGrid::Vector{Float64}
                               ,horzScale::Matrix{Float64}
                               ,vertScale::Matrix{Float64}
                               ,mask::Matrix{Bool}
                               ,Epsilon::Float64
                               ,horzCorrLength::Vector{Float64}
                               ,vertCorrLength::Vector{Float64}
                               ,numLengthVals::Int64=5
                               ,numEpsilonVals::Int64=5
                               ,methodVal::Int64=0)
    # Perform DIVA cross validation
    if horzCoordinate == "longitude"
        latLon = matchLonConventions(horzGrid,latLon)
    elseif horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\"
        or \"latitude\"")
    end

    P_grid,L_grid = ndgrid(vertGrid,horzGrid)
    varMean, varAnom = splitMeanAnom(variable)
    corrLengths = corrLengthVectorToMatrix(horzCorrLength,vertCorrLength,vertGrid,horzGrid)

    bestFactorLen,bestFactorε, cvval,cvvalues, x2Ddata,y2Ddata,cvinter,xi2D,yi2D =
    DIVAnd_cv(mask, (vertScale,horzScale), (P_grid, L_grid), (vertVar,latLon),
     varAnom, corrLengths, Epsilon, numLengthVals, numEpsilonVals, methodVal)

    return bestFactorLen, bestFactorε
end
