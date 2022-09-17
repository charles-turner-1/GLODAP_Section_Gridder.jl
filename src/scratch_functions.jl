function splitMeanAnomAlt(;obsVariable::Vector{Float64},obsPres::Vector{Float64} #kwargs might not be necessary since this is at the bottom of the call stack
                      ,obsLatLon::Vector{Float64},bgField::Matrix{Float64}
                      ,vertGrid::Vector{Float64},horzGrid::Vector{Float64})

    obsXvals = unique(obsLatLon)
    varAnom = fill(NaN,length(obsVariable))

    for xVal in obsXvals
        # Run through each station and find the closest horzGrid value. Then run
        # through the station and find the closest background field cell to each
        # observation and subtract the background field value there.

        colIdx = argmin(abs.(horzGrid .- xVal))

        stnPres = obsPres[obsLatLon.==xVal]
        stnVar = obsVariable[obsLatLon.==xVal]

        @threads for pVal in obsPres
            rowIdx = argmin(abs.(vertGrid .- pVal))
            bgVal = bgField[rowIdx,colIdx]

            # We now need the index of the observation we are looking at.
            obsIdxP = findall(obsPres .== pVal)
            obsIdxX = findall(obsLatLon .== xVal)

            obsIdx = intersect(obsIdxP,obsIdxX)

            length(obsIdx) > 0 ? varAnom[obsIdx[1]] = obsVariable[obsIdx[1]] - bgVal : nothing
        end
    end
    return varAnom
end

function findNaNIndices(variable::Vector)
    # Very simple function to find NaN indices of a vector
    goodIdx = findall(convert(Vector{Bool},fill(0,size(variable)) + isnan.(variable)))
    return goodIdx
end

function checkPartialCruise(;glodapDir::String="/home/ct/MATLAB/GLODAP"
                            ,goshipDir::String="/home/ct/MATLAB/GO_SHIP"
                            ,maskMatFile::String="/home/ct/Julia/GLODAP_Easy_Ocean/GOSHIP_MaskStruct.mat"
                            ,sectionName::String
                            ,horzCoordinate::String
                            ,expocode::AbstractString)

    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    # Since we just want to check if we have merged cruises, we can pass longitude
    # (or latitude) as our variable of interest.
    isAnException = testExpocodeException(expocode=expocode,variableName="G2longitude",maskMatfile=maskMatFile)

    gridDir, repDir, expocodeDir = load_GOSHIP_Directories(goshipDir)
    horzGrid, prGrid, sectionMask = loadSectionInfo(sectionName,MASK_MATFILE,gridDir)
    griddingVars = loadGLODAPVariables(glodapDir,["G2longitude","G2latitude"],GLODAP_expocode=expocode)

    if !isAnException
        griddingVars = loadGLODAPVariables(glodapDir,["G2longitude","G2latitude"],GLODAP_expocode=expocode)
    else
        (_, griddingVars) = loadExceptionData(expocode=expocode
        ,variableName="G2longitude",maskMatfile=maskMatFile)
    end

    lon = griddingVars["G2longitude"]
    lat = griddingVars["G2latitude"]

    if horzCoordinate == "longitude"
         obsXval = lon
         lon = matchLonConventions(horzGrid,lon)
    else
          obsXval = lat
    end

    maxGridVal = maximum(horzGrid)
    minGridVal = minimum(horzGrid)

    maxObsVal = maximum(obsXval)
    minObsVal = minimum(obsXval)

    ΔminVal = sign(minGridVal)*(minGridVal - minObsVal) # If minGrid << minObs, this will be +ve
    println("ΔminVal = " * string(ΔminVal))
    ΔmaxVal = sign(minGridVal)*(maxGridVal - maxObsVal) # If maxGrid << maxObs, this will be +ve
    println("ΔmaxVal = " * string(ΔmaxVal))

    isPartialCruise = false
    if ΔmaxVal > 2 || ΔminVal > 2
        isPartialCruise = true
    end

    return isPartialCruise
end
