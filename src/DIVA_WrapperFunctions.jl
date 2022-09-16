#= The aim of this file is to create a bunch of functions which we are going to
use in order to generalise `Diva_RunScript_Simple.jl` to a more functional,
general purpose tool for interpolating any GLODAP data onto a regular grid, as
defined in the GO-SHIP Easy Ocean toolbox. As such, this file itself is going to
be just a bunch of functions which can be called when you know what variables
you want to interpolate from what sections.    =#

# For now, lets just focus on obtaining a single occupation. Then we can add in
# some of that handy optional typing in order to obtain a series more easily

using MATLAB
using DIVAnd
using CSV
using DelimitedFiles
using NCDatasets
using DataFrames
using Statistics
using Interpolations
using FiniteDifferences
using Polynomials
using Base.Threads

function loadGoShipDirectories(GOSHIP_DIR)
    # Give the paths of the GO_SHIP data and return the files we
    # need to load

    GRID_DIR = joinpath(GOSHIP_DIR , "go_ship_clean_ctd/gridded/" )
    REP_DIR  = joinpath(GOSHIP_DIR , "go_ship_clean_ctd/reported/")
    CONV_DIR = joinpath(GOSHIP_DIR , "codeConversions/SectionFiles") # Should really move this to live with GLODAP_Easy_Ocean

    return GRID_DIR, REP_DIR, CONV_DIR
end

function listAvailableMasks(Mask_MatFile::String)
    # Reads the section mask file and returns the available masks.
    SectionMaskFile = MatFile(Mask_MatFile)
    maskDict = jdict(get_mvariable(SectionMaskFile,"maskStruct"))
    println("\nAvailable WOCE Section Masks:")
    for key in maskDict
        println(key.first)
    end
end

function loadSectionInfo(Mask_MatFile::String, sectionName::String,
    Grid_Directory::String = "go_ship_clean_ctd/gridded/")
    # Loads the grid data for the section, ie. lat/lon grid, pressure grid, mask.
    SectionMaskFile = MatFile(Mask_MatFile)
    maskDict = jdict(get_mvariable(SectionMaskFile,"maskStruct"))
    sectionMaskName = maskNameFromSectionName(sectionName)
    mask = maskDict[sectionMaskName]
    # Convert to the type we need
    mask = convert(Array{Bool},(mask .== 1))

    # Now we need to also return the grid.
    sectionName = replace(sectionName,"_" => "-")
    basin = sectionName[1]
    matName  = lowercase(sectionName) * ".mat"
    # Julia has no switch (WTF)
    if basin == 'A'
        maskPath = joinpath("atlantic/",sectionName,matName)
    elseif basin == 'P'
        maskPath = joinpath("pacific/",sectionName,matName)
    elseif basin == 'I'
        maskPath = joinpath("indian/",sectionName,matName)
    elseif basin == 'S'
        maskPath = joinpath("southern/",sectionName,matName)
    end

    gridFile = MatFile(joinpath(Grid_Directory,maskPath))
    pr_grid = vec(get_variable(gridFile,"pr_grid"))
    ll_grid = vec(get_variable(gridFile,"ll_grid"))

    if basin != 'P'
        ll_grid[ll_grid .> 180] .-=360 # So GLODAP and GO-SHIP data are using the
        # same longitudes. Probably will want to update this to be basin dependent
        # or we might run into problems in the Pacific
    end

    return ll_grid, pr_grid, mask
end

function matchLonConventions(lonGrid::Vector{Float64}, lonVals::Vector{Float64})
    # GO-SHIP Easy Ocean toolbox and GLODAP use different longitude conventions
    # ie. lon ⋹ [-180,180] or lon ⋹ [0,360]. This will ensure tat the match up.
    if maximum(lonGrid) > 180
        lonVals[lonVals .< 0] .+= 360
    end
    return lonVals
end

function maskNameFromSectionName(sectionName::String)
    # Gets the mask name from the name of the section we're using.
    maskName = "mask" * sectionName
    maskName = replace(maskName,"-" => "_")
    return maskName
end

function listSectionExpocodes(sectionName::String,convDir::String)
    # Lists out all the expocodes of cruises occupying a given section
    expocodes = joinpath(convDir,sectionName) * ".csv"
    expocodes = CSV.read(expocodes,DataFrame)
    return expocodes
end


function listAvailableGLODAPVariables(GLODAP_DIR::String,
    GLODAP_DATAFILE::String="GLODAPv2.2021_Merged_Master_File.mat")
    # Lists all variables contained in GLODAP

    GLODAP_DATAFILE = joinpath(GLODAP_DIR,GLODAP_DATAFILE)
    GLODAP_Data = MatFile(GLODAP_DATAFILE)
    println(" ");println("Available GLODAP variables:")
    for variable in variable_names(GLODAP_Data)
        println(variable)
    end
end

function loadGLODAPVariable(GLODAP_DIR::String,GLODAP_VariableName::String
    ;GLODAP_DATA_FILENAME::String="GLODAPv2.2021_Merged_Master_File.mat"
    ,GLODAP_expocode::Union{String,String15,Nothing} = nothing)
    # Loads a single variable from GLODAP, for either a given expocode or the 
    # entirety of the GLODAP dataset.

    GLODAP_DATAFILE = joinpath(GLODAP_DIR,GLODAP_DATA_FILENAME)
    GLODAP_Data = MatFile(GLODAP_DATAFILE)

    variable = get_variable(GLODAP_Data,GLODAP_VariableName)
    expocodes = get_variable(GLODAP_Data,"expocode")
    expocodeno = get_variable(GLODAP_Data,"expocodeno")
    cruise = get_variable(GLODAP_Data,"G2cruise")

    # We need to find the index where expocodes == GLODAP_Expocode, and then the
    # expocodeno of that index. We then use CRUISE == expocodeno[index]
    if GLODAP_expocode !== nothing
        cruiseNo = convert(Int32,first(expocodeno[idx]))
        dataIdx = findall(cruise .== cruiseNo)
        variable = variable[dataIdx]
    end
    return variable
end


function loadGLODAPvariables(GLODAP_DIR::String,GLODAP_VariableNames::Vector{String}
    ,GLODAP_expocode::Union{String,String15,Nothing}=nothing
    ,GLODAP_DATA_FILENAME::String="GLODAPv2.2021_Merged_Master_File.mat")
    # Loads multiple variables from GLODAP, for eiher a single cruise, or the 
    # entire GLODAP dataset

    GLODAP_DATAFILE = joinpath(GLODAP_DIR,GLODAP_DATA_FILENAME)
    GLODAP_Data = MatFile(GLODAP_DATAFILE)

    variables = Dict{String,Vector{Float64}}()
    for variableName in GLODAP_VariableNames
        variable = get_variable(GLODAP_Data,variableName)
        expocodes = get_variable(GLODAP_Data,"expocode")
        expocodeno = get_variable(GLODAP_Data,"expocodeno")
        cruise = get_variable(GLODAP_Data,"G2cruise")

        # We need to find the index where expocodes == GLODAP_Expocode, and then the
        # expocodeno of that index. We then use CRUISE == expocodeno[index]
        if GLODAP_expocode !== nothing
            idx = first(findall(expocodes .== GLODAP_expocode))
            cruiseNo = convert(Int32,first(expocodeno[idx]))
            dataIdx = findall(cruise .== cruiseNo)
            variable = variable[dataIdx]
        end

        variables[variableName] = variable
    end
    return variables
end

function loadGLODAPvariables(GLODAP_DIR::String,GLODAP_VariableNames::Vector{String}
    ,GLODAP_expocodes::Union{Vector{String},Vector{String15}}
    ,GLODAP_DATA_FILENAME::String="GLODAPv2.2021_Merged_Master_File.mat")

    GLODAP_DATAFILE = joinpath(GLODAP_DIR,GLODAP_DATA_FILENAME)
    GLODAP_Data = MatFile(GLODAP_DATAFILE)

    variables = Dict{String,Vector{Float64}}()
    idxDict = Dict{String,Vector{Int64}}()
    # Now run through each expocode and extract the indices of observations at
    # the stations we are after.
    for expocode in GLODAP_expocodes
        expocodeList = get_variable(GLODAP_Data,"expocode")
        expocodeno = get_variable(GLODAP_Data,"expocodeno")
        cruise = get_variable(GLODAP_Data,"G2cruise")

        if length(findall(expocodeList .== expocode)) > 0
            idx = first(findall(expocodeList .== expocode))
            cruiseNo = convert(Int64,first(expocodeno[idx]))
            dataIdx = findall(cruise .== cruiseNo)
            idxDict[string(expocode)] = dataIdx
            # Wrapping this in an if statement means that it won't break if no
            # data has been found for that cruise
        end

    end
    goodIdx = reduce(vcat,[indices[2] for indices in idxDict])

    # Create a dict containing all the glodap variables and populate the names,
    # but initialise with empty vectors
    variables = Dict{String,Vector{Float64}}()
    for variableName in GLODAP_VariableNames
        variables[variableName] = get_variable(GLODAP_Data,variableName)[goodIdx]
    end

    return variables
end


function loadGLODAPVariables(GLODAP_DIR::String,GLODAP_VariableNames::Vector{String}
    ;GLODAP_DATA_FILENAME::String="GLODAPv2.2021_Merged_Master_File.mat"
    ,GLODAP_expocodes::Union{Vector{String},Vector{String15},Nothing}=nothing
    ,GLODAP_expocode::Union{String,String15,Nothing}=nothing)
    # We can't use multiple dispatch with keyword arguments, so I'm going to put
    # this wrapper around loadGLODAPVariables to sort the issue. 

    # COME BACK AND FIX THIS AT SOME POINT BECAUSE ITS NOT VERY CLEAN

    if GLODAP_expocodes !== nothing # Structuring logic this way will load all expocodes if both are passed
        println("GLODAP_expocodes is a vector")
        variables = loadGLODAPvariables(GLODAP_DIR,GLODAP_VariableNames
                                       ,GLODAP_expocodes,GLODAP_DATA_FILENAME)
    else
        println("GLODAP_expocode is a string")
        variables = loadGLODAPvariables(GLODAP_DIR,GLODAP_VariableNames
                                       ,GLODAP_expocode,GLODAP_DATA_FILENAME)
    end
    return variables

end

function centralDiff(v::AbstractVector)
    dvF  = diff(v)
    dvB  = reverse(-diff(reverse(v)))
    dx = Vector{AbstractFloat}(undef,length(dvF)+1)
    dx[1] = dvF[1]; dx[end] = dvB[end]
    dx[2:end-1] = (dvF[1:end-1] + dvB[2:end]) / 2

    return dx
end

function gridHorzDistance(GLODAP_latitudes::Vector{Float64}
                        ,GLODAP_longitudes::Vector{Float64}
                        ,latlonGrid)
    # Compute the distance between each station
    dLon = centralDiff(GLODAP_longitudes)
    dLat = centralDiff(GLODAP_latitudes)

    dLat_m = dLat * 111.2
    dLon_m = dLon * 111.2 .* cos.(GLODAP_latitudes*pi/180)

    horzDist_kilometres = sum(sqrt.(dLat_m.^2 + dLon_m.^2)) / length(latlonGrid)
    horzDist_kilometres = fill(horzDist_kilometres,size(latlonGrid))

    return horzDist_kilometres
end

function gridVertDistance(pressureGrid
                         ,GLODAP_pressures=nothing)
    # Because GO-SHIP Easy Ocean grids everything onto a 10m vertical grid, all
    # we need to do is return a distance of 10m
    vertDist_metres = fill(10,size(pressureGrid))
    return vertDist_metres
end

# For some reason which escapes me, I can't type either of the grid___Distance
# functions without breaking them. Something to understand

function createSigmaGrid(sigmaVals::Vector{Float64},numLevels::Int64=600)

    sigmaVals = unique(sort(filter(!isnan,sigmaVals)))
    sigmaStep = convert(Int64,ceil(length(sigmaVals) / numLevels))
    sigmaGrid = sigmaVals[1:sigmaStep:end]
    return sigmaGrid
end

function gridSigDistance(sigmaGrid::Vector{Float64})

    sigmaMeanDist = mean(centralDiff(sigmaGrid))
    sigmaMeanDist = fill(sigmaMeanDist, size(sigmaGrid))

    return sigmaMeanDist
end

function calcScaleFactors(verticalDistance::Vector,horizontalDistance::Vector)

     dP_grid, dL_grid = ndgrid(verticalDistance,horizontalDistance)
     scaleVert = ones(size(dP_grid)) ./ dP_grid
     scaleHorz = ones(size(dL_grid)) ./ dL_grid
     println("Horizontal scale factor: ", mean(scaleHorz))
     println("Vertical scale factor: ", mean(scaleVert))
     return scaleVert, scaleHorz
end

function splitMeanAnom(variable::Vector{Float64})
    varAnom =  variable .- mean(filter(!isnan,variable))
    varMean = mean(filter(!isnan,variable))
    return varMean, varAnom
end


function splitMeanAnom(;obsVariable::Vector{Float64},obsPres::Vector{Float64} #kwargs might not be necessary since this is at the bottom of the call stack
                      ,obsLatLon::Vector{Float64},bgField::Matrix{Float64}
                      ,vertGrid::Vector{Float64},horzGrid::Vector{Float64})

    lowVert = vertGrid[1:end-1]; highVert = vertGrid[2:end]
    lowHorz = horzGrid[1:end-1]; highHorz = horzGrid[2:end]

    varMeanGrid = fill(NaN,length(vertGrid),length(horzGrid))
    varAnom = fill(NaN,length(obsVariable))

    for vertVal in eachindex(highVert)
        @threads for horzVal in eachindex(highHorz) # Still too slow. How do we make this faster
                minVertVal = lowVert[vertVal]; maxVertVal = highVert[vertVal]
                idxVert = findall(minVertVal .< obsPres .< maxVertVal)

                if length(idxVert) > 0

                    minHorzVal = lowHorz[horzVal]; maxHorzVal = highHorz[horzVal]
                    idxHorz = findall(minHorzVal .< obsLatLon .< maxHorzVal)

                    idx = intersect(idxHorz,idxVert)

                    varInBox = obsVariable[idx]
                    goodIdx = findNonNaNIndices(varInBox); varInBox = varInBox[goodIdx]
                    varAnom[idx] = obsVariable[idx] .- bgField[vertVal,horzVal]
                end
        end
    end
    return varAnom
end

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
    goodIdx = findall(convert(Vector{Bool},fill(0,size(variable)) + isnan.(variable)))
    return goodIdx
end

function findNonNaNIndices(variable::Vector)
    goodIdx = findall(convert(Vector{Bool},fill(1,size(variable)) - isnan.(variable)))
    return goodIdx
end

function findNonNaNIndices(variable::Vector,gamma::Vector)
    varGoodIdx = findall(convert(Vector{Bool},fill(1,size(variable)) - isnan.(variable)))
    gamGoodIdx = findall(convert(Vector{Bool},fill(1,size(gamma)) - isnan.(gamma)))
    goodIdx = intersect(varGoodIdx,gamGoodIdx)
    return goodIdx
end

function horzCorrDistanceKilometres(horzLengthDegrees::Vector{Float64}
                               ;latitudes::Union{Vector{Float64},Nothing}=nothing
                               ,meanLatitude::Union{Float64,Nothing}=nothing
                               ,horzCoordinate::Union{String,Nothing}=nothing)

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
    else
        error("\"meanValue\" must be specified to be either \"horzMean\",\"climatology\" or \"scalar\"")
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

function hasDataFlags(variableName::String)
    flagsExist = false
    flaggedVariables = ["G2aou","G2c13","G2c14","G2ccl4","G2cfc113","G2cfc11"
            ,"G2cfc12","G2chla","G2doc","G2fco2","G2h3","G2he3","G2neon"
            ,"G2nitrate","G2nitrite","G2o18","G2oxygen","G2phosphate"
            ,"G2phtsinsitutp","G2salinity","G2sf6","G2silicate","G2talk"
            ,"G2tco2","G2tdn","G2toc"]

    numHits = findall(flaggedVariables.==variableName)
    if length(numHits) > 0
        flagsExist = true
    end
    return flagsExist
end




function removeFlaggedData(variable::Vector{Float64}, variable_fFlag::Vector{Float64})

    outputVar = copy(variable)
    outputVar[variable_fFlag .!= 2] .= NaN
    if length(filter(!isnan,outputVar)) == 0
        outputVar = copy(variable)
        goodIdx = [x==0 || x==2 for x in variable_fFlag]
        outputVar[goodIdx.==0] .= NaN
    end

    return outputVar
end

function createDensitySpaceMask(obsSigma::Vector{Float64}, obsLonLat::Vector{Float64}
                               ,lonLatGrid::Vector{Float64}, sigmaGrid::Vector{Float64}
                               ,pressureSpaceMask::Matrix{Bool})

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

function gridCruisePipeline(;glodapDir::String="/Users/ct6g18/MATLAB/GLODAP"
                            ,goshipDir::String="/Users/ct6g18/MATLAB/GO_SHIP"
                            ,maskMatFile::String="/Users/ct6g18/Julia/glodap_easy_ocean/GLODAP_Easy_Ocean/GOSHIP_MaskStruct.mat"
                            ,sectionName::String
                            ,horzCoordinate::String
                            ,expocode::String
                            ,gridding::String="isobaric"
                            ,meanValue::String="scalar"
                            ,epsilonVal::Float64=0.1
                            ,variableName::String
                            ,plotResults::Bool=false
                            ,autoTruncateMask::Bool=false
                            ,crossValidate::Bool=false
                            ,crossValidationNum::Int64=5)

    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    if gridding != "isobaric" && gridding != "isopycnic"
        error("\"gridding\" must be specified to be either \"isobaric\" or \"isopycnic\"")
    end

    if meanValue != "horzMean" && meanValue != "scalar" && meanValue != "climatology"
        error("\"meanValue\" must be specified to be either \"horzMean\", \"scalar\" or \"climatology\"")
    end

    meanValue == "climatology" ? bgField = readBackgroundField(
    sectionName=sectionName,variableName=variableName) : bgField = nothing

    if crossValidate == false
        crossValidationNum = 0
        bestFactorLen = 1
        bestFactorε = 1
    end

    gridDir, repDir, convDir = loadGoShipDirectories(goshipDir)
    llGrid, prGrid, sectionMask = loadSectionInfo(maskMatFile,sectionName,gridDir)

    isAnException = testExpocodeException(expocode=expocode,variableName=variableName,maskMatfile=maskMatFile)
    println(isAnException)

    if !isAnException
        variable = loadGLODAPVariable(glodapDir,variableName,GLODAP_expocode=expocode)
        griddingVars = loadGLODAPVariables(glodapDir,["G2longitude","G2latitude","G2pressure","G2gamma","G2station"],GLODAP_expocode=expocode)

        varNeedsFlags = hasDataFlags(variableName)
        if varNeedsFlags == true
            println("Removing flagged data")
            varFlagsName::String = variableName * "f"
            variableFlags = loadGLODAPVariable(glodapDir,varFlagsName,GLODAP_expocode=expocode)
            variable = removeFlaggedData(variable,variableFlags)
        end

    else
        println("Cruise "* expocode * " is an exception: loading exception data")
        (variable, griddingVars) = loadExceptionData(expocode=expocode,variableName=variableName,maskMatfile=maskMatFile)
    end

    if variableName == "G2tco2"
        tco2Adj = findGLODAPtco2Adjustment(GLODAP_DIR,expocode=expocode)
        variable .+= tco2Adj
        tco2Adj == 0 ? println("No DIC adjustment necessary") : println("Adjusting DIC by " * string(tco2Adj) * "μmol/kg")
    end

    lon = griddingVars["G2longitude"]
    lat = griddingVars["G2latitude"]
    pressure = griddingVars["G2pressure"]
    sigma = griddingVars["G2gamma"]
    station = griddingVars["G2station"]

    goodIdx = checkVariableExceptions(expocode=expocode,variableName=variableName
    ,maskMatfile=maskMatFile,glodapDir=glodapDir,variable=variable,station=station,pressure=pressure)

    lon = lon[goodIdx]; lat = lat[goodIdx]; pressure = pressure[goodIdx];
    sigma = sigma[goodIdx]; station = station[goodIdx]; variable = variable[goodIdx]


    if horzCoordinate == "longitude"
        XDIR = lon
    elseif horzCoordinate == "latitude"
        XDIR = lat
    end

    if autoTruncateMask
        isPartialCruise = checkPartialCruise(llGrid
                                              ;horzCoordinate=horzCoordinate
                                              ,obsLat=lat,obsLon=lon)
        if isPartialCruise
            println("Cruise did not occcupy full section, truncating mask")
            sectionMask = maskPartialCruise(sectionMask,obsLat=lat,obsLon=lon
                                            ,horzGrid=llGrid,horzCoordinate=horzCoordinate)
        end
    end

    if variableName in ["G2longitude","G2latitude","G2year","G2month"]
        sectionMask .= true
    end

    horzDist = gridHorzDistance(lat,lon,llGrid)
    vertDist = gridVertDistance(prGrid)

    scaleVert, scaleHorz = calcScaleFactors(vertDist,horzDist)

    lenxFactor = checkHorzLenFactor(expocode=expocode,variableName=variableName
                               ,maskMatfile=maskMatFile, griddingType=gridding
                               ,glodapDir=glodapDir)

    len = calcCorrLengths(variable,obsLat=lat,obsLon=lon,obsPres=pressure
    ,presGrid=prGrid,pressureStepNumber=10,verticalSearchRange=100,lenxFactor=lenxFactor)

    P_grid,L_grid = ndgrid(prGrid,llGrid)

    if gridding == "isobaric"
        if crossValidate
            bestFactorLen, bestFactorε = easyDIVACrossValidate(variable=variable
            ,vertVar=pressure,latLon=XDIR,vertGrid=prGrid,horzGrid=llGrid
            ,horzScale=scaleHorz,vertScale=scaleVert,mask=sectionMask
            ,Epsilon=epsilonVal,horzCorrLength=len[2],vertCorrLength=len[1]
            ,horzCoordinate=horzCoordinate)
            # As of now, I haven't passed CrossValidate any extra parameters - instead
            # I'm just going to let it do 5 (default) of each just to keep things simple.
        end
        griddedVarEasyPipeline = easyDIVAGrid(variable=variable,vertVar=pressure
        ,latLon=XDIR,vertGrid=prGrid,horzGrid=llGrid,horzScale=scaleHorz
        ,vertScale=scaleVert,mask=sectionMask,Epsilon=bestFactorε*epsilonVal
        ,horzCorrLength=bestFactorLen*len[2],vertCorrLength=bestFactorLen*len[1]
        ,horzCoordinate=horzCoordinate,meanValue=meanValue,bgField=bgField)
    elseif gridding == "isopycnic"
        lenxPrescribed = checkHorzLenFactor(expocode=expocode
                                                 ,variableName=variableName
                                                 ,maskMatfile=maskMatFile
                                                 ,griddingType=gridding
                                                 ,glodapDir=glodapDir)

        griddedVarEasyPipeline = easyDIVAisopycnal(obsVariable=variable,
        obsSigma=sigma,obsPressure=pressure,obsLat=lat,obsLon=lon,latLon=XDIR
        ,vertGrid=prGrid,horzGrid=llGrid,horzScale=scaleHorz,vertScale=scaleVert
        ,pressureMask=sectionMask,Epsilon=epsilonVal,horzCorrLength=len[2]
        ,vertCorrLength=len[1],horzCorrLengthPrescribed=lenxPrescribed)
    end

    if plotResults == true
        display(heatmap(vec(llGrid),vec(prGrid),griddedVarEasyPipeline, c=:jet1
        ,yflip=true, title=(variableName[3:end] * ", cruise:" * expocode * ", gridding: "* gridding)
        ,xlabel=horzCoordinate,ylabel="Pressure"))
    end

    return griddedVarEasyPipeline
end

function gridSectionPipeline(;glodapDir::String="/Users/ct6g18/MATLAB/GLODAP"
                            ,goshipDir::String="/Users/ct6g18/MATLAB/GO_SHIP"
                            ,maskMatFile::String="/Users/ct6g18/Julia/ExcessHeat/GLODAP_Easy_Ocean/GOSHIP_MaskStruct.mat"
                            ,sectionName::String
                            ,horzCoordinate::String
                            ,convDir::Union{Nothing,String}=nothing
                            ,gridding::String="isobaric"
                            ,fieldMeanVal::String="scalar"
                            ,epsilonVal::Float64=0.01
                            ,autoTruncateMask::Bool=false
                            ,variableName::String)

                            # Check that oxygen has flagging
    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    if gridding != "isobaric" && gridding != "isopycnic"
        error("\"gridding\" must be specified to be either \"isobaric\" or \"isopycnic\"")
    end

    if fieldMeanVal != "horzMean" && fieldMeanVal != "scalar" && fieldMeanVal != "climatology"
        error("\"meanValue\" must be specified to be either \"horzMean\", \"scalar\" or \"climatology\"")
    end

    fieldMeanVal == "climatology" ? bgField = readBackgroundField(
    sectionName=sectionName,variableName=variableName) : bgField = nothing

    if convDir === nothing
        gridDir, repDir, convDir = loadGoShipDirectories(GOSHIP_DIR)
    else
        gridDir, repDir, _ = loadGoShipDirectories(GOSHIP_DIR)
    end

    llGrid, prGrid, sectionMask = loadSectionInfo(MASK_MATFILE,sectionName,gridDir)

    expocodeInfo = listSectionExpocodes(sectionName,convDir)
    expocodes = expocodeInfo[!,"GLODAP Expocode"]

    outputArray = fill(NaN,size(sectionMask,1),size(sectionMask,2),length(expocodes))

    goodExpocodeIdx = fill(1,0)
    for expocode in enumerate(expocodes)
        try
            isAnException = testExpocodeException(expocode=expocode[2],variableName=variableName,maskMatfile=maskMatFile)
            if !isAnException
                variable = loadGLODAPVariable(glodapDir,variableName,GLODAP_expocode=string(expocode[2]))
            else
                (variable, griddingVars) = loadExceptionData(expocode=expocode[2]
                ,variableName=variableName,maskMatfile=maskMatFile)
            end
            goodExpocodeIdx = vcat(goodExpocodeIdx,expocode[1])
        catch
        end
    end
    goodExpocodeIdx = convert(Vector{Int64},goodExpocodeIdx)
    expocodes = expocodes[goodExpocodeIdx,:]

    for expocode in enumerate(expocodes)
        println("Starting on cruise " * expocode[2])

        isAnException = testExpocodeException(expocode=expocode[2],variableName=variableName,maskMatfile=maskMatFile)

        if !isAnException
            variable = loadGLODAPVariable(glodapDir,variableName,GLODAP_expocode=expocode[2])
            griddingVars = loadGLODAPVariables(glodapDir,["G2longitude","G2latitude","G2pressure"
            ,"G2gamma","G2station"],GLODAP_expocode=expocode[2])

            varNeedsFlags = hasDataFlags(variableName)
            if varNeedsFlags == true
                println("Removing flagged data")
                varFlagsName::String = variableName * "f"
                variableFlags = loadGLODAPVariable(glodapDir,varFlagsName,GLODAP_expocode=expocode[2])
                variable = removeFlaggedData(variable,variableFlags)
            end
            if variableName == "G2tco2"
                tco2Adj = findGLODAPtco2Adjustment(GLODAP_DIR,expocode=expocode[2])
                tco2Adj == 0 ? println("No DIC adjustment necessary") : println("Adjusting DIC by " * string(tco2Adj) * "μmol/kg")
                variable .+= tco2Adj
            end

        else
            println("Cruise " * expocode[2]* " is an exception: loading data")
            (variable, griddingVars) = loadExceptionData(expocode=expocode[2]
            ,variableName=variableName,maskMatfile=maskMatFile)
        end

        lon = griddingVars["G2longitude"]
        lat = griddingVars["G2latitude"]
        pressure = griddingVars["G2pressure"]
        sigma = griddingVars["G2gamma"]
        station = griddingVars["G2station"]

        goodIdx = checkVariableExceptions(expocode=expocode[2],variableName=variableName
        ,maskMatfile=maskMatFile,glodapDir=glodapDir,variable=variable,station=station,pressure=pressure)

        lon = lon[goodIdx]; lat = lat[goodIdx]; pressure = pressure[goodIdx];
        sigma = sigma[goodIdx]; station = station[goodIdx]; variable = variable[goodIdx]

        if horzCoordinate == "longitude"
            XDIR = lon
        elseif horzCoordinate == "latitude"
            XDIR = lat
        end

        if autoTruncateMask
            _, _, sectionMask = loadSectionInfo(MASK_MATFILE,sectionName,gridDir)
            # Need to overwrite sectionMask on each iteration or we will have problems
            isPartialCruise = checkPartialCruise(llGrid
                                                  ;horzCoordinate=horzCoordinate
                                                  ,obsLat=lat,obsLon=lon)
            if isPartialCruise
                println("Cruise did not occcupy full section, truncating mask")
                sectionMask = maskPartialCruise(sectionMask,obsLat=lat,obsLon=lon
                                                ,horzGrid=llGrid,horzCoordinate=horzCoordinate)
            end
        end

        horzDist = gridHorzDistance(lat,lon,llGrid)
        vertDist = gridVertDistance(prGrid)

        scaleVert, scaleHorz = calcScaleFactors(vertDist,horzDist)

        lenxFactor = checkHorzLenFactor(expocode=expocode[2],variableName=variableName
                                          ,maskMatfile=maskMatFile, griddingType=gridding
                                          ,glodapDir=glodapDir)

        len = calcCorrLengths(variable,obsLat=lat,obsLon=lon,obsPres=pressure
        ,presGrid=prGrid,pressureStepNumber=10,verticalSearchRange=100,lenxFactor=lenxFactor)

        P_grid,L_grid = ndgrid(prGrid,llGrid)

        if gridding == "isobaric"
            griddedVarEasyPipeline = easyDIVAGrid(variable=variable,vertVar=pressure,latLon=XDIR
            ,horzCoordinate=horzCoordinate,vertGrid=prGrid,horzGrid=llGrid,horzScale=scaleHorz
            ,vertScale=scaleVert,mask=sectionMask,Epsilon=epsilonVal,horzCorrLength=len[2]
            ,vertCorrLength=len[1],meanValue=fieldMeanVal,bgField=bgField)
        elseif gridding == "isopycnic"
            griddedVarEasyPipeline = easyDIVAisopycnal(obsVariable=variable,obsSigma=sigma,obsPressure=pressure
            ,pressureMask=sectionMask,obsLat=lat,obsLon=lon,latLon=XDIR,vertGrid=prGrid,horzGrid=llGrid
            ,horzScale=scaleHorz,vertScale=scaleVert,Epsilon=epsilonVal,horzCorrLength=len[2],vertCorrLength=len[1])
        end
    outputArray[:,:,expocode[1]] = griddedVarEasyPipeline
    end

    println(goodExpocodeIdx)
    return outputArray
end


function testExpocodeException(;expocode::Union{String,String15},variableName::String
                               ,maskMatfile::String)
    # Take an expocode, check whether that expocode is one of our exceptions.
    # If it is, then check the variable we are after is in the list of variables
    # contained in the exception data.
    # If it is, return true. If not, return false
    exceptionDir, _  = splitdir(maskMatfile)
    expocodes = joinpath(exceptionDir,"Exceptions/ExpocodeList.csv")
    expocodeList = Vector{String}(readdlm(expocodes,',',String;header=false)[2:end,2])

    if length(intersect(expocodeList,[expocode])) > 0
        isException = true
    else
        isException = false
    end

    if isException == true # Now open the .mat file associated with our expocode
        # see if we can find the variable we're after
        exceptionData = joinpath(exceptionDir,"Exceptions/ExceptionData",expocode*".mat")
        SectionFile = MatFile(exceptionData)
        exceptionDataVariables = variable_names(SectionFile)

        if length(intersect(exceptionDataVariables,[variableName])) < 1
            isException = false
        end
    end

    return isException
end

function loadExceptionData(;expocode::Union{String,String15},variableName::String
                               ,maskMatfile::String)

    exceptionDir, _  = splitdir(maskMatfile)
    exceptionData = joinpath(exceptionDir,"Exceptions/ExceptionData",expocode*".mat")
    SectionFile = MatFile(exceptionData)
    exceptionDataVariables = variable_names(SectionFile)

    variable = get_variable(SectionFile,variableName)

    gridding_varNames = ["G2longitude","G2latitude","G2pressure","G2gamma","G2station"]
    griddingVars= Dict{String,Vector{Float64}}()
    for variableName in gridding_varNames
        griddingVars[variableName]  = get_variable(SectionFile,variableName)
    end

    return (variable, griddingVars)

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

function checkVariableExceptions(;expocode::Union{String,String15},variableName::String
                               ,maskMatfile::String
                               ,variable::Union{Vector{Float64},Nothing} = nothing
                               ,station::Union{Vector{Float64},Nothing} = nothing
                               ,pressure::Union{Vector{Float64},Nothing} = nothing
                               ,glodapDir::String="/home/ct/MATLAB/GLODAP")

    if variable === nothing
        vars = loadGLODAPVariables(glodapDir,[variableName,"G2station","G2pressure"]
        ,GLODAP_expocode=expocode)
        variable = vars[variableName]
        station = vars["G2station"]
        pressure = vars["G2pressure"]
    end

    exceptionDir, _  = splitdir(maskMatfile)
    variableExceptionData = joinpath(exceptionDir,"Exceptions/variableExceptions.csv")

    expocodeList = Vector{String}(readdlm(variableExceptionData,',';header=false)[2:end,1])
    variableList = Vector{String}(readdlm(variableExceptionData,',';header=false)[2:end,2])
    stationList = Vector{Int64}(readdlm(variableExceptionData,',';header=false)[2:end,3])
    minPressureList = Vector{Int64}(readdlm(variableExceptionData,',';header=false)[2:end,4])
    maxPressureList = Vector{Int64}(readdlm(variableExceptionData,',';header=false)[2:end,5])

    requiresFlagging = false
    expocodeIdx = findall(expocodeList .== [expocode])
    variableIdx = findall(variableList .== [variableName])

    colIdx = intersect(expocodeIdx,variableIdx) # colIdx tells us where to look
    # for data to remove, and whether we  need to return anything
    if length(colIdx) > 0
        println("Cruise " * expocode * " contains excluded data for variable " * variableName)
        requiresFlagging = true
    else
        return collect(1:length(variable)) # Return all bottle numbers
    end

    goodIdx = collect(1:length(variable))

    for idx in colIdx
        criteriaMatched = fill(0,size(goodIdx))
        stnVal = stationList[idx]
        maxPrsVal = maxPressureList[idx]; minPrsVal = minPressureList[idx]
        criteriaMatched[station .== stnVal] .+=1
        criteriaMatched[pressure .< maxPrsVal] .+=1
        criteriaMatched[pressure .> minPrsVal] .+=1

        goodIdx[criteriaMatched .== 3] .= -1
    end
    goodIdx = goodIdx[goodIdx .> 0]
    # Since we know we should have removed something, we can test its worked and
    # return an error if not pretty trivially
    if length(goodIdx) == length(variable)
        error("checkVariableExceptions is not behaving as expected")
    end

    return goodIdx

end

function findGLODAPtco2Adjustment(GLODAP_DIR::String;
    GLODAP_AdjTable::String="AdjustmentTable.csv",expocode::Union{String,String15})
    
    # This will check whether the expocode we are looking at is in Jens Muller's 
    # recommended adjustments and adjust up here.
    JensList = ["316N19941201","316N19950124","316N19950310","316N19950423",
    "316N19950611","316N19950715","316N19950829","316N19951111","316N19951202"]

    if expocode in JensList
        println("Expocode in Jens Muller's recommended adjustment list, returning 1.7")
        return 1.7
    end

    GLODAP_AdjTable = joinpath(GLODAP_DIR,GLODAP_AdjTable)

    tableCSV = CSV.read(GLODAP_AdjTable,DataFrame)

    expocodeTable = tableCSV[!,"cruise_expocode"]
    tco2AdjTable = tableCSV[!,"tco2_adj"]

    if length(findall(expocodeTable .== expocode)) == 0
        println("No adjustment found, returning 0.0")
        return 0.0
    end

    adjustment = tco2AdjTable[findall(expocodeTable .== expocode)][1]
    if adjustment < -665
        adjustment = 0.0
    end
    return adjustment
end

function checkHorzLenFactor(;expocode::Union{String,String15},variableName::String
                               ,maskMatfile::String, griddingType::String
                               ,glodapDir::String="/Users/ct6g18/MATLAB/GLODAP")

    exceptionDir, _  = splitdir(maskMatfile)
    variableExceptionData = joinpath(exceptionDir,"Exceptions/horzLenExceptions.csv")

    exceptionDataFrame = CSV.read(variableExceptionData,DataFrame)
    EDFsubset = exceptionDataFrame[(exceptionDataFrame.Expocode .== expocode) .&
    (exceptionDataFrame.Variable .== variableName) .&
    (exceptionDataFrame.Gridding .== griddingType), :]
    if prod(size(EDFsubset)) > 0
        println("Expocode contains horizonal correlation length exception for selected variable and gridding")
        factor = EDFsubset[!,"Factor"]
        return factor[1]
    elseif griddingType == "isobaric"
        return 1
    elseif griddingType == "isopycnic"
        return 1.0
    end
end

function checkPartialCruise(horzGrid::Vector{Float64};horzCoordinate::String
                            ,obsLat::Vector{Float64},obsLon::Vector{Float64})

    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    if horzCoordinate == "longitude"
         obsXval = obsLon
         obsLon = matchLonConventions(horzGrid,obsLon)
    else
         obsXval = obsLat
    end

    maxGridVal = maximum(horzGrid)
    minGridVal = minimum(horzGrid)

    maxObsVal = maximum(obsXval)
    minObsVal = minimum(obsXval)

    ΔminVal = minGridVal - minObsVal # If minGrid << minObs, this will be +ve
    println("ΔminVal = " * string(ΔminVal))
    ΔmaxVal = maxGridVal - maxObsVal # If maxGrid << maxObs, this will be +ve
    println("ΔmaxVal = " * string(ΔmaxVal))

    isPartialCruise = false
    if abs(ΔmaxVal) > 2 || abs(ΔminVal) > 2 # Abs because we might get NH/SH problems
        isPartialCruise = true
    end

    return isPartialCruise
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

    gridDir, repDir, convDir = loadGoShipDirectories(goshipDir)
    horzGrid, prGrid, sectionMask = loadSectionInfo(maskMatFile,sectionName,gridDir)
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

function maskPartialCruise(mask::Matrix{Bool};obsLat::Vector{Float64}
                           ,obsLon::Vector{Float64},horzGrid::Vector{Float64}
                           ,horzCoordinate::String)


    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    if horzCoordinate == "longitude"
         obsXval = obsLon
         obsLon = matchLonConventions(horzGrid,obsLon)
    else
         obsXval = obsLat
    end

    maxObsVal = maximum(obsXval)
    minObsVal = minimum(obsXval)

    llGrid = Matrix{Float64}(undef,size(mask'))
    [llGrid[:,i] = horzGrid for i in 1:size(mask,1)]
    llGrid = llGrid' # Transposing for performance. Probably unimportant & counterproductive

    truncatedMask = copy(mask)
    truncatedMask[llGrid .< minObsVal] .= false
    truncatedMask[llGrid .> maxObsVal] .= false

    return truncatedMask
end


function maskPartialSectionPipeline(;glodapDir::String="/home/ct/MATLAB/GLODAP"
                                  ,goshipDir::String="/home/ct/MATLAB/GO_SHIP"
                                  ,maskMatFile::String="/home/ct/Julia/GLODAP_Easy_Ocean/GOSHIP_MaskStruct.mat"
                                  ,sectionName::String
                                  ,horzCoordinate::String
                                  ,variable::Array{Float64})

    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    gridDir, repDir, convDir = loadGoShipDirectories(GOSHIP_DIR)
    llGrid, prGrid, sectionMask = loadSectionInfo(MASK_MATFILE,sectionName,gridDir)

    expocodeInfo = listSectionExpocodes(sectionName,convDir)
    expocodes = expocodeInfo[!,"GLODAP Expocode"]

    outputArray = fill(NaN,size(sectionMask,1),size(sectionMask,2),length(expocodes))

    goodExpocodeIdx = fill(1,0)
    for expocode in enumerate(expocodes)
        try
            variable = loadGLODAPVariable(glodapDir,variableName,GLODAP_expocode=string(expocode[2]))
            goodExpocodeIdx = vcat(goodExpocodeIdx,expocode[1])
        catch
        end
    end
    goodExpocodeIdx = convert(Vector{Int64},goodExpocodeIdx)
    expocodes = expocodes[goodExpocodeIdx,:]

    for expocode in enumerate(expocodes)
        println("Starting on cruise " * expocode[2])

        isAnException = testExpocodeException(expocode=expocode[2],variableName=variableName,maskMatfile=maskMatFile)

        if !isAnException
            griddingVars = loadGLODAPVariables(glodapDir,["G2longitude","G2latitude"],GLODAP_expocode=expocode[2])
        else
            println("Cruise " * expocode[2]* " is an exception: loading data")
            (_, griddingVars) = loadExceptionData(expocode=expocode[2]
            ,variableName=variableName,maskMatfile=maskMatFile)
        end

        lon = griddingVars["G2longitude"]
        lat = griddingVars["G2latitude"]


        _, _, sectionMask = loadSectionInfo(MASK_MATFILE,sectionName,gridDir)
        # Need to overwrite sectionMask on each iteration or we will have problems
        isPartialCruise = checkPartialCruise(llGrid;horzCoordinate=horzCoordinate
                                            ,obsLat=lat,obsLon=lon)
        if isPartialCruise
            println("Cruise did not occcupy full section, truncating outputs")
            sectionMask = maskPartialCruise(sectionMask,obsLat=lat,obsLon=lon
                                            ,horzGrid=lonGrid,horzCoordinate=horzCoordinate)
            variable[:,:,expocode[1]] = sectionMask .* variable[:,:,expocode[1]]
        end

    end

    return variable
end

function checkSectionPartialCruises(;glodapDir::String="/home/ct/MATLAB/GLODAP"
                                  ,goshipDir::String="/home/ct/MATLAB/GO_SHIP"
                                  ,maskMatFile::String="/home/ct/Julia/GLODAP_Easy_Ocean/GOSHIP_MaskStruct.mat"
                                  ,sectionName::String
                                  ,horzCoordinate::String)

    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    gridDir, repDir, convDir = loadGoShipDirectories(GOSHIP_DIR)
    llGrid, prGrid, sectionMask = loadSectionInfo(MASK_MATFILE,sectionName,gridDir)

    expocodeInfo = listSectionExpocodes(sectionName,convDir)
    expocodes = expocodeInfo[!,"GLODAP Expocode"]


    goodExpocodeIdx = fill(1,0)
    for expocode in enumerate(expocodes)
        try
            variable = loadGLODAPVariable(glodapDir,"G2longitude",GLODAP_expocode=string(expocode[2]))
            goodExpocodeIdx = vcat(goodExpocodeIdx,expocode[1])
        catch
        end
    end
    goodExpocodeIdx = convert(Vector{Int64},goodExpocodeIdx)
    expocodes = expocodes[goodExpocodeIdx,:]

    truthArray = fill(false,length(expocodes))

    for expocode in enumerate(expocodes)
        println("Checking for partial occupation during cruise " * expocode[2])

        isAnException = testExpocodeException(expocode=expocode[2],variableName="G2longitude"
                                             ,maskMatfile=maskMatFile)
        if !isAnException
            griddingVars = loadGLODAPVariables(glodapDir,["G2longitude","G2latitude"],GLODAP_expocode=expocode[2])
        else
            println("Cruise " * expocode[2]* " is an exception: loading data")
            (_, griddingVars) = loadExceptionData(expocode=expocode[2]
            ,variableName="G2longitude",maskMatfile=maskMatFile)
        end

        lon = griddingVars["G2longitude"]
        lat = griddingVars["G2latitude"]


        _, _, sectionMask = loadSectionInfo(MASK_MATFILE,sectionName,gridDir)
        # Need to overwrite sectionMask on each iteration or we will have problems
        isPartialCruise = checkPartialCruise(llGrid;horzCoordinate=horzCoordinate
                                            ,obsLat=lat,obsLon=lon)
        if isPartialCruise
            truthArray[expocode[1]] = true
        end

    end

    return truthArray
end

function concatenateCruises(expocodes::Vector{String};GLODAP_DIR::String="/Users/ct6g18/MATLAB/GLODAP"
                           ,GLODAP_DATA_FILENAME::String="GLODAPv2.2021_Merged_Master_File.mat"
                           ,stationRanges::Union{Vector{Vector{Int64}},Nothing}=nothing)
    #= Take two (or more) expocodes, merge them together into a single dict which
    will contain all the observations in each cruise, or a range of the observations
    in each of the cruises based on station numbers. This can then be passed to
    the writeCruiseException function to create an exception .mat file for general
    use.

    This can also be used to take a range of stations from a single cruise in order
    to write an exception.
    =#

     if stationRanges !== nothing && length(stationRanges) != length(expocodes)
         error("If station ranges are specified, you must give the same number of
         station ranges and expocodes")
     end



    GLODAP_DATAFILE = joinpath(GLODAP_DIR,GLODAP_DATA_FILENAME)
    GLODAP_Data = MatFile(GLODAP_DATAFILE)


    # Now run through each expocode and extract the indices of observations at
    # the stations we are after.
    idxDict = Dict{String,Vector{Int64}}()
    for expocode in enumerate(expocodes)
        expocodeList = get_variable(GLODAP_Data,"expocode")
        expocodeno = get_variable(GLODAP_Data,"expocodeno")
        cruise = get_variable(GLODAP_Data,"G2cruise")

        idx = first(findall(expocodeList .== expocode[2]))
        cruiseNo = convert(Int64,first(expocodeno[idx]))
        dataIdx = findall(cruise .== cruiseNo)

        idxDict[string(expocode[2])] = dataIdx

        if stationRanges !== nothing
            stations = get_variable(GLODAP_Data,"G2station")
            cruiseStations = stations[dataIdx]
            incStn = findall(convert(Vector{Bool},sum([cruiseStations .== stnVal for stnVal in stationRanges[expocode[1]]])))
            idxDict[string(expocode[2])] = dataIdx[incStn]
        end

    end
    goodIdx = reduce(vcat,[indices[2] for indices in idxDict])

    # Create a dict containing all the glodap variables and populate the names,
    # but initialise with empty vectors
    variables = Dict{String,Vector{Float64}}()
    for variableName in variable_names(GLODAP_Data)[2:end-2] #Ignore DOI, expocode, expocodeno
        variables[variableName] = get_variable(GLODAP_Data,variableName)[goodIdx]
    end


    return variables
end



function writeCruiseException(;outputExpocode::String,cruiseDict::Dict{String, Vector{Float64}}
                             ,yearString::String="0000" # Can be manually fixed afterwards if needs be
                             ,exceptionMatFileLocation::String="/Users/ct6g18/Julia/ExcessHeat/glodap_easy_ocean/Exceptions/ExceptionData"
                             ,exceptionListCSVFile::String="/Users/ct6g18/Julia/ExcessHeat/glodap_easy_ocean/Exceptions/ExpocodeList.csv")
    #=  When we join cruises, we are going to create a matfile which contains
    every observation from each of the cruises. That way, when we test for an
    exception, no matter the choice of variable, we will load in the joined
    variables =#

    outputFilename = outputExpocode * ".mat"
    outputDatafile = joinpath(exceptionMatFileLocation,outputFilename)

    mf = MatFile(outputDatafile,"w")
    for varName in keys(cruiseDict)
        tmpVar = cruiseDict[varName]
        put_variable(mf,varName,tmpVar)
    end
    close(mf)

    # Now we need to add the expocode we've written to the ExpocodeList.csv file
    # which contains our exceptions.


    expocodeList = Vector{String}(readdlm(exceptionListCSVFile,',';header=false)[2:end,2])

    if sum(expocodeList .== outputExpocode) < 1 # If our expocode isnt found in the list of exceptions
        expocodeListStr = yearString * "," * outputExpocode * "," * outputExpocode
        println("Appending \""* expocodeListStr * "\" to " * exceptionListCSVFile)
        expectionListIOStream = open(exceptionListCSVFile,"a")
        write(expectionListIOStream, expocodeListStr);
        close(expectionListIOStream)
    end
    # May need to force a new line character to the end of the str

    return nothing
end

function createBackgroundField(;variable::String,sectionName::String
                              ,goshipDir::String="/home/ct/MATLAB/GO_SHIP"
                              ,glodapDir::String="/home/ct/MATLAB/GLODAP"
                              ,horzCoordinate::String
                              ,horzLenFactor::Number=100
                              ,vertLenFactor::Number=1
                              ,maskMatFile::String="/home/ct/Julia/GLODAP_Easy_Ocean/GOSHIP_MaskStruct.mat")

    if variable == "G2tco2" || variable == "DIC" || variable == "TCO2"
        climatology = Dataset("/home/ct/Julia/GLODAP_Easy_Ocean/Climatologies/GLODAPv2.2016b.TCO2.nc")
        clim_var = climatology["TCO2"][:]
    elseif variable == "G2theta" || variable == "Temperature" || variable == "TMP"
        climatology = Dataset("/home/ct/Julia/GLODAP_Easy_Ocean/Climatologies/GLODAPv2.2016b.temperature.nc")
        clim_var = climatology["temperature"][:]
    elseif variable == "G2salinity" || variable == "Salinity" || variable == "SAL"
        climatology = Dataset("/home/ct/Julia/GLODAP_Easy_Ocean/Climatologies/GLODAPv2.2016b.salinity.nc")
        clim_var = climatology["salinity"][:]
    end

    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    replace!(clim_var,missing=>NaN)

    clim_lon = climatology["lon"][:]
    clim_lat = climatology["lat"][:]
    clim_dep = climatology["Depth"][:]

    llGrid, prGrid, mask = loadSectionInfo(MASK_MATFILE,sectionName,gridDir)
    _, _, convDir = loadGoShipDirectories(GOSHIP_DIR)
    expocodes = listSectionExpocodes(sectionName,convDir)[:,"GLODAP Expocode"]

    vars = loadGLODAPVariables(GLODAP_DIR,["G2theta","G2tco2","G2salinity"
    ,"G2gamma","G2latitude","G2longitude","G2pressure"]
    ,GLODAP_expocodes=expocodes)



    lon = vars["G2longitude"]
    lat = vars["G2latitude"]

    if sectionName[1] == 'P'
         lon = matchLonConventions(clim_lon,lon)
         lonOffset = 0
    else
         lonOffset = 200
    end

    ll = unique!(DataFrame(lon=lon,lat=lat)) # Find unique locations, ie. stations
    obsLonVals = ll[!,"lon"]
    obsLatVals = ll[!,"lat"]

    cLonGrid = fill(0.0,length(clim_lon),length(obsLonVals))
    [cLonGrid[:,i] = clim_lon .- obsLonVals[i] .- lonOffset for i in eachindex(obsLonVals)]

    k = abs.(cLonGrid)
    lonIdx = fill(0,length(obsLonVals))
    [lonIdx[i] = argmin(k[:,i]) for i in 1:length(obsLonVals)]

    cLatGrid = fill(0.0,length(clim_lat),length(obsLatVals))
    [cLatGrid[:,i] = clim_lat .- obsLatVals[i] for i in eachindex(obsLatVals)]

    k = abs.(cLatGrid)
    latIdx = fill(0,length(obsLatVals))
    [latIdx[i] = argmin(k[:,i]) for i in eachindex(obsLatVals)]


    llG = unique!(DataFrame(lon=lonIdx,lat=latIdx)) # Find unique locations, ie. grid cells
    lons = llG[!,"lon"]
    lats = llG[!,"lat"]
    # These are the same length so we can continue to just use lons as a counter

    # Grid is on 33 levels
    bgField = convert(Matrix{Union{Float64,Missing}},fill(NaN,length(lons),33))

    [bgField[i,:] = clim_var[lons[i],lats[i],:] for i in eachindex(lons)]

    #return (clim_lon[lons],clim_lat[lats],bgField)

    oLonGrid = fill(NaN,size(bgField))
    oLatGrid = fill(NaN,size(bgField))
    oPrsGrid = fill(NaN,size(bgField))

    [oLonGrid[:,i] = clim_lon[lons] for i in 1:length(clim_dep)]
    [oLatGrid[:,i] = clim_lat[lats] for i in 1:length(clim_dep)]
    [oPrsGrid[i,:] = clim_dep for i in 1:length(lons)]

    zVar = bgField[:];  goodIdx = findNonNaNIndices(zVar)
    zVar = convert(Vector{Float64},zVar)[goodIdx]

    w = oLatGrid[:][goodIdx]
    x = oLonGrid[:][goodIdx]
    y = oPrsGrid[:][goodIdx]

    P_grid,L_grid = ndgrid(clim_dep,lons)

    horzDist = gridHorzDistance(w,x,llGrid)
    vertDist = gridVertDistance(prGrid)
    scaleVert, scaleHorz = calcScaleFactors(vertDist,horzDist)

    (clVert, clHorz) = calcCorrLengths(zVar,obsLat=w,obsLon=x
    ,obsPres=y,presGrid=prGrid)

    # There should be some sort of call to matchLonConventions here I think
    # instead of just subtracting 200 from all longitudes carelessly to undo our fuckery.
    if horzCoordinate == "longitude"
        divaGridClim = easyDIVAGrid(variable=zVar,vertVar=y,latLon=x.-lonOffset
        ,horzCoordinate=horzCoordinate,meanValue="scalar",vertGrid=prGrid,horzGrid=llGrid
        ,horzScale=scaleHorz,vertScale=scaleVert,mask=mask,Epsilon=0.1,horzCorrLength=horzLenFactor*clHorz
        ,vertCorrLength = clVert*vertLenFactor)
    elseif horzCoordinate == "latitude"
        divaGridClim = easyDIVAGrid(variable=zVar,vertVar=y,latLon=w
        ,horzCoordinate=horzCoordinate,meanValue="scalar",vertGrid=prGrid,horzGrid=llGrid
        ,horzScale=scaleHorz,vertScale=scaleVert,mask=mask,Epsilon=0.1,horzCorrLength=horzLenFactor*clHorz
        ,vertCorrLength = clVert*vertLenFactor)
    end

    return divaGridClim

end

function writeBackgroundField(;variable::Matrix{Float64},sectionName::String
                             ,variableName::String
                             ,outputPathStr::String="/home/ct/Julia/GLODAP_Easy_Ocean/BackgroundFields/")
    # This function takes a background field on the same grid as the mask we
    # plan on using and writes it to disk.

    # I've only implemented functionality to use a climatological background field
    # for DIC, salinity and temperature. As a result, we can only write if the
    # variable name is one of these. For convenience, I'm going to force them to
    # be written in the GLODAP nomenculature ie. G2tco2.
    if variableName != "G2tco2"  && variableName != "G2theta" && variableName != "G2salinity"
        error("variableName must be either \"G2tco2\", \"G2theta\" or \"G2salinity\"")
    end

    fileName = sectionName * "_" * variableName * ".nc"

    fileNameFull = joinpath(outputPathStr,fileName)

    ds = Dataset(fileNameFull,"c")

    defDim(ds,"Vertical_Cell_No",size(variable,1))
    defDim(ds,"Horizontal_Cell_No",size(variable,2))

    titleStr= "Background field for " * variableName * " for " * sectionName
    ds.attrib["Title"] = titleStr

    v1 = defVar(ds,variableName,Float64,("Horizontal_Cell_No","Vertical_Cell_No"))

    v1[:,:] = variable'

    ds.attrib["Comments"] = "Generated using the GLODAP Climatology"
    ds.attrib["Creator"] = "Charles Turner"

    close(ds)

    return nothing

end

function readBackgroundField(;sectionName::String,variableName::String)
    # This function will go looking for the background field and load it
    if variableName != "G2tco2"  && variableName != "G2theta" && variableName != "G2salinity"
        error("variableName must be either \"G2tco2\", \"G2theta\" or \"G2salinity\"")
    end

    fileName = sectionName * "_" * variableName * ".nc"
    pathStr = "/home/ct/Julia/GLODAP_Easy_Ocean/BackgroundFields/"

    fileNameFull = joinpath(pathStr,fileName)

    ds = Dataset(fileNameFull)

    field = ds[variableName][:]
    field = convert(Matrix{Float64},field')

    return field

end


function writeOutputField(;variable::Array{Float64,3},sectionName::String
                         ,variableName::String
                         ,glodapDir::String="/home/ct/MATLAB/GLODAP"
                         ,horzCoordinate::String
                         ,outputPath::String="/home/ct/Julia/RawOutputFields"
                         ,GLODAP_Datafile::String="GLODAPv2.2021_Merged_Master_File.mat")
    # This function takes an output field and writes it to disk. It is based on
    # the background field functionality

    fileName = sectionName * "_" * variableName * ".nc"

    fileNameFull = joinpath(outputPath,fileName)


    isfile(fileNameFull) ?
    error("File you wish to write to already exists") : nothing

    ds = Dataset(fileNameFull,"c")

    defDim(ds,horzCoordinate,size(variable,1))
    defDim(ds,"Pressure",size(variable,2))
    defDim(ds,"Repeat_Number",size(variable,3))

    titleStr= "Gridded " * variableName * " field for " * sectionName
    ds.attrib["Title"] = titleStr

    v1 = defVar(ds,variableName,Float64,(horzCoordinate,"Pressure","Repeat_Number"))

    v1[:,:,:] = variable

    ds.attrib["Comments"] = "Generated using the GLODAP dataset and gridding data from GO-SHIP Easy Ocean"
    ds.attrib["Creator"] = "Charles Turner"

    close(ds)

    return nothing

end

function writeOutputFields(;variables::Vector{Array{Float64,3}},sectionName::String
                         ,variableNames::Vector{String}
                         ,glodapDir::String="/home/ct/MATLAB/GLODAP"
                         ,horzCoordinate::String
                         ,outputPath::String="/home/ct/Julia/RawOutputFields"
                         ,GLODAP_Datafile::String="GLODAPv2.2021_Merged_Master_File.mat")
    # This function takes an output field and writes it to disk. It is based on
    # the background field functionality

    fileName = sectionName * ".nc"

    fileNameFull = joinpath(outputPath,fileName)


    isfile(fileNameFull) ?
    error("File you wish to write to already exists") : nothing

    ds = Dataset(fileNameFull,"c")

    defDim(ds,horzCoordinate,size(variables[1],1))
    defDim(ds,"Pressure",size(variables[1],2))
    defDim(ds,"Repeat_Number",size(variables[1],3))

    titleStr= "Gridded fields for " * sectionName
    ds.attrib["Title"] = titleStr

    for variable in enumerate(variableNames)
        v = defVar(ds,variable[2],Float64,(horzCoordinate,"Pressure","Repeat_Number"))
        v[:,:,:] = variables[variable[1]]
    end

    ds.attrib["Comments"] = "Generated using the GLODAP dataset and gridding data from GO-SHIP Easy Ocean."
    ds.attrib["Creator"] = "Charles Turner"

    close(ds)

    return nothing

end

function calcGLODAP_Date(glodapYear::Vector{Float64},glodapMonth::Vector{Float64}
                        ,glodapDay::Vector{Float64})

    glodapMonth .-= 1
    
    glodapDate = glodapYear + (glodapMonth ./ 12) + (glodapDay ./ 30)

    return glodapDate
end

function splitMeanAnom(;obsVariable::Vector{Float64},obsPres::Vector{Float64} #kwargs might not be necessary since this is at the bottom of the call stack
                      ,bgField::Vector{Float64},vertGrid::Vector{Float64})

    lowVert = vertGrid[1:end-1]; highVert = vertGrid[2:end]

    varMean= fill(NaN,length(obsVariable))
    varAnom = fill(NaN,length(obsVariable))

    for val in eachindex(highVert)
        bgVal = bgField[val]
        minVal = lowVert[val]; maxVal = highVert[val]
        idx = findall(minVal .< obsPres .< maxVal)

        varInBox = obsVariable[idx]
        goodIdx = findNonNaNIndices(varInBox); varInBox = varInBox[goodIdx]
        varAnom[idx] .= obsVariable[idx] .- bgVal
        varMean[idx] .= bgVal
    end

    return varMean, varAnom
end

function expocodeFromG2cruise(;GLODAPdir::String= "/Users/ct6g18/MATLAB/GLODAP"
    ,GLODAPdataFilename::String="GLODAPv2.2021_Merged_Master_File.mat"
    ,G2cruise::Union{Float64,Vector{Float64}})

    if length(G2cruise) > 1
        G2cruise = G2cruise[1]
    end

    GLODAP_DATAFILE = joinpath(GLODAPdir,GLODAPdataFilename)
    GLODAP_Data = MatFile(GLODAP_DATAFILE)

    expocodes = get_variable(GLODAP_Data,"expocode")
    expocodeno = get_variable(GLODAP_Data,"expocodeno")

    expocodeIdx = findall(expocodeno .== G2cruise)

    expocode = expocodes[expocodeIdx]
    return expocode[1]

end

function matchLonConventions(lonGrid::Vector{Float64}, lonVal::Float64)
    if maximum(lonGrid) > 180 && lonVal < 0
        lonVal += 360
    end
    return lonVal
end

