function checkPartialCruise(horzGrid::Vector{Float64};horzCoordinate::String
                            ,obsLat::Vector{Float64},obsLon::Vector{Float64})
    # Check if the cruise we are gridding didn't occupy the full section. This 
    # is a simple version which will just look at the recorded latitude and 
    # longitude values
    obsXval = legacy_obs_xvals(horzCoordinate, horzGrid, obsLat, obsLon)
    deltas = legacy_partial_cruise_deltas(horzGrid, obsXval)

    println("ΔminVal = " * string(deltas.ΔminVal))
    println("ΔmaxVal = " * string(deltas.ΔmaxVal))

    isPartialCruise = false
    if abs(deltas.ΔmaxVal) > 2 || abs(deltas.ΔminVal) > 2 # Abs because we might get NH/SH problems
        isPartialCruise = true
    end

    return isPartialCruise
end

function maskPartialCruise(mask::Matrix{Bool};obsLat::Vector{Float64}
                           ,obsLon::Vector{Float64},horzGrid::Vector{Float64}
                           ,horzCoordinate::String)
    # If a cruise is a partial cruise, then mask out the regions where we don't 
    # have any observations: otherwise it gets filled in with nonsense
    obsXval = legacy_obs_xvals(horzCoordinate, horzGrid, obsLat, obsLon)
    deltas = legacy_partial_cruise_deltas(horzGrid, obsXval)

    return legacy_truncate_mask(mask, horzGrid, deltas.minObsVal, deltas.maxObsVal)
end


function maskPartialSectionPipeline(;GLODAP_DIR::Union{String,Nothing}=nothing
                                    ,GOSHIP_DIR::Union{String,Nothing}=nothing
                                    ,MASK_MATFILE::Union{String,Nothing}=nothing
                                    ,EXCEPTIONS_DIR::Union{String,Nothing}=nothing
                                    ,EXCEPTIONS_FILENAME::Union{String,Nothing}=nothing
                                    ,sectionName::String
                                    ,horzCoordinate::String
                                    ,variable::Array{Float64})
    # Chain together all the functions needed to mask a partial cruise out
    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end


    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GOSHIP_DIR === nothing ? GOSHIP_DIR = readDefaults()["GOSHIP_DIR"] : nothing
    MASK_MATFILE === nothing ? MASK_MATFILE = readDefaults()["MASK_MATFILE"] : nothing
    EXCEPTIONS_DIR === nothing ? EXCEPTIONS_DIR = readDefaults()["EXCEPTIONS_DIR"] : nothing
    EXCEPTIONS_FILENAME === nothing ? EXCEPTIONS_FILENAME = readDefaults()["EXCEPTIONS_FILENAME"] : nothing


    gridDir, repDir, expocodeDir = load_GOSHIP_Directories(GOSHIP_DIR)
    llGrid, prGrid, sectionMask = loadSectionInfo(sectionName,gridDir,MASK_MATFILE)

    expocodeInfo = listSectionExpocodes(sectionName,expocodeDir)
    expocodes = expocodeInfo[!,"GLODAP Expocode"]

    outputArray = fill(NaN,size(sectionMask,1),size(sectionMask,2),length(expocodes))

    goodExpocodeIdx = fill(1,0)
    for expocode in enumerate(expocodes)
        try
            variable = loadGLODAPVariable(GLODAP_DIR,variableName,GLODAP_expocode=string(expocode[2]))
            goodExpocodeIdx = vcat(goodExpocodeIdx,expocode[1])
        catch
        end
    end
    goodExpocodeIdx = convert(Vector{Int64},goodExpocodeIdx)
    expocodes = expocodes[goodExpocodeIdx,:]

    for expocode in enumerate(expocodes)
        println("Starting on cruise " * expocode[2])

        isAnException = testExpocodeException(expocode=expocode[2]
                            ,variableName=variableName
                            ,EXCEPTIONS_FILENAME=EXCEPTIONS_FILENAME
                            ,EXCEPTIONS_DIR=EXCEPTIONS_DIR)

        if !isAnException
            griddingVars = load_glodap_vars(["G2longitude","G2latitude"],GLODAP_DIR,GLODAP_expocode=expocode[2])
        else
            println("Cruise " * expocode[2]* " is an exception: loading data")
            (_, griddingVars) = loadExceptionData(expocode=expocode[2]
            ,variableName=variableName,EXCEPTIONS_DIR=EXCEPTIONS_DIR)
        end

        lon = griddingVars["G2longitude"]
        lat = griddingVars["G2latitude"]


        _, _, sectionMask = loadSectionInfo(sectionName,gridDir,MASK_MATFILE)
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

function checkSectionPartialCruises(;GLODAP_DIR::Union{String,Nothing}=nothing
                                    ,GOSHIP_DIR::Union{String,Nothing}=nothing
                                    ,MASK_MATFILE::Union{String,Nothing}=nothing
                                    ,sectionName::String
                                    ,horzCoordinate::String)
    # Check a section for partial cruises.
    if horzCoordinate != "longitude" && horzCoordinate != "latitude"
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GOSHIP_DIR === nothing ? GOSHIP_DIR = readDefaults()["GOSHIP_DIR"] : nothing
    MASK_MATFILE === nothing ? MASK_MATFILE = readDefaults()["MASK_MATFILE"] : nothing

    gridDir, repDir, expocodeDir = load_GOSHIP_Directories(GOSHIP_DIR)
    llGrid, prGrid, sectionMask = loadSectionInfo(sectionName,gridDir,MASK_MATFILE)

    expocodeInfo = listSectionExpocodes(sectionName,expocodeDir)
    expocodes = expocodeInfo[!,"GLODAP Expocode"]


    goodExpocodeIdx = fill(1,0)
    for expocode in enumerate(expocodes)
        try
            variable = loadGLODAPVariable(GLODAP_DIR,"G2longitude",GLODAP_expocode=string(expocode[2]))
            goodExpocodeIdx = vcat(goodExpocodeIdx,expocode[1])
        catch
        end
    end
    goodExpocodeIdx = convert(Vector{Int64},goodExpocodeIdx)
    expocodes = expocodes[goodExpocodeIdx,:]

    truthArray = fill(false,length(expocodes))

    for expocode in enumerate(expocodes)
        println("Checking for partial occupation during cruise " * expocode[2])

        isAnException = testExpocodeException(expocode=expocode[2]
                            ,variableName="G2longitude"
                            ,EXCEPTIONS_FILENAME=EXCEPTIONS_FILENAME
                            ,EXCEPTIONS_DIR=EXCEPTIONS_DIR)
        if !isAnException
            griddingVars = load_glodap_vars(["G2longitude","G2latitude"],GLODAP_DIR,GLODAP_expocode=expocode[2])
        else
            println("Cruise " * expocode[2]* " is an exception: loading data")
            (_, griddingVars) = loadExceptionData(expocode=expocode[2]
            ,variableName=variableName,EXCEPTIONS_DIR=EXCEPTIONS_DIR)
        end

        lon = griddingVars["G2longitude"]
        lat = griddingVars["G2latitude"]


        _, _, sectionMask = loadSectionInfo(sectionName,gridDir,MASK_MATFILE)
        # Need to overwrite sectionMask on each iteration or we will have problems
        isPartialCruise = checkPartialCruise(llGrid;horzCoordinate=horzCoordinate
                                            ,obsLat=lat,obsLon=lon)
        if isPartialCruise
            truthArray[expocode[1]] = true
        end

    end

    return truthArray
end
