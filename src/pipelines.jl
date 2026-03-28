using DIVAnd

include("./structs.jl")
include("./data_loading.jl")
include("./correlation_lengths.jl")
include("./simple_functionality.jl")

function gridCruisePipeline(;GLODAP_DIR::Union{String,Nothing}=nothing
                            ,GOSHIP_DIR::Union{String,Nothing}=nothing
                            ,MASK_MATFILE::Union{String,Nothing}=nothing
                            ,EXCEPTIONS_DIR::Union{String,Nothing}=nothing
                            ,EXCEPTIONS_FILENAME::Union{String,Nothing}=nothing
                            ,VARIABLE_EXCEPTIONS::Union{String,Nothing}=nothing
                            ,HORZLEN_EXCEPTIONS::Union{String,Nothing}=nothing
                            ,sectionName::String
                            ,horzCoordinate::String
                            ,expocode::String
                            ,variableName::String
                            ,gridding::String="isobaric"
                            ,meanValue::String="scalar"
                            ,epsilonVal::Float64=0.1
                            ,plotResults::Bool=false
                            ,autoTruncateMask::Bool=false
                            ,crossValidate::Bool=false
                            ,crossValidationNum::Int64=5)
    # Calls all the functions necessary to grid a single cruise.
    check_gridding_vars(horzCoordinate,gridding,meanValue)

    # Now load all the defaults if needs be
    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GOSHIP_DIR === nothing ? GOSHIP_DIR = readDefaults()["GOSHIP_DIR"] : nothing
    MASK_MATFILE === nothing ? MASK_MATFILE = readDefaults()["MASK_MATFILE"] : nothing
    EXCEPTIONS_FILENAME === nothing ? EXCEPTIONS_FILENAME = readDefaults()["EXCEPTIONS_FILENAME"] : nothing
    VARIABLE_EXCEPTIONS === nothing ? VARIABLE_EXCEPTIONS = readDefaults()["VARIABLE_EXCEPTIONS"] : nothing
    HORZLEN_EXCEPTIONS === nothing ? HORZLEN_EXCEPTIONS = readDefaults()["HORZLEN_EXCEPTIONS"] : nothing


    meanValue == "climatology" ? bgField = read_bgfield(
    sectionName=sectionName,variableName=variableName) : bgField = nothing

    if crossValidate == false
        crossValidationNum = 0
        bestFactorLen = 1
        bestFactorε = 1
    end

    llGrid, prGrid, sectionMask = loadSectionInfo(sectionName,GOSHIP_DIR,MASK_MATFILE)

    isAnException = testExpocodeException(expocode=expocode,variableName=variableName
         ,EXCEPTIONS_FILENAME=EXCEPTIONS_FILENAME,EXCEPTIONS_DIR=EXCEPTIONS_DIR)

    if !isAnException
        variable = loadGLODAPVariable(variableName,GLODAP_DIR,GLODAP_expocode=expocode)
        griddingVars = load_glodap_vars(["G2longitude","G2latitude","G2pressure","G2gamma","G2station"],GLODAP_DIR,GLODAP_expocode=expocode)

        varNeedsFlags = has_dataflags(variableName)
        if varNeedsFlags == true
            println("Removing flagged data")
            varFlagsName::String = variableName * "f"
            variableFlags = loadGLODAPVariable(varFlagsName,GLODAP_DIR,GLODAP_expocode=expocode)
            variable = rm_flagged_data(variable,variableFlags)
        end

    else
        println("Cruise "* expocode * " is an exception: loading exception data")
        (variable, griddingVars) = loadExceptionData(expocode=expocode
            ,variableName=variableName,EXCEPTIONS_DIR=EXCEPTIONS_DIR)
    end

    if variableName == "G2tco2"
        tco2Adj = adjust_tco2(expocode=expocode)
        variable .+= tco2Adj
        tco2Adj == 0 ? println("No DIC adjustment necessary") : println("Adjusting DIC by " * string(tco2Adj) * "μmol/kg")
    end

    lon = griddingVars["G2longitude"]
    lat = griddingVars["G2latitude"]
    pressure = griddingVars["G2pressure"]
    sigma = griddingVars["G2gamma"]
    station = griddingVars["G2station"]

    goodIdx = checkVariableExceptions(expocode=expocode,variableName=variableName
    ,EXCEPTIONS_FILENAME=VARIABLE_EXCEPTIONS,GLODAP_DIR=GLODAP_DIR,variable=variable
    ,station=station,pressure=pressure,EXCEPTIONS_DIR=EXCEPTIONS_DIR)

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

    horzDist = grid_horz_dist(lat,lon,llGrid)
    vertDist = grid_vert_dist(prGrid)

    scales = calc_scale_factors(vertDist,horzDist)
    scaleVert, scaleHorz = scales.vert, scales.horz

    lenxFactor = checkHorzLenFactor(expocode=expocode,variableName=variableName
                   ,griddingType=gridding,HORZLEN_EXCEPTIONS=HORZLEN_EXCEPTIONS)

    len = calcCorrLengths(variable,obsLat=lat,obsLon=lon,obsPres=pressure
    ,presGrid=prGrid,pressureStepNumber=10,verticalSearchRange=100,lenxFactor=lenxFactor)

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
        lenxPrescribed = checkHorzLenFactor(expocode=expocode,variableName=variableName
                   ,griddingType=gridding,HORZLEN_EXCEPTIONS=HORZLEN_EXCEPTIONS)

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

function gridSectionPipeline(;sectionName::String
                             ,variableName::String
                             ,horzCoordinate::String
                             ,gridding::String="isobaric"
                             ,meanValue::String="scalar"
                             ,epsilonVal::Float64=0.1
                             ,autoTruncateMask::Bool=false
                             ,GLODAP_DIR::Union{String,Nothing}=nothing
                             ,GOSHIP_DIR::Union{String,Nothing}=nothing
                             ,MASK_MATFILE::Union{String,Nothing}=nothing
                             ,EXCEPTIONS_DIR::Union{String,Nothing}=nothing
                             ,EXCEPTIONS_FILENAME::Union{String,Nothing}=nothing
                             ,EXPOCODE_DIR::Union{Nothing,String}=nothing
                             ,VARIABLE_EXCEPTIONS::Union{String,Nothing}=nothing
                             ,HORZLEN_EXCEPTIONS::Union{String,Nothing}=nothing)
    # Calls all the functions necessary to grid all the cruises from an entire 
    # section
    check_gridding_vars(horzCoordinate,gridding,meanValue)

    meanValue == "climatology" ? bgField = read_bgfield(
    sectionName=sectionName,variableName=variableName) : bgField = nothing

    # Now load all the defaults if needs be
    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GOSHIP_DIR === nothing ? GOSHIP_DIR = readDefaults()["GOSHIP_DIR"] : nothing
    MASK_MATFILE === nothing ? MASK_MATFILE = readDefaults()["MASK_MATFILE"] : nothing
    #EXCEPTIONS_DIR === nothing ? EXCEPTIONS_DIR = readDefaults()["EXCEPTIONS_DIR"] : nothing
    EXCEPTIONS_FILENAME === nothing ? EXCEPTIONS_FILENAME = readDefaults()["EXCEPTIONS_FILENAME"] : nothing
    VARIABLE_EXCEPTIONS === nothing ? VARIABLE_EXCEPTIONS = readDefaults()["VARIABLE_EXCEPTIONS"] : nothing
    HORZLEN_EXCEPTIONS === nothing ? HORZLEN_EXCEPTIONS = readDefaults()["HORZLEN_EXCEPTIONS"] : nothing

    if EXPOCODE_DIR === nothing
        gridDir, repDir, EXPOCODE_DIR = load_GOSHIP_Directories(GOSHIP_DIR)
    else
        gridDir, repDir, _ = load_GOSHIP_Directories(GOSHIP_DIR)
    end

    llGrid, prGrid, sectionMask = loadSectionInfo(sectionName,GOSHIP_DIR,MASK_MATFILE)

    expocodeInfo = listSectionExpocodes(sectionName,EXPOCODE_DIR)
    expocodes = expocodeInfo[!,"GLODAP Expocode"]

    outputArray = fill(NaN,size(sectionMask,1),size(sectionMask,2),length(expocodes))

    goodExpocodeIdx = fill(1,0)
    for expocode in enumerate(expocodes)
        try
    isAnException = testExpocodeException(expocode=expocode[2],variableName=variableName
         ,EXCEPTIONS_FILENAME=EXCEPTIONS_FILENAME,EXCEPTIONS_DIR=EXCEPTIONS_DIR)
            if !isAnException
                variable = loadGLODAPVariable(variableName,GLODAP_DIR,GLODAP_expocode=expocode[2])
            else
                (variable, griddingVars) = loadExceptionData(expocode=expocode[2]
            ,variableName=variableName,EXCEPTIONS_DIR=EXCEPTIONS_DIR)
            end
            goodExpocodeIdx = vcat(goodExpocodeIdx,expocode[1])
        catch
        end
    end
    goodExpocodeIdx = convert(Vector{Int64},goodExpocodeIdx)
    expocodes = expocodes[goodExpocodeIdx,:]

    for expocode in enumerate(expocodes)
        println("Starting on cruise " * expocode[2])

        isAnException = testExpocodeException(expocode=expocode[2]
            ,variableName=variableName,EXCEPTIONS_FILENAME=EXCEPTIONS_FILENAME
            ,EXCEPTIONS_DIR=EXCEPTIONS_DIR)

        if !isAnException
            variable = loadGLODAPVariable(variableName,GLODAP_DIR,GLODAP_expocode=expocode[2])
            griddingVars = load_glodap_vars(["G2longitude","G2latitude","G2pressure"
            ,"G2gamma","G2station"],GLODAP_DIR,GLODAP_expocode=expocode[2])

            varNeedsFlags = has_dataflags(variableName)
            if varNeedsFlags == true
                println("Removing flagged data")
                varFlagsName::String = variableName * "f"
                variableFlags = loadGLODAPVariable(varFlagsName,GLODAP_DIR,GLODAP_expocode=expocode[2])
                variable = rm_flagged_data(variable,variableFlags)
            end
            if variableName == "G2tco2"
                tco2Adj = adjust_tco2(expocode=expocode[2])
                tco2Adj == 0 ? println("No DIC adjustment necessary") : println("Adjusting DIC by " * string(tco2Adj) * "μmol/kg")
                variable .+= tco2Adj
            end

        else
            println("Cruise " * expocode[2]* " is an exception: loading data")
            (variable, griddingVars) = loadExceptionData(expocode=expocode[2]
            ,variableName=variableName,EXCEPTIONS_DIR=EXCEPTIONS_DIR)
        end

        lon = griddingVars["G2longitude"]
        lat = griddingVars["G2latitude"]
        pressure = griddingVars["G2pressure"]
        sigma = griddingVars["G2gamma"]
        station = griddingVars["G2station"]

        goodIdx = checkVariableExceptions(expocode=expocode[2],variableName=variableName
    ,EXCEPTIONS_FILENAME=VARIABLE_EXCEPTIONS,GLODAP_DIR=GLODAP_DIR,variable=variable
    ,station=station,pressure=pressure,EXCEPTIONS_DIR=EXCEPTIONS_DIR)

        lon = lon[goodIdx]; lat = lat[goodIdx]; pressure = pressure[goodIdx];
        sigma = sigma[goodIdx]; station = station[goodIdx]; variable = variable[goodIdx]

        if horzCoordinate == "longitude"
            XDIR = lon
        elseif horzCoordinate == "latitude"
            XDIR = lat
        end

        if autoTruncateMask
            _, _, sectionMask = loadSectionInfo(sectionName,GOSHIP_DIR,MASK_MATFILE)
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

        horzDist = grid_horz_dist(lat,lon,llGrid)
        vertDist = grid_vert_dist(prGrid)

        scales = calc_scale_factors(vertDist,horzDist)
        scaleVert, scaleHorz = scales.vert, scales.horz

        lenxFactor = checkHorzLenFactor(expocode=expocode[2],variableName=variableName
                   ,griddingType=gridding,HORZLEN_EXCEPTIONS=HORZLEN_EXCEPTIONS)

        len = calcCorrLengths(variable,obsLat=lat,obsLon=lon,obsPres=pressure
        ,presGrid=prGrid,pressureStepNumber=10,verticalSearchRange=100,lenxFactor=lenxFactor)

        if gridding == "isobaric"
            griddedVarEasyPipeline = easyDIVAGrid(variable=variable,vertVar=pressure,latLon=XDIR
            ,horzCoordinate=horzCoordinate,vertGrid=prGrid,horzGrid=llGrid,horzScale=scaleHorz
            ,vertScale=scaleVert,mask=sectionMask,Epsilon=epsilonVal,horzCorrLength=len[2]
            ,vertCorrLength=len[1],meanValue=meanValue,bgField=bgField)
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

function gridExceptionPipeline(;GLODAP_DIR::Union{String,Nothing}=nothing
                            ,GOSHIP_DIR::Union{String,Nothing}=nothing
                            ,MASK_MATFILE::Union{String,Nothing}=nothing
                            ,EXCEPTIONS_DIR::Union{String,Nothing}=nothing
                            ,horzLenFactor::Float64=1.0
                            ,sectionName::String
                            ,horzCoordinate::String
                            ,expocode::String
                            ,variableName::String
                            ,gridding::String="isobaric"
                            ,meanValue::String="scalar"
                            ,epsilonVal::Float64=0.1
                            ,plotResults::Bool=false
							,autoTruncateMask::Bool=false)

    # Calls all the functions necessary to grid a single cruise.
    check_gridding_vars(horzCoordinate,gridding,meanValue)

    # Now load all the defaults if needs be
    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GOSHIP_DIR === nothing ? GOSHIP_DIR = readDefaults()["GOSHIP_DIR"] : nothing
    MASK_MATFILE === nothing ? MASK_MATFILE = readDefaults()["MASK_MATFILE"] : nothing

    meanValue == "climatology" ? bgField = read_bgfield(
    sectionName=sectionName,variableName=variableName) : bgField = nothing


    llGrid, prGrid, sectionMask = loadSectionInfo(sectionName,GOSHIP_DIR,MASK_MATFILE)

    (variable, griddingVars) = loadExceptionData(expocode=expocode
            ,variableName=variableName,EXCEPTIONS_DIR=EXCEPTIONS_DIR)

    lon = griddingVars["G2longitude"]
    lat = griddingVars["G2latitude"]
    pressure = griddingVars["G2pressure"]
    sigma = griddingVars["G2gamma"]
    station = griddingVars["G2station"]


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

    horzDist = grid_horz_dist(lat,lon,llGrid)
    vertDist = grid_vert_dist(prGrid)

    scales = calc_scale_factors(vertDist,horzDist)
    scaleVert, scaleHorz = scales.vert, scales.horz

    lenxFactor = horzLenFactor

    len = calcCorrLengths(variable,obsLat=lat,obsLon=lon,obsPres=pressure
    ,presGrid=prGrid,pressureStepNumber=10,verticalSearchRange=100,lenxFactor=lenxFactor)

    if gridding == "isobaric"
        griddedVarEasyPipeline = easyDIVAGrid(variable=variable,vertVar=pressure
        ,latLon=XDIR,vertGrid=prGrid,horzGrid=llGrid,horzScale=scaleHorz
        ,vertScale=scaleVert,mask=sectionMask,Epsilon=epsilonVal
        ,horzCorrLength=len[2],vertCorrLength=len[1]
        ,horzCoordinate=horzCoordinate,meanValue=meanValue,bgField=bgField)
    elseif gridding == "isopycnic"
        lenxPrescribed = checkHorzLenFactor(expocode=expocode,variableName=variableName
                   ,griddingType=gridding,HORZLEN_EXCEPTIONS=HORZLEN_EXCEPTIONS)

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


"""
grid_cruise(expocode::String, section_name::String, varname::String, gridding::String="isobaric")::GriddedCruise

Grids cruise data for a specified section, variable, and gridding method using the DIVAnd interpolation framework.

# Arguments
- `expocode::String`: The unique identifier for the cruise.
- `section_name::String`: The name of the section to be gridded.
- `varname::String`: The variable name to be gridded (e.g., temperature, salinity).
- `gridding::String="isobaric"`: The gridding method to use. Defaults to "isobaric".

# Returns
- `GriddedCruise`: A structured object containing the gridded data and associated metadata.

# Workflow
1. Logs the gridding process and validates the gridding method.
2. Loads section information, including coordinates and masks.
3. Loads GLODAP data for the specified cruise and variables.
4. Removes the scalar mean from the variable data.
5. Performs an initial DIVAnd fitting with long correlation lengths.
6. Computes residuals and fits correlation lengths based on the data.
7. Performs a second DIVAnd fitting using computed correlation lengths.
8. Combines the mean, first fit, and second fit to produce the final gridded data.
9. Returns a `GriddedCruise` object containing the gridded data, residuals, and metadata.

# Notes
- The function assumes the availability of GLODAP data at the specified file path.
- The correlation lengths for the second DIVAnd fitting are computed dynamically from the data.
- The `pmn` parameter for the second fitting is currently hardcoded and may require adjustment.

# Example
```julia
gridded_cruise = grid_cruise("33RO20230101", "A05", "G2theta")
"""
function grid_cruise(
    expocode::String,
    section_name::String,
    varname::String,
    gridding::String="isobaric",
)::GriddedCruise

    @info "Gridding cruise $expocode for section $section_name with variable $varname using $gridding gridding method."
    check_gridding_vars(gridding,"scalar")

    @info "Loading section information for $section_name"
    section_info = load_section_info(section_name)

    ll_grid, pr_grid, section_mask, horz_coordinate = section_info.ll_grid, section_info.pr_grid, section_info.mask, section_info.horz_coord

    G2horz_coord = "G2$horz_coordinate"

    secondary_horz_coord = horz_coordinate == "longitude" ? "G2latitude" : "G2longitude"

    @info "Loading GLODAP data for $expocode, variables $varname, $G2horz_coord, G2pressure"
    vars = load_glodap_vars(
        [varname,"G2pressure",G2horz_coord, secondary_horz_coord],
        expocode,
        "/Users/u1166368/GLODAP/GLODAPv2.2023_Merged_Master_File.csv"
    )
    
    var = vars[!, varname]
    pr_vals = vars[!, "G2pressure"]
    ll_vals = vars[!, G2horz_coord]

    var_mean, var_anom = remove_scalar_mean(var)

    pmn = ones(size(section_mask)), ones(size(section_mask))

    @info "Performing initial (long correlation) DIVAnd fitting"
    fi, s = DIVAndrun(section_mask, pmn, (ndgrid(ll_grid, pr_grid)), (ll_vals, pr_vals), 
                      var_anom, (200,200), 0.1)

                      
    data_residual = DIVAnd_residual(s, fi)

    
    lenx, lenz = fit_lengths(vars, data_residual, pr_grid, search_z_func)
    
    @info "Performing second (computed correlation) DIVAnd fitting on residuals"

    lenx = repeat(lenx, 1, length(ll_grid))'
    lenz = repeat(lenz, 1, length(ll_grid))'

    #= 
    *** NEED TO CORRECTLY FIT THIS FROM THE DATA ***
    =#
    pmn = ones(size(mask)), 100 * ones(size(mask))

    fi2, s2 = DIVAndrun(mask, pmn, (ndgrid(ll_grid, pr_grid)), (ll_vals, pr_vals), 
    data_residual, (lenz,lenx), 0.1)
    
    gridded_data = var_mean .+ fi .+ fi2

    data_residual2 = DIVAnd_residual(s2, fi2)

    residuals = data_residual .+ data_residual2
    
    return GriddedCruise(
        section_name,
        expocode,
        gridding,
        varname,
        horz_coordinate,
        ll_vals,
        pr_vals,
        residuals,
        ll_grid,
        pr_grid,
        mask,
        gridded_data'
        )
        
    end
    