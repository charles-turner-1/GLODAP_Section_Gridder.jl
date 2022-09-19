function load_GOSHIP_Directories(GOSHIP_DIR::Union{String,Nothing}=nothing)
    # Give the paths of the GO_SHIP data and return the files we
    # need to load
    GOSHIP_DIR === nothing ? GOSHIP_DIR = readDefaults()["GOSHIP_DIR"] : nothing

    GRID_DIR = joinpath(GOSHIP_DIR , "go_ship_clean_ctd/gridded/" )
    REP_DIR  = joinpath(GOSHIP_DIR , "go_ship_clean_ctd/reported/")
    CONV_DIR = joinpath(root,"data/SectionExpocodes/")

    return GRID_DIR, REP_DIR, CONV_DIR
end

function loadSectionInfo(sectionName::String
                        ,GRID_INFO_DIR::Union{String,Nothing}=nothing
                        ,MASK_MATFILE::Union{String,Nothing}=nothing)
    # Loads the grid data for the section, ie. lat/lon grid, pressure grid, mask.
    # This needs extending so that a user can add their own section if they want
    # one.

    GRID_INFO_DIR === nothing ? GOSHIP_DIR = readDefaults()["GOSHIP_DIR"] : nothing
    MASK_MATFILE === nothing ? MASK_MATFILE = readDefaults()["MASK_MATFILE"] : nothing

    SectionMaskFile = MatFile(MASK_MATFILE)
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

    gridFile = MatFile(joinpath(GOSHIP_DIR,"go_ship_clean_ctd/gridded/",maskPath))
    pr_grid = vec(get_variable(gridFile,"pr_grid"))
    ll_grid = vec(get_variable(gridFile,"ll_grid"))

    if basin != 'P'
        ll_grid[ll_grid .> 180] .-=360 # So GLODAP and GO-SHIP data are using the
        # same longitudes. Probably will want to update this to be basin dependent
        # or we might run into problems in the Pacific
    end

    return ll_grid, pr_grid, mask
end

function loadGLODAPVariable(GLODAP_VariableNames::String
    ,GLODAP_expocodes::Union{Vector{String},Vector{String15}}
    ,GLODAP_DIR::Union{String,Nothing}=nothing
    ,GLODAP_FILENAME::Union{String,Nothing}=nothing)
    # Loads multiple variables from GLODAP, for a number of cruises

    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GLODAP_FILENAME === nothing ? GLODAP_FILENAME = readDefaults()["GLODAP_FILENAME"] : nothing

    GLODAP_DATAFILE = joinpath(GLODAP_DIR,GLODAP_FILENAME)
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

function loadGLODAPVariable(GLODAP_VariableName::String
    ,GLODAP_DIR::Union{String,Nothing}=nothing
    ;GLODAP_FILENAME::Union{String,Nothing}=nothing
    ,GLODAP_expocode::Union{String,String15,Nothing} = nothing)
    # Loads a single variable from GLODAP, for either a given expocode or the 
    # entirety of the GLODAP dataset.
    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GLODAP_FILENAME === nothing ? GLODAP_FILENAME = readDefaults()["GLODAP_FILENAME"] : nothing

    GLODAP_DATAFILE = joinpath(GLODAP_DIR,GLODAP_FILENAME)
    GLODAP_Data = MatFile(GLODAP_DATAFILE)

    variable = get_variable(GLODAP_Data,GLODAP_VariableName)
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
    return variable
end


function loadGLODAPvariables(GLODAP_VariableNames::Vector{String}
    ,GLODAP_DIR::Union{String,Nothing}=nothing
    ,GLODAP_expocode::Union{String,String15,Nothing}=nothing
    ,GLODAP_FILENAME::Union{String,Nothing}=nothing)
    # Loads multiple variables from GLODAP, for eiher a single cruise, or the 
    # entire GLODAP dataset

    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GLODAP_FILENAME === nothing ? GLODAP_FILENAME = readDefaults()["GLODAP_FILENAME"] : nothing
    GLODAP_DATAFILE = joinpath(GLODAP_DIR,GLODAP_FILENAME)
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


function loadGLODAPVariables(GLODAP_VariableNames::Vector{String}
    ,GLODAP_DIR::Union{String,Nothing}=nothing
    ;GLODAP_FILENAME::Union{String,Nothing}=nothing
    ,GLODAP_expocodes::Union{Vector{String},Vector{String15},Nothing}=nothing
    ,GLODAP_expocode::Union{String,String15,Nothing}=nothing)
    # We can't use multiple dispatch with keyword arguments, so I'm going to put
    # this wrapper around loadGLODAPVariables to sort the issue. This function 
    # will allow you to load either a single or multiple variables for a single 
    # cruise, multiple cruises, or the whole GLODAP dataset.

    # COME BACK AND FIX THIS AT SOME POINT BECAUSE ITS NOT VERY CLEAN
    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing

    if GLODAP_expocodes !== nothing # Structuring logic this way will load all expocodes if both are passed
        println("GLODAP_expocodes is a vector")
        variables = loadGLODAPvariables(GLODAP_VariableNames,GLODAP_DIR
                                       ,GLODAP_expocodes,GLODAP_FILENAME)
    else
        println("GLODAP_expocode is a string")
        variables = loadGLODAPvariables(GLODAP_VariableNames,GLODAP_DIR
                                       ,GLODAP_expocode,GLODAP_FILENAME)
    end
    return variables
end

function hasDataFlags(variableName::String)
    # Check whether we are gridding a varible which has flagged data
    flaggedVariables = ["G2aou","G2c13","G2c14","G2ccl4","G2cfc113","G2cfc11"
            ,"G2cfc12","G2chla","G2doc","G2fco2","G2h3","G2he3","G2neon"
            ,"G2nitrate","G2nitrite","G2o18","G2oxygen","G2phosphate"
            ,"G2phtsinsitutp","G2salinity","G2sf6","G2silicate","G2talk"
            ,"G2tco2","G2tdn","G2toc"]

    if variableName in flaggedVariables; return true 
    else; return false 
    end
end

function removeFlaggedData(variable::Vector{Float64}, variable_fFlag::Vector{Float64})
    # Remove data which has been flagged as being bad
    outputVar = copy(variable)
    outputVar[variable_fFlag .!= 2] .= NaN
    if length(filter(!isnan,outputVar)) == 0
        outputVar = copy(variable)
        goodIdx = [x==0 || x==2 for x in variable_fFlag]
        outputVar[goodIdx.==0] .= NaN
    end

    return outputVar
end

function findGLODAPtco2Adjustment(;GLODAP_AdjTable::String="AdjustmentTable.csv"
                                 ,expocode::Union{String,String15})
    # Adjust tCO2 using the GLODAP adjustment table
    
    # This will check whether the expocode we are looking at is in Jens Muller's 
    # recommended adjustments and adjust up here.
    JensList = ["316N19941201","316N19950124","316N19950310","316N19950423",
    "316N19950611","316N19950715","316N19950829","316N19951111","316N19951202"]

    if expocode in JensList
        println("Expocode in Jens Muller's recommended adjustment list, returning 1.7")
        println("This is temporary behaviour until the next GLODAP update.")
        return 1.7
    end

    GLODAP_AdjTable = joinpath(root,"data",GLODAP_AdjTable)

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

function expocodeFromG2cruise(;GLODAP_DIR::Union{String,Nothing}=nothing
                              ,GLODAP_FILENAME::Union{String,Nothing}=nothing
                              ,G2cruise::Union{Float64,Vector{Float64}})
    # Get an expocode back from a value of G2cruise

    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GLODAP_FILENAME === nothing ? GLODAP_FILENAME = readDefaults()["GLODAP_FILENAME"] : nothing

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