function readDefaults()
    # This function will look at the defaults in defaults.toml and save them. It
    # will also print them unless told not to.
    defaultsDict = TOML.parsefile("./defaults.toml")
    return defaultsDict
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

function listSectionExpocodes(sectionName::String,convDir::String)
    # Lists out all the expocodes of cruises occupying a given section
    expocodes = joinpath(convDir,sectionName) * ".csv"
    expocodes = CSV.read(expocodes,DataFrame)
    return expocodes
end


function listAvailableGLODAPVariables(GLODAP_DIR::Union{String,Nothing}=nothing,
                                      GLODAP_DATAFILE::Union{String,Nothing}=nothing)
    # Lists all variables contained in GLODAP

    GLODAP_DIR === nothing ? GLODAP_DIR = readDefaults()["GLODAP_DIR"] : nothing
    GLODAP_DATAFILE === nothing ? GLODAP_DATAFILE = readDefaults()["GLODAP_FILENAME"] : nothing

    GLODAP_DATAFILE = joinpath(GLODAP_DIR,GLODAP_DATAFILE)
    GLODAP_Data = MatFile(GLODAP_DATAFILE)
    println(" ");println("Available GLODAP variables:")
    for variable in variable_names(GLODAP_Data)
        println(variable)
    end
end