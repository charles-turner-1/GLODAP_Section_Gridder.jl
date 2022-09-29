function readDefaults()
    # This function will look at the defaults in defaults.toml and save them. It
    # will also print them unless told not to.
    return TOML.parsefile("./defaults.toml")
end

function changeDefaults()
    @warn("This function will overwrite your defaults.toml file. If you want to
    restore it to the default settings, use restoreDefaults().")
    defaults = TOML.parsefile("./defaults.toml")
    changeFile = true
    while changeFile == true
        display(defaults)
        println("Do you wish to change default file? Y/n")
        if lowercase(readline()) == "n" 
            changeFile = false; continue
        end
        println("Which field do you wish to change?")
        fieldName = readline()
        println("What do you wish to change \"" * fieldName * "\" to?")
        fieldVal = readline()
        haskey(defaults,fieldName) ? defaults[fieldName] = fieldVal :
        @warn("Field \"" * fieldName * "\" not found")
    end
    @info("Writing new defaults")
    open("./defaults.toml", "w") do io
        TOML.print(io,defaults)
    end
    return nothing
end

function restoreDefaults()
    @warn("This function will restore your \"defaults.toml\" file to its default
    state (haha), if you have overwritten it and somehow buggered it up. Use with
    caution")
    println("Are you sure you wish to proceed? y/N")
    lowercase(readline()) == "y" ? nothing : return nothing
    cp(joinpath(root,".defaults.toml"),joinpath(root,"defaults.toml"),force=true)
    return nothing
end

function listAvailableMasks(MASK_MATFILE::Union{String,Nothing}=nothing)
    # Reads the section mask file and returns the available masks.
    MASK_MATFILE === nothing ? MASK_MATFILE = readDefaults()["MASK_MATFILE"] : nothing
    SectionMaskFile = MatFile(MASK_MATFILE) 
    maskDict = jdict(get_mvariable(SectionMaskFile,"maskStruct"))
    println("\nAvailable WOCE Section Masks:")
    for key in maskDict
        println(key.first)
    end
end

function listSectionExpocodes(sectionName::String
                             ,expocodeDir::Union{String,Nothing}=nothing)
    # Lists out all the expocodes of cruises occupying a given section
    if expocodeDir === nothing
        expocodes = joinpath(root,"data/SectionExpocodes",sectionName) * ".csv"
    else # Allow manual specification of expocodeDir so user can specify something weird if they want
        expocodes = joinpath(expocodeDir,sectionName) * ".csv"
    end
    return CSV.read(expocodes,DataFrame)
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
