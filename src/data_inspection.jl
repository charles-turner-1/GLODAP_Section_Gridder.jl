function readDefaults()
    # This function will look at the defaults in defaults.toml and save them. It
    # will also print them unless told not to.
    return TOML.parsefile(joinpath(root,"defaults.toml"))
end

function changeDefaults()
    @warn("This function will overwrite your defaults.toml file. If you want to
    restore it to the default settings, use restoreDefaults().")
    defaults = TOML.parsefile(joinpath(root,"defaults.toml"))
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
    open(joinpath(root,"defaults.toml"), "w") do io
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

# function listAvailableMasks(MASK_MATFILE::Union{String,Nothing}=nothing)
function list_masks(MASK_MATFILE::Union{String,Nothing}=nothing)
    # Reads the section mask file and returns the available masks.
    MASK_MATFILE === nothing ? MASK_MATFILE = readDefaults()["MASK_MATFILE"] : nothing
    SectionMaskFile = MatFile(MASK_MATFILE) 
    maskDict = jdict(get_mvariable(SectionMaskFile,"maskStruct"))
    println("\nAvailable WOCE Section Masks:")
    for key in maskDict
        println(key.first)
    end
end

# function listSectionExpocodes(sectionName::String ,expocodeDir::Union{String,Nothing}=nothing)
function list_section_expocodes(
    section_name::AbstractString,
    expocode_dir::Union{AbstractString,Nothing}=nothing,
)::DataFrame
    # Lists out all the expocodes of cruises occupying a given section
    if expocode_dir === nothing
        expocodes = joinpath(root,"data/SectionExpocodes","$sectionName.csv")
    else # Allow manual specification of expocodeDir so user can specify something weird if they want
        expocodes = joinpath(expocode_dir,"$sectionName.csv")
    end
    return CSV.read(expocodes,DataFrame)
end


# function listAvailableGLODAPVariables(GLODAP_DIR::Union{String,Nothing}=nothing, GLODAP_DATAFILE::Union{String,Nothing}=nothing)
function list_glodap_vars(
    glodap_dir::Union{String,Nothing}=nothing,
    glodap_datafile::Union{String,Nothing}=nothing
)::Nothing
    # Lists all variables contained in GLODAP

    glodap_dir === nothing ? glodap_dir = readDefaults()["GLODAP_DIR"] : nothing
    glodap_datafile === nothing ? glodap_datafile = readDefaults()["GLODAP_FILENAME"] : nothing

    glodap_datafile = joinpath(glodap_dir,glodap_datafile)
    GLODAP_Data = MatFile(GLODAP_DATAFILE) 
    # Need to excise this & replace with reading the CSV headers or whatever
    println("\nAvailable GLODAP variables:")
    for variable in variable_names(GLODAP_Data)
        println("\t--> $variable")
    end
end
