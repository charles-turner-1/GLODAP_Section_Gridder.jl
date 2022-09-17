function writeOutputField(;variable::Array{Float64,3},sectionName::String
                         ,variableName::String
                         ,GLODAP_DIR::String="/home/ct/MATLAB/GLODAP"
                         ,horzCoordinate::String
                         ,outputPath::String="/home/ct/Julia/RawOutputFields"
                         ,GLODAP_FILENAME::String="GLODAPv2.2021_Merged_Master_File.mat")
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
                         ,GLODAP_DIR::String="/home/ct/MATLAB/GLODAP"
                         ,horzCoordinate::String
                         ,outputPath::String="/home/ct/Julia/RawOutputFields"
                         ,GLODAP_FILENAME::String="GLODAPv2.2021_Merged_Master_File.mat")
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
