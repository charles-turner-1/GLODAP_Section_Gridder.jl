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
    # Load exception data from one of our exception datafiles
    exceptionDir, _  = splitdir(maskMatfile)
    exceptionData = joinpath(exceptionDir,"Exceptions/ExceptionData",expocode*".mat")
    SectionFile = MatFile(exceptionData)

    variable = get_variable(SectionFile,variableName)

    gridding_varNames = ["G2longitude","G2latitude","G2pressure","G2gamma","G2station"]
    griddingVars= Dict{String,Vector{Float64}}()
    for variableName in gridding_varNames
        griddingVars[variableName]  = get_variable(SectionFile,variableName)
    end

    return (variable, griddingVars)

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

function checkVariableExceptions(;expocode::Union{String,String15},variableName::String
                               ,maskMatfile::String
                               ,variable::Union{Vector{Float64},Nothing} = nothing
                               ,station::Union{Vector{Float64},Nothing} = nothing
                               ,pressure::Union{Vector{Float64},Nothing} = nothing
                               ,GLODAP_DIR::String="/home/ct/MATLAB/GLODAP")
    # Check if there are exception data for the variable and cruise we are looking 
    # at
    if variable === nothing
        vars = loadGLODAPVariables([variableName,"G2station","G2pressure"]
        ,GLODAP_DIR,GLODAP_expocode=expocode)
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

function checkHorzLenFactor(;expocode::Union{String,String15},variableName::String
                               ,maskMatfile::String, griddingType::String
                               ,GLODAP_DIR::String="/Users/ct6g18/MATLAB/GLODAP")
    # Check whether we have a manual horizonal correlation length exception. If 
    # so, apply it in the calculations.
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
