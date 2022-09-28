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

function centralDiff(v::AbstractVector)
    # Very simple central difference funciton
    dvF  = diff(v)
    dvB  = reverse(-diff(reverse(v)))
    dx = Vector{AbstractFloat}(undef,length(dvF)+1)
    dx[1] = dvF[1]; dx[end] = dvB[end]
    dx[2:end-1] = (dvF[1:end-1] + dvB[2:end]) / 2
    return dx
end

function splitMeanAnom(variable::Vector{Float64})
    # Remove the mean value from a vector.
    varAnom =  variable .- mean(filter(!isnan,variable))
    varMean = mean(filter(!isnan,variable))
    return varMean, varAnom
end


# Not sure this should be in here
function splitMeanAnom(;obsVariable::Vector{Float64},obsPres::Vector{Float64} #kwargs might not be necessary since this is at the bottom of the call stack
                      ,obsLatLon::Vector{Float64},bgField::Matrix{Float64}
                      ,vertGrid::Vector{Float64},horzGrid::Vector{Float64})
    # Remove the mean value from a vector, removing the horizontal mean value at 
    # each depth
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

function findNonNaNIndices(variable::Vector)
    # Very simple function to find indices where vector isn't NaN
    goodIdx = findall(convert(Vector{Bool},fill(1,size(variable)) - isnan.(variable)))
    return goodIdx
end

function findNonNaNIndices(variable::Vector,sigma::Vector)
    # Find all indices where our vector isn't NaN, and the sigma vector isn't NaN
    varGoodIdx = findall(convert(Vector{Bool},fill(1,size(variable)) - isnan.(variable)))
    sigGoodIdx = findall(convert(Vector{Bool},fill(1,size(sigma)) - isnan.(sigma)))
    goodIdx = intersect(varGoodIdx,sigGoodIdx)
    return goodIdx
end

function calcGLODAP_Date(glodapYear::Vector{Float64},glodapMonth::Vector{Float64}
                        ,glodapDay::Vector{Float64})
    # Calculate a mean date for a cruise
    glodapMonth .-= 1
    glodapDate = glodapYear + (glodapMonth ./ 12) + (glodapDay ./ 30)
    return glodapDate
end

function splitMeanAnom(obsVariable::Vector{Float64},obsPres::Vector{Float64}
                      ,bgField::Vector{Float64},vertGrid::Vector{Float64})
    # Remove the mean value from observations at each index of our vertical grid
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

function matchLonConventions(lonGrid::Vector{Float64}, lonVal::Float64)
    # Match a single longitude value up with conventions
    if maximum(lonGrid) > 180 && lonVal < 0
        lonVal += 360
    end
    return lonVal
end

function checkGriddingVariables(horzCoordinate::String,gridding::String,meanValue::String)
    # Checks that we are using defined variables. 
    if horzCoordinate ∉ ["longitude","latitude"]
        error("\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\"")
    end

    if gridding ∉ ["isobaric","isopycnic"]
        error("\"gridding\" must be specified to be either \"isobaric\" or \"isopycnic\"")
    end

    if meanValue ∉ ["horzMean","scalar","climatology","calculated"]
        error("\"meanValue\" must be specified to be either \"horzMean\", \"calculated\" ,\"scalar\" or \"climatology\"")
    end
end