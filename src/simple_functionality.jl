using Base.Threads
using Statistics
using Dates


function modulo_lon(lon_grid::AbstractVector{<:Real}, lon_vals::AbstractVector{<:Real})::AbstractVector{<:Real}
    # GO-SHIP Easy Ocean toolbox and GLODAP use different longitude conventions
    # ie. lon ⋹ [-180,180] or lon ⋹ [0,360]. This will ensure tat the match up.
    if maximum(lon_grid) > 180
        lon_vals .= mod.(lon_vals, 360) 
    end
    return lon_vals
end

function modulo_lon(lon_grid::AbstractVector{<:Real}, lon_val::Real)::Real
    # Match a single longitude value up with conventions
    if maximum(lon_grid) > 180
        return mod(lon_val, 360)
    end
    return lon_val
end

function to_mask_name(section_name::String)::String
    # Gets the mask name from the name of the section we're using.
    return replace("mask$section_name","-" => "_")
end

function central_diff(v::AbstractVector{<:Real})::AbstractVector{<:Real}
    # Very simple central difference funciton
    dv_fwds  = diff(v)
    dv_bwds  = reverse(-diff(reverse(v)))
    dx = Vector{AbstractFloat}(undef,length(v))
    dx[1] = dvF[1]; dx[end] = dvB[end]
    dx[2:end-1] = (dv_fwds[1:end-1] + dv_bwds[2:end]) / 2
    return dx
end


# function splitMeanAnom(variable::Vector{Float64})
function remove_scalar_mean(var::AbstractVector{Float64})::Tuple{Float64,AbstractVector{Float64}}
    # Remove the mean value from a vector.
    var_mean = mean(filter(!isnan,var))
    var_anom =  var .- var_mean
    return var_mean, var_anom
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
        goodIdx = non_nan_indices(varInBox); varInBox = varInBox[goodIdx]
        varAnom[idx] .= obsVariable[idx] .- bgVal
        varMean[idx] .= bgVal
    end

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
                    goodIdx = non_nan_indices(varInBox); varInBox = varInBox[goodIdx]
                    varAnom[idx] = obsVariable[idx] .- bgField[vertVal,horzVal]
                end
        end
    end
    return varAnom
end

function non_nan_indices(variable::AbstractVector{<:Real})
    # Very simple function to find indices where vector isn't NaN
    return findall(convert(Vector{Bool},fill(1,size(variable)) - isnan.(variable)))
end

function non_nan_index_pairs(variable::AbstractVector{<:Real},sigma::AbstractVector{<:Real})
    # Find all indices where our vector isn't NaN, and the sigma vector isn't NaN
    var_idx = findall(convert(Vector{Bool},fill(1,size(variable)) - isnan.(variable)))
    sigma_idx = findall(convert(Vector{Bool},fill(1,size(sigma)) - isnan.(sigma)))
    return intersect(var_idx,sigma_idx)
end

function calc_decimal_year(year::AbstractVector{<:Real},month::AbstractVector{<:Real},day::AbstractVector{<:Real})
    # Calculate a mean date for a cruise

    date = Date.(year,month,day)
    doy = dayofyear.(date)
    year_len = daysinyear.(year)

    return year .+ ((doy .-1) ./ year_len)
    
end

struct InvalidHorzCoordError <: Exception
    msg::String
end

struct InvalidGriddingError <: Exception
    msg::String
end

struct InvalidMeanValueError <: Exception
    msg::String
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
    goodIdx = non_nan_indices(varInBox); varInBox = varInBox[goodIdx]
    varAnom[idx] .= obsVariable[idx] .- bgVal
    varMean[idx] .= bgVal
end

return varMean, varAnom
end

function check_gridding_vars(horz_coord::String,gridding::String,mean_val::String)
    # Checks that we are using defined variables. 
    if !(horz_coord in Set(["longitude","latitude"]))
        throw(InvalidHorzCoordError(
            "\"horzCoordinate\" must be specified to be either \"longitude\" or \"latitude\""
            ))
    end

    if !(gridding in Set(["isobaric","isopycnic"]))
        throw(InvalidGriddingError(
            "\"gridding\" must be specified to be either \"isobaric\" or \"isopycnic\""
        ))
    end

    if !(mean_val in Set(["horzMean","scalar","climatology","calculated"]))
        throw(InvalidMeanValueError(
            "\"meanValue\" must be specified to be either \"horzMean\", \"calculated\" ,\"scalar\" or \"climatology\""
        ))
    end
end
