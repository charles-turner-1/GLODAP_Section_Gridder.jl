"""
    fit_lengths(
        vars::DataFrame,
        data_residual::Vector{Float64},
        pr_grid::Vector{Real}, 
        search_z_func::Function
    )::Tuple{Vector{Float64}, Vector{Float64}}

Fits horizontal and vertical correlation lengths for a given dataset using residuals and pressure grid.

# Arguments
- `vars::DataFrame`: A DataFrame containing longitude, latitude, and pressure.
- `data_residual::Vector{Float64}`: The residuals of the data to be used for fitting.
- `pr_grid::Vector{Real}`: The pressure grid over which the correlation lengths are computed.
- `search_z_func::Function`: A function used to search for the vertical correlation length.

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: A tuple containing:
  - `lenx::Vector{Float64}`: The horizontal correlation lengths.
  - `lenz::Vector{Float64}`: The vertical correlation lengths.

# Workflow
1. Computes horizontal correlation lengths using the `fithorzlen` function from DIVAnd.
2. Computes vertical correlation lengths using the `fitvertlen` function from DIVAnd.

# Notes
- The horizontal correlation lengths are computed based on longitude, latitude, and pressure.
- The vertical correlation lengths are constrained to a specific range using a limiting function.
- The `search_z_func` parameter customizes the search behavior for vertical correlation lengths.
- It is assumed that only residual data, not the original data, is passed to this function. 
    If original data is used, correlation length fitting may fail.
- Applies a limiting function to constrain vertical correlation lengths between 10 and 1000.
"""
function fit_lengths(
    vars::DataFrame,
    data_residual::Vector{Float64},
    pr_grid::Vector{Real}, 
    search_z_func::Function
)::Tuple{Vector{Float64}, Vector{Float64}}
    @info "Computing correlation lengths: horizontal"
    _x = (vars[!, "G2longitude"], vars[!, "G2latitude"], vars[!, "G2pressure"])

    lenx, _dbinfo = fithorzlen(
        _x, 
        data_residual, 
        pr_grid; 
        searchz=search_z_func
    )

    @info "Computing correlation lengths: vertical"
    lenz, _dbinfo = fitvertlen(
        _x, 
        data_residual, 
        pr_grid;
        searchz=search_z_func, 
        limitfun= (z, len) -> max(min(len, 1000), 10)
    )

    return lenx, lenz
end


"""
Try to fit with a multiplier value. If it fails, return this function, with the
multiplier incremented by one. In theory, this should keep trying to fit with an
increasing search window until it succeeds.
"""
function fithorzlen_increasing_search(_x, data_residual, pr_grid; search_z_func, multiplier =1 )
    
    try
        return fithorzlen(_x, data_residual, pr_grid; searchz=multiplier*search_z_func)
    catch
        @warn "fithorzlen failed with multiplier $multiplier; retrying with multiplier $(multiplier + 1)"
        return fithorzlen_increasing_search(
            _x,
            data_residual,
            pr_grid;
            search_z_func=search_z_func,
            multiplier = multiplier + 1
        )
    end
end