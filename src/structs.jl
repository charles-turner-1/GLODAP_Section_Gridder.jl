"""
    ScaleFactors

Container for the vertical and horizontal inverse-distance scale grids used by
DIVAnd-style interpolation steps.
"""
struct ScaleFactors
    vert::Matrix{Float64}
    horz::Matrix{Float64}
end

"""
    GriddedCruise

Container for the inputs, intermediate residuals, and final gridded output returned by `grid_cruise`.
"""
struct GriddedCruise
    expocode::String
    section_name::String
    gridding::String
    varname::String
    horz_coordinate::String
    obs_ll::Vector{Real}
    obs_pr::Vector{Real}
    residuals::Vector{Float64}
    horz_grid::Vector{Real}
    vert_grid::Vector{Real}
    mask::Matrix{Bool}
    gridded_data::Matrix{Float64}
end
