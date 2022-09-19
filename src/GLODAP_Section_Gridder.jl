module GLODAP_Section_Gridder

root = dirname(@__FILE__)[1:end-4]
# Little bit hacky but it's probably easier to use the package root rather than 
# src.

using MATLAB
using DIVAnd
using CSV
using DelimitedFiles
using NCDatasets
using DataFrames
using Interpolations
using Base.Threads
using TOML
using Plots
# Write your package code here.
include("data_inspection.jl") # Checkec
include("data_loading.jl") # Checkec
include("simple_functionality.jl")
include("exception_handling.jl")
include("partial_cruises.jl")
include("distances_scales.jl")
include("correlation_lengths.jl")
include("DIVA_WrapperFunctions.jl")
include("pipelines.jl")
include("background_fields.jl")
include("output_fields.jl")

end
