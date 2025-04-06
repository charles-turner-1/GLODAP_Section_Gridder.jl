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
using Statistics
# Write your package code here.

include("data_inspection.jl") # Checked
include("data_loading.jl") # Checked
include("simple_functionality.jl") # Checked
include("exception_handling.jl") # Checked
include("partial_cruises.jl") # Checked
include("distances_scales.jl") # Checked, still requires some work but the repackaging is fine.
include("correlation_lengths.jl") # Checked
include("DIVA_WrapperFunctions.jl") # Checked
include("pipelines.jl") # Checked
include("background_fields.jl") # Not sure if I want these functions?
include("output_fields.jl") # Checked

end
