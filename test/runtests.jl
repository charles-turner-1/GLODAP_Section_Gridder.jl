using Test
import GLODAP_Section_Gridder

const GSG = GLODAP_Section_Gridder

@testset "GLODAP_Section_Gridder.jl" begin
    include("test_simple_functionality.jl")
    include("test_distances_scales.jl")
    include("test_data_loading.jl")
end
