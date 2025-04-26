using Test
include("../src/distances_scales.jl")


#=
@testset "grid_horz_dist Tests" begin
    lon_grid = [0.0, 90.0, 180.0, 270.0]
    lon = -90.0

    @test grid_vert_dist(lon_grid, lon) == 270.0

end
=#


@testset "grid_vert_dist Tests" begin
    pr_grid = [0.0, 90.0, 180.0, 270.0]

    @test grid_vert_dist(pr_grid) == [10,10,10,10]

end


@testset "create_sigma_grid Tests" begin
    # Tests with an exact fraction
    x = [v for v in 0:6000]
    gridded_sig = create_sigma_grid(x)

    expected = [v for v in 0:10:6000]

    @test gridded_sig == expected

    gridded_sig = create_sigma_grid(x,1200)
    expected = [v for v in 0:5:6000]

    @test gridded_sig == expected

    # Tests with a non-exact fraction
    x = [v for v in 0:6001]
    gridded_sig = create_sigma_grid(x)
    expected = [v for v in 0:10:6000]

    @test gridded_sig == expected

    # Test with a float
    x = [v for v in 0:0.1:59.9]
    expected = [v for v in 0:0.1:59.9]
    gridded_sig = create_sigma_grid(x)

    @test gridded_sig == expected
end


@testset "grid_sigma_distance Tests" begin
    lon_grid = [0.0, 90.0, 180.0, 270.0]
    lon = -90.0

    @test modulo_lon(lon_grid, lon) == 270.0

end


@testset "calc_scale_factors Tests" begin
    lon_grid = [0.0, 90.0, 180.0, 270.0]
    lon = -90.0

    @test modulo_lon(lon_grid, lon) == 270.0

end

