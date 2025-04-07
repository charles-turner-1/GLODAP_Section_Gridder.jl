using Test
include("../src/simple_functionality.jl")


@testset "modulo_lon Tests" begin
    lon_grid = [0.0, 90.0, 180.0, 270.0]
    lon = -90.0

    @test modulo_lon(lon_grid, lon) == 270.0

    lon_grid = [0.0, 90.0, 180.0, 270.0]
    lon = 90.0

    @test modulo_lon(lon_grid, lon) == 90.0

    lon_grid = [-180.0, -90.0, 0.0, 90.0]
    lon = -90.0
    @test modulo_lon(lon_grid, lon) == -90.0

    lon_grid = [-180.0, -90.0, 0.0, 90.0]
    lon = 90.0
    @test modulo_lon(lon_grid, lon) == 90.0

    lon_grid = [0, 90, 180, 270]
    lon = -90
    @test modulo_lon(lon_grid, lon) == 270

    lon_grid = [0, 90, 180, 270]
    lon = 90.0
    @test modulo_lon(lon_grid, lon) == 90.0

    lon_grid = [0, 90, 180, 270]
    lon = 90
    @test modulo_lon(lon_grid, lon) == 90

    lon_grid = [0.0, 90.0, 180.0, 270.0]
    lon = [-90.0, 90]

    @test modulo_lon(lon_grid, lon) == [270.0, 90]

end

@testset "check_gridding_vars Tests" begin
    gridding = "isobaric"
    mean_val = "horzMean"

    @test isnothing(check_gridding_vars( gridding, mean_val))

    gridding = "isopycnic"
    mean_val = "scalar"

    @test isnothing(check_gridding_vars( gridding, mean_val))

    gridding = "isotonic"
    mean_val = "horzMean"

    @test_throws InvalidGriddingError check_gridding_vars( gridding, mean_val)

    gridding = "isobaric"
    mean_val = "invalid"

    @test_throws InvalidMeanValueError check_gridding_vars( gridding, mean_val)
end

@testset "to_mask_name Tests" begin
    section_name = "P01-2019"
    mask_name = to_mask_name(section_name)
    @test mask_name == "maskP01_2019"
end

@testset "remove_scalar_mean" begin
    obs_variable = [1.0, 2.0, 3.0, 4.0, NaN]

    var_mean, var_anom = remove_scalar_mean(obs_variable)

    @test isapprox(var_anom ,[-1.5, -0.5, 0.5, 1.5, NaN], nans=true)
    @test var_mean == 2.5
end

@testset "calc_decimal_year" begin
    year = [2000, 2001, 2002, 2004]
    month = [1, 1, 1, 2,]
    day = [1, 2, 3, 4]

    decimal_year = calc_decimal_year(year, month, day)

    @test isapprox(decimal_year , [2000.003, 2001.005, 2002.008, 2004.093], atol=5e-3)
end